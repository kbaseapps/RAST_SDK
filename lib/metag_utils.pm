# -*- perl -*-

package metag_utils;

########################################################################
# This module is built to handle the annotation of an Assembly
# (KBaseGenomeAnnotations.Assembly/KBaseMetagenomes.AnnotatedMetagenomeAssembly)
# or Genome (KBaseGenomes/KBaseGenomeAnnotations.GenomeAnnotation) object
# and create a KBaseMetagenomes.AnnotatedMetagenomeAssembly object
########################################################################

use strict;
use warnings;

use Bio::KBase::kbaseenv;
use Config::IniFiles;
use Data::Dumper;
use File::Spec::Functions qw(catfile);
use File::Copy;
use Carp qw( croak );
use Bio::SeqIO;
use Bio::Perl;
use Bio::Tools::CodonTable;


use installed_clients::GenomeAnnotationAPIClient;
use installed_clients::AssemblyUtilClient;
use installed_clients::GenomeFileUtilClient;
use installed_clients::WorkspaceClient;

use gjoseqlib;


my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = Config::IniFiles->new(-file=>$config_file);
my $ws_url = $config->{'workspace-url'};
my $call_back_url = $ENV{ SDK_CALLBACK_URL };
my $rast_scratch = $config->val('RAST_SDK', 'scratch');


#-------------------------Reference from prodigal command line-------------------
#Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
#                 [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
#                 [-p mode] [-q] [-s start_file] [-t training_file] [-v]
#
#         -a:  Write protein translations to the selected file.
#         -c:  Closed ends.  Do not allow genes to run off edges.
#         -d:  Write nucleotide sequences of genes to the selected file.
#         -f:  Select output format (gbk, gff, or sco).  Default is gbk.
#         -g:  Specify a translation table to use (default 11).
#         -h:  Print help menu and exit.
#         -i:  Specify FASTA/Genbank input file (default reads from stdin).
#         -m:  Treat runs of N as masked sequence; don't build genes across them.
#         -n:  Bypass Shine-Dalgarno trainer and force a full motif scan.
#         -o:  Specify output file (default writes to stdout).
#         -p:  Select procedure (single or meta).  Default is single.
#         -q:  Run quietly (suppress normal stderr output).
#         -s:  Write all potential genes (with scores) to the selected file.
#         -t:  Write a training file (if none exists); otherwise, read and use
#              the specified training file.
#         -v:  Print version number and exit.
#--------------------------------------------------------------------------------
sub _build_prodigal_params {
    my ($fasta_file, $trans_file, $nuc_file, $output_file,
        $mode, $start_file, $training_file) = @_;

    my $prd_params = {
        input_file => $fasta_file,       # -i (FASTA/Genbank file)
        trans_fle => $trans_file,        # -a (Write protein translations to $trans_file)
        nuc_file => $nuc_file,           # -d (Write nucleotide sequences of genes to $nuc_file)
        output_type => $mode,            # -f (gbk, gff, or sco, default to gbk)
        output_file => $output_file,     # -o (Specify output file (default writes to stdout).) 
        closed_ends => 1,                # -c (Closed ends.  Do not allow genes to run off edges.)
        N_as_masked_seq => 1,            # -m (Treat runs of N as masked sequence; don't build genes across them.)
        trans_table => 11,               # -g (Specify a translation table to use (default 11).)
        procedure => 'meta',             # -p (Select procedure (single or meta).  Default is single.)
        quiet => 0,                      # -q (Run quietly)
        start_file => $start_file,       # -s (Write all potential genes (with scores) to $start_file)
        training_file => $training_file  # -t (Write a training file (if none exists); otherwise, read and use $training_file)
    };
    return $prd_params;
}

sub _run_prodigal {
    my $prodigal_params = @_;

    my @cmd = ('/kb/runtime/prodigal');
    if ($prodigal_params->{input_file}) {
        push @cmd, '-i $prodigal_params->{input_file}';
    }
    else {
        die "An input FASTA/Genbank file is required for Prodigal to run.\n";
    }
    my $out_type = 'gbk';
    if ($prodigal_params->{output_type}) {
        $out_type = $prodigal_params->{output_type};
    }
    push @cmd, '-f $out_type';
    if ($prodigal_params->{output_file}) {
        push @cmd, '-o $prodigal_params->{output_file}.' . $out_type;
    }
    else {
        push @cmd, '-o prodigal_output.' . $out_type;
    }
    if ($prodigal_params->{trans_file}) {
        push @cmd, '-a $prodigal_params->{trans_file}';
    }
    if ($prodigal_params->{nuc_file}) {
        push @cmd, '-d $prodigal_params->{nuc_file}';
    }
    if ($prodigal_params->{closed_ends} != 0) {
        push @cmd, '-c';
    }
    if ($prodigal_params->{trans_table} != 11) {
        push @cmd, '-g $prodigal_params->{trans_table}';
    }
    if ($prodigal_params->{procedure} != 'single') {
        push @cmd, '-p $prodigal_params->{procedure}';
    }
    if ($prodigal_params->{start_file}) {
        push @cmd, '-s $prodigal_params->{start_file}';
    }
    if ($prodigal_params->{training_file}) {
        push @cmd, '-t $prodigal_params->{training_file}';
    }
    if ($prodigal_params->{quiet} != 0) {
        push @cmd, '-q';
    }
    return system(@cmd);  # success if return 0
}

#--From https://github.com/hyattpd/prodigal/wiki/understanding-the-prodigal-output---
# By default, Prodigal produces one output file, which consists of gene coordinates
# and some metadata associated with each gene.
# However, the program can produce four more output files at the user's request:
#       protein translations (with the -a option),
#       nucleotide sequences (with the -d option),
#       a complete listing of all start/stop pairs along with score information (with the -s option),
#       a summary of statistical information about the genome or metagenome (with the -w option).
#
# By default, Prodigal produces a Genbank-like feature table; however, the user can
# specify some other output types via the -f option:
#       gbk:  Genbank-like format (Default)
#       gff:  GFF format
#       sqn:  Sequin feature table format
#       sco:  Simple coordinate output
#
# In the following parsing method, we assume "sco" is the output_type together
# with a translation file (protein translations) and a nucleotide file (sequences).
#
# An sco file has a head portion like:
#-----------------------------------------------------------
# # Sequence Data: seqnum=1;seqlen=4641652;seqhdr="Escherichia coli str. K-12 substr. MG1655, complete genome."
# Model Data: version=Prodigal.v2.6.3;run_type=Single;model="Ab initio";gc_cont=50.79;transl_table=11;uses_sd=1
# >1_337_2799_+
# >2_2801_3733_+
# >3_3734_5020_+
# >4_5234_5530_+
# >5_5683_6459_-
# >6_6529_7959_-
# >7_8238_9191_+
# ...
#
# ##gff-version  3
# Sequence Data: seqnum=1;seqlen=4639675;seqhdr="NC_000913 # Escherichia coli str. K-12 substr. MG1655, complete genome."
# Model Data: version=Prodigal.v3.0.0-devel.1.0;run_type=Normal;model="Ab initio";gc_cont=50.79;transl_table=11;uses_sd=1
# NC_000913	Prodigal_v3.0.0-devel.1.0	CDS	337	2799	338.7	+	0	ID=1_2;partial=00;start_type=ATG;stop_type=TGA;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.531;conf=99.99;score=338.70;cscore=322.16;sscore=16.54;rscore=11.24;uscore=1.35;tscore=3.95;
# ...
#
# Protein Translations:
# The protein translation file consists of all the proteins from all the sequences
# in multiple FASTA format. The FASTA header begins with a text id consisting of
# the first word of the original FASTA sequence header followed by an underscore
# followed by the ordinal ID of the protein. This text id is not guaranteed to be
# unique (it depends on the FASTA headers supplied by the user), which is why we
# recommend using the "ID" field in the final semicolon-delimited string instead.
#
# A translation file has an entry like;
#-----------------------------------------------------------
# >Escherichia_2 # 2801 # 3733 # 1 # ID=1_2;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;gc_cont=0.563
# MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEP
# RENIVYQCWERFCQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLND
# TRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGI
# KVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRERLLP
# GFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRL
# DTAGARVLEN*
# ...
#--------------------------------------------------------------------------------
sub _parse_prodigal_results {
    my ($trans_file, $output_file, $mode) = @_;
    my %transH;
    if ($mode == 'gff') {
        my ($gff_contents, $attr_del) = _parse_gff($output_file, '=');
        return $gff_contents;
    }
    elsif ($mode == 'sco') {
        %transH = _parse_translation($trans_file);
        my $sco_tbl = _parse_sco($output_file, %transH);
        return $sco_tbl;
    }
    return [];
}

sub _parse_translation {
    my ($trans_file) = @_;

    my %transH;
    my ($fh_trans, $trans_id, $comment, $seq);
    $fh_trans = _openRead($trans_file);

    while (($trans_id, $comment, $seq) = &gjoseqlib::read_next_fasta_seq($fh_trans)) {
        my ($contig_id) = ($trans_id =~ m/^(\S+)_\d+$/o);
        my ($left, $right, $strand, $left_trunc, $right_trunc) 
            = ($comment =~ m/^\#\s+(\d+)\s+\#\s+(\d+)\s+\#\s+(-?1)\s+\#.*partial=([01])([01])/o);
	
        if ($contig_id && $left && $right && $strand && defined($left_trunc) && defined($right_trunc)) {
            $seq =~ s/\*$//o;
            $strand = ($strand == 1) ? q(+) : q(-);
            $transH{"$contig_id\t$left\t$right\t$strand"} = [ $seq, $left_trunc, $right_trunc ];
        }
        else {
            die ("Could not parse record:\n",
                 "trans_id=$trans_id\n",
                 "comment=$comment",
                 "left=$left\n",
                 "right=$right\n",
                 "strand=$strand\n",
                 "left_trunc=$left_trunc\n",
                 "right_trunc=$right_trunc\n",
            );
        }
    }
    return %transH;
}

sub _parse_sco {
    my ($sco_file, %transH) = @_;

    my $encoded_tbl = [];
    my $fh_sco = openRead($sco_file);
    while (defined(my $line = <$fh_sco>)) {
        chomp $line;
        my $contig_id;
        if ($line =~ m/^\# Sequence Data:.*seqhdr=\"([^\"]+)\"/o) {
            $contig_id = $1;
            next;
        }
        
        if ($line =~ m/^\# Model Data/o) {
            next;
        }
        
        if (my ($num, $left, $right, $strand) = ($line =~ m/^\>(\d+)_(\d+)_(\d+)_([+-])/o)) {
            my ($beg, $end, $trunc_flag);
            
            if (my ($seq, $trunc_left, $trunc_right) = @{$transH{"$contig_id\t$left\t$right\t$strand"}}) {
                my $len = 1 + $right - $left;
            
                if ($strand eq q(+)) {
                    ($beg, $end) = ($left, $right);
                    $trunc_flag = "$trunc_left,$trunc_right";
                }
                else {
                    ($beg, $end) = ($right, $left);
                    $trunc_flag = "$trunc_right,$trunc_left";
                }
            
                push @$encoded_tbl, [$contig_id, $beg, $end, $strand, $len, $seq, $trunc_flag];
            }
            else {
                warn "No translation found for \"$sco_file\" line: $line\n";
            }
        }
        else {
            warn "Could not parse calls for \"$sco_file\" line: $line\n";
        }
    }
    close $fh_sco;
    return $encoded_tbl;
}



sub _get_fasta_from_assembly {
    my $assembly_ref = @_;

    my $au = new installed_clients::AssemblyUtilClient($call_back_url);
    my $output = $au->get_assembly_as_fasta({"ref" => $assembly_ref});
    return $output->{path};
}

sub _write_fasta_from_metagenome {
    my ($fasta_filename, $input_obj_ref) = @_;

    my $ws_client = new installed_clients::WorkspaceClient($ws_url, token => $token);
    my $genome_obj = $ws_client->get_objects([{ref=>$input_obj_ref}])->[0]->{data};
    my $fa_file = _get_fasta_from_assembly($genome_obj->{assembly_ref});
    copy($fa_file, $fasta_filename);

    unless (-s $fasta_filename) {print "Fasta file is empty.";}
    return $fasta_filename;
}

sub _write_gff_from_metagenome {
    my ($gff_filename, $genome_ref) = @_;

    my $gfu = new installed_clients::GenomeFileUtilClient($call_back_url);
    my $gff_result = $gfu.metagenome_to_gff({"genome_ref" => $genome_ref});
    copy($gff_result->{file_path}, $gff_filename);

    unless (-s $gff_filename) {print "GFF is empty ";}
    return $gff_filename;
}

#----FILE IO ----#
sub _openWrite {
    # Open a file for writing
    my ($fn) = @_;
    open my $fh, qw(>), $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub _openRead {
    # Open a file for reading
    my ($fn) = @_;
    open my $fh, qw(<), $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

##----subs for parsing GFF by Seaver----##
sub _parse_gff {
    my ($gff_filename, $attr_delimiter) = shift;

    # Open $gff_filename to read into an array
    my $fh = _openRead($gff_filename);
    my @gff_lines=();
    chomp(@gff_lines = <$fh>);
    close($fh);

    my @gff_contents=();
    foreach my $current_line (@gff_lines){
        my ($contig_id, $source_id, $feature_type, $start, $end,
            $score, $strand, $phase, $attributes) = split("\t",$current_line);

        #Some lines in a GFF can have an empty attributes column
        if(!defined($attributes) || $attributes =~ /^s\s*$/) {
            push(@gff_contents,[split("\t",$current_line)]);
            next;
        }

        # Populating with attribute key-value pair
        # This is where the feature id is from
        my %ftr_attributes=();
        my @attr_order=();
        foreach my $attribute (split(";",$attributes)) {
            chomp $attribute;

            # Sometimes empty string
            next if $attribute =~ /^\s*$/;

            # Use of 1 to limit split as '=' character can also be made available later
            # Sometimes lack of "=", assume spaces instead
            my ($key,$value)=(undef,undef);
            $attr_delimiter="=";
            if ($attribute =~ /=/) {
                ($key, $value) = split("=", $attribute, 2);
            } elsif($attribute =~ /\s/) {
                ($key, $value) = split(" ", $attribute, 2);
                $attr_delimiter=" ";
            }

            if(!defined($key)) {
                print "Warning: $attribute not parsed right\n";
            }

            #Force to lowercase in case changes in lookup
            $ftr_attributes{lc($key)}=$value;
            push(@attr_order,lc($key));
        }

        my $gff_line_contents = [$contig_id, $source_id, $feature_type, $start, $end,
                                 $score, $strand, $phase, \%ftr_attributes];
        push(@gff_contents, $gff_line_contents);
    }
    return (\@gff_contents, $attr_delimiter);
}


sub _update_gff_functions_from_features {
    my ($gff_contents, $features) = @_;

    #Feature Lookup Hash
    my %ftrs_function_lookup = ();
    foreach my $ftr (@$features){
        next if !exists($ftr->{'functions'});
        $ftrs_function_lookup{$ftr->{'id'}}=join(" / ",@{$ftr->{'functions'}});
 
        # Use these lines if the feature is an old type using singular 'function' field
        #next if !exists($ftr->{'function'});
        #$ftrs_function_lookup{$ftr->{'id'}}=$ftr->{'function'};
    }

    my @new_gff_contents=();
    foreach my $gff_line (@$gff_contents) {
        if(scalar(@$gff_line) < 9){
            #No attributes column
            push(@new_gff_contents,$gff_line);
            next;
        }

        my $ftr_attributes = $gff_line->[8];

        # According to the genome loading code, the function must be added to the produce attribute (go figure)
        # https://github.com/kbaseapps/GenomeFileUtil/blob/master/lib/GenomeFileUtil/core/FastaGFFToGenome.py#L665-L666
        # Note that this overwrites anything that was originally in the 'product' field it it previously existed
        # Also note that I'm forcing every feature to have at least an empty product field
	#
	# If the gene already has a function and RAST fails to reannotate it, keep the original function.
	if (!defined($ftr_attributes{'product'})) {
            $ftr_attributes{'product'} = "";
        }

        #Look for, and add function
        if(exists($ftrs_function_lookup{$ftr_attributes->{'id'}})) {
            $ftr_attributes->{'product'}=$ftrs_function_lookup{$ftr_attributes->{'id'}};
        }

        $gff_line->[8] = $ftr_attributes;
        push(@new_gff_contents,$gff_line);
    }

    return \@new_gff_contents;
}

sub _write_gff {
    my ($gff_contents, $gff_filename, $attr_delimiter) = @_;

    # Open $gff_filename to write the @readin_array back to the file
    my $fh = _openWrite($gff_filename);

    # Loop over the array
    foreach (@$gff_contents) {
        if(scalar(@$_)>=9){
            # Reform the attributes string
            my $attributes = $_->[8];
            my @attributes = ();
            foreach my $key ( sort keys %$attributes ) {
                my $attr=$key.$attr_delimiter.$attributes->{$key};
                push(@attributes,$attr);
            }
            $_->[8] = join(";",@attributes);
        }

	print $fh join("\t", @$_)."\n";
    }
    close $fh;
}


sub _parse_fasta {
    my ($fasta_filename) = @_;

    my @fasta_lines = ();
    # Open $fasta_filename to read into an array
    my $fh = openRead($fasta_filename);
    chomp(@fasta_lines = <$fh>);
    close($fh);

    # Read in the fasta
    my %fasta_contents = ();
    my $seqio = Bio::SeqIO->new( -file => $fasta_filename, -format => 'fasta' );
    while(my $seq_obj = $seqio->next_seq) {
        $fasta_contents{$seq_obj->id}=$seq_obj->seq;
    }
    return \%fasta_contents;
}


sub _extract_cds_sequences_from_fasta {
    my ($fasta_contents, $gff_contents) = @_;

    my %gene_seqs = ();
    foreach my $gff_line (@$gff_contents) {
        #First, there are gff "lines" that are empty
        if(scalar(@$gff_line)<9) {
            next;
        }

        #Secondly, and this is critical, this is assuming that we are looking at the CDS subtype
        #You really must check with the MAG folks about the subtypes in the GFF file
        next if $gff_line->[2] !~ /CDS/i;

        #Thirdly, you have to check the contig id
        my $contig_id = $gff_line->[0];
        if(!exists($fasta_contents->{$contig_id})) {
            print "Warning: the contig id $contig_id in the gff file doesn't exist in the fasta file!\n";
            next;
        }
        my $contig_seq = $fasta_contents->{$contig_id};

        my ($start, $end, $strand) = ($gff_line->[3],$gff_line->[4],$gff_line->[6]);

        my $gene_length = $end - $start + 1;
        my $gene_seq = substr($contig_seq, $start-1, $gene_length);

        # Before 'storing' the sequence, check strand and return reverse complement if necessary
        # Lifted from BioPerl
            # https://github.com/bioperl/bioperl-live/blob/master/lib/Bio/PrimarySeqI.pm#L421-L422
        if($strand eq '-') {
            $gene_seq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
            $gene_seq = scalar reverse $gene_seq;
        }

        $gene_seqs{$gff_line->[8]{'id'}}=$gene_seq;
    }
    return \%gene_seqs;
}

sub _translate_gene_to_protein_sequences {
    my $gene_seqs = shift;
    my $codon_table  = Bio::Tools::CodonTable -> new ( -id => 'protein' );

    my %protein_seqs = ();
    foreach my $gene (sort keys %$gene_seqs){
        my $gene_seq = $gene_seqs->{$gene};
            my $protein_seq  = $codon_table->translate($gene_seq);
        $protein_seqs{$gene}=$protein_seq;
    }
    return \%protein_seqs;
}

##----end from Seaver----##

sub rast_metagenome {
    my $params = @_;
    my $input_obj_ref = $params->{object_ref};
    my $inputgenome = {
        features => []
    };
    my $out_metag_name = $params -> {output_metagenome_name};
    my $ws = $params -> {output_workspace};

    my $gff_filename = catfile($rast_scratch, 'genome.gff');
    my $new_gff_file = catfile($rast_scratch, 'new_genome.gff');
    my $input_fasta_file = catfile($rast_scratch, 'input.fasta');
    my $trans_file = catfile($rast_scratch, 'protein_translation');
    my $nuc_file = catfile($rast_scratch, 'nucleotide_seq');
    my $output_file = catfile($rast_scratch, 'prodigal_out');
    my $start_file = catfile($rast_scratch, 'start_file');
    my $training_file = '';  # catfile($rast_scratch, 'training_file');

    my $ws_client = new installed_clients::WorkspaceClient($ws_url, token => $token);
    my $info = $ws_client->get_object_info([{ref=>$input_obj_ref}],0);

    my ($fasta_contents, $gff_contents, $attr_delimiter) = ([], [], "=");
    # Check if input is an assembly, if so run Prodigal and parse for proteins
    if ($info->[0]->[2] =~ /Assembly/) {
	my $out_file = _get_fasta_from_assembly($input_obj_ref);
	copy($out_file, $input_fasta_file) || die "Could not find file: ".$out_file;

        my $mode = 'gff';
        my $prodigal_params = _build_prodigal_params($input_fasta_file,
                                                     $trans_file,
                                                     $nuc_file,
                                                     $output_file,
                                                     $mode,
                                                     $start_file,
                                                     $training_file);

        if (_run_prodigal($prodigal_params) == 0) {
            # Prodigal finished run, files are written into $output_file/$trans_file/$nuc_file/$start_file/$training_file
            my $prodigal_result = _parse_prodigal_results($trans_file, $output_file, $mode);

            my $count = @$prodigal_result;
            my $cur_id_suffix = 1;
            foreach my $entry (@$prodigal_result) {
                # print Data::Dumper->Dump($entry)."\n";
                my ($contig, $beg, undef, $strand, $length, $translation) = @$entry;

                my $id = join(".", "peg", $cur_id_suffix);
                $cur_id_suffix++;
                push(@{$inputgenome->{features}}, {
                        id                  => $id,
                        type                => 'CDS',
                        location            => [[ $contig, $beg, $strand, $length ]],
                        annotator           => 'prodigal',
                        annotation          => 'Add feature called by PRODIGAL',
                        protein_translation => $translation
                });
            }
            $gff_filename = $output_file;  # assuming Prodigal generates a GFF file
            $fasta_contents = _parse_fasta($input_fasta_file);
            ($gff_contents, $attr_delimiter) = _parse_gff($gff_filename, $attr_delimiter);
        }
        else {
            die "Prodigal run failed.";
        }
    }
    else {# input is a (meta)genome, get its protein sequences and gene IDs
        # generating the fasta and gff files for saving the annotated metagenome
        $input_fasta_file = _write_fasta_from_metagenome($input_fasta_file, $input_obj_ref);
        $gff_filename = _write_gff_from_metagenome($gff_filename, $input_obj_ref);

        # 2) fetch protein sequences and gene IDs from the above fasta and gff files
        my $fasta_contents = _parse_fasta($input_fasta_file);
        ($gff_contents, $attr_delimiter) = _parse_gff($gff_filename, $attr_delimiter);

        my $gene_seqs = _extract_cds_sequences_from_fasta($fasta_contents, $gff_contents);
        my $protein_seqs = _translate_gene_to_protein_sequences($gene_seqs);

        my %gene_id_index=();
        my $i=1;
        foreach my $gene (sort keys %$protein_seqs){
            push(@{$inputgenome->{features}},{
                id => "peg".$i,
                protein_translation => $protein_seqs->{$gene}
            });
            $gene_id_index{$i}=$gene;
            $i++;
		}
    }
    # Call RAST to annotate the proteins/genome
    my $rast_client = Bio::kbase::kbaseenv::ga_client();
    my $rasted_genome = $rast_client->run_pipeline($inputgenome,
            {stages => [{name => "annotate_proteins_kmer_v2", kmer_v2_parameters => {}},
                        {name => "annotate_proteins_similarity",
                         similarity_parameters => { annotate_hypothetical_only => 1 }}]}
    );

    my $ftrs = $rasted_genome->{features};
    my $updated_gff_contents = _update_gff_functions_from_features($gff_contents, $ftrs);

    # call $gfu->fasta_gff_to_metagenome() to save the annotated (meta)genome
    _write_gff($updated_gff_contents, $new_gff_file, $attr_delimiter);

    my $gfu = new installed_clients::GenomeFileUtilClient($call_back_url);
    my $annotated_metag = $gfu->fasta_gff_to_metagenome ({
            "fasta_file" => {'path' => $input_fasta_file},
            "gff_file" => {'path' => $new_gff_file},
            "genome_name" => $out_metag_name,
            "workspace_name" => $ws,
            "generate_missing_genes" => 1
    });

    return $annotated_metag->{genome_ref};
}

