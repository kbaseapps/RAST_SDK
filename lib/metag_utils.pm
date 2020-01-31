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
my $ws_client = new installed_clients::WorkspaceClient($ws_url, token => $token);
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
sub build_prodigal_params {
    my ($ref, $fasta_file, $trans_file, $nuc_file, $output_file, $start_file, $training_file) = @_;

    my $out_file = get_fasta_from_assembly($ref);
    copy($out_file, $fasta_file) || die "Could not find file: ".$out_file;

    my $prd_params = {
        input_file => $fasta_file,       # -i (FASTA/Genbank file)
        trans_fle => $trans_file,        # -a (Write protein translations to $trans_file)
        nuc_file => $nuc_file,           # -d (Write nucleotide sequences of genes to $nuc_file)
        output_type => 'sco',            # -f (gbk, gff, or sco, default to gbk)
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

sub run_prodigal {
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
sub parse_prodigal_results {
    my ($trans_file, $nuc_file, $output_file) = @_;

    my %transH;
    my ($fh_trans, $trans_id, $comment, $seq);
    open($fh_trans, q(<), $trans_file) or die qq(Could not read-open \"$trans_file\");

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
    
    my $fh_sco;
    my $encoded_tbl = [];
    open($fh_sco, q(<), $output_file) or die qq(Could not read-open sco_file=\"$output_file\");
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
            
            if (my ($seq, $trunc_left, $trunc_right) = @ { $transH{"$contig_id\t$left\t$right\t$strand"} })
            {
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
                warn "No translation found for \"$output_file\" line: $line\n";
            }
        }
        else {
            warn "Could not parse calls for \"$output_file\" line: $line\n";
        }
    }
    close($fh_sco);
    
    return $encoded_tbl;
}

sub get_fasta_from_assembly {
    my $assembly_ref = @_;

    my $au = new installed_clients::AssemblyUtilClient($call_back_url);
    my $output = $au->get_assembly_as_fasta({"ref" => $assembly_ref});
    return $output->{path};
}

sub write_fasta_from_metagenome {
    my ($fasta_filename, $input_obj_ref) = @_;

    my $genome_obj = $ws_client->get_objects([{ref=>$input_obj_ref}])->[0]->{data};
    my $fa_file = get_fasta_from_assembly($genome_obj->{assembly_ref});
    copy($fa_file, $fasta_filename);

    unless (-s $fasta_filename) { print "Fasta file is empty.";}
    return $fasta_filename;
}

sub write_gff_from_metagenome {
    my ($gff_filename, $genome_ref) = @_;

    my $gfu = new installed_clients::GenomeFileUtilClient($call_back_url);
    my $gff_result = $gfu.metagenome_to_gff({"genome_ref" => $genome_ref});
    move($gff_result->{file_path}, $gff_filename);

    unless (-s $gff_filename) {print "GFF is empty ";}
    return $gff_filename;
}

##------ subs for converting fasta and gff files into protein sequences
sub revcompl {
    my ($seq) = @_;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return scalar reverse $seq;
}

sub fasta_cut {
    #-----
    # Cut up a fasta sequence
    #
    my ($fa_str, $prot, $line_wrap) =  @_;

    # translate if need be
    if(0 != $prot)
    {
        my $codon_table  = Bio::Tools::CodonTable -> new ( -id => $prot );
        $fa_str = $codon_table->translate($fa_str);
    }

    # wrap the line if need be
    if(0 != $line_wrap)
    {
        my $return_str = "";
        my $len = length $fa_str;
        my $start = 0;
        while($start < $len)
        {
            $return_str .= substr $fa_str, $start, $line_wrap;
            $start += $line_wrap;
        }
        return $return_str;
    }
    return "$fa_str";
}

sub parse_proteins_from_gff_fasta {
    my ($f_gff, $f_fasta, $f_out, $ftr) = @_;
    my %gff_orfs = ();
    my $reject_length = 50;
    my $keep_non_orfs = 0;
    my $protein_code = 0;
    my $include_nulls = 0;
    my $seq_desc = 0;
    if((0 != $keep_non_orfs) and (0 != $protein_code))
    {
        print "**WARNING: $0 : non_orfs flag is invalid when translating into protein space --> ignoring\n";
        $keep_non_orfs = 0;
    }
    my $line_wrap = 80;

    my @proteins = ();

    # Open $f_gff to read in the gff into %gff_orfs
    my $gff_fh;
    unless (open( $gff_fh, q(<), $f_gff )) {
        croak "Could not open file $f_gff $!";

    while(<$gff_fh>){
        next if ($_ =~ m/^#/);
        my ($seqid, undef, $feature, $start, $end,
            undef, $strand, undef, $attributes) = split;
        ($attributes) = ($attributes =~ m/(ID|Alias)=([^;]*)/i)[1];
        push @{$gff_orfs{$seqid}}, [$start, $end, $strand, $attributes, $feature];
    }
    close $gff_fh;

    # open the output file
    my $out_fh;
    unless (open( $out_fh, q(>), $f_out )) {
        croak "Could not open file $f_out $!";
    }

    # Read in the fasta
    my $seqio = Bio::SeqIO->new( -file => $f_fasta, -format => 'fasta' ) or croak "**ERROR: Could not open FASTA file: $f_fasta $!\n";
    while(my $sobj = $seqio->next_seq)
    {
        my $seqid = $sobj->id;
        my $seq = $sobj->seq;
        my $seq_length = $sobj->length;
        my $pro_entry;

        # make sure this guy has an annotation
        if(defined($gff_orfs{$seqid})) {
            for(@{$gff_orfs{$seqid}}) {
                my ($start, $end, $strand, $attributes, $feature) = @$_;

                next if (defined $ftr && $ftr !~ m/$feature/);

                # work out the length of the sub string
                my $length = $end - $start + 1;

                # check if he's long enough
                if($length < $reject_length)
                {
                    print "Rejecting: $seqid -> ($start, $end) on $strand. $length is shorter than cutoff!\n";
                    next;
                }

                my $this_seq = substr($seq, $start-1, $length);
                if($strand eq "+")
                {
                    if ($seq_desc) {
                        $pro_entry = ">$attributes | $seqid:$start-$end FORWARD";
                    }
                    else {
                        $pro_entry = ">$seqid"."_$start"."_$end"."_F";
                    }
                    push @proteins, $pro_entry;
                    print $out_fh $pro_entry."\n";
                    $pro_entry = fasta_cut($this_seq, $protein_code, $line_wrap);
                    push @proteins, $pro_entry;
                    print $out_fh $pro_entry."\n";
                }
                else {
                    if ($seq_desc) {
                        $pro_entry = ">$attributes | $seqid:$start-$end REVERSE";
                    }
                    else {
                        $pro_entry = ">$seqid"."_$start"."_$end"."_R";
                    }
                    push @proteins, $pro_entry;
                    print $out_fh $pro_entry."\n";
                    $pro_entry = fasta_cut(revcompl($this_seq), $protein_code, $line_wrap);
                    push @proteins, $pro_entry;
                    print $out_fh $pro_entry."\n";
                }
            }
        }
        elsif($include_nulls)
        {
            # include anyway
            $pro_entry =">$seqid"."_1_$seq_length"."_X\n".fasta_cut($seq, $protein_code, $line_wrap);
            push @proteins, $pro_entry;
            print $out_fh $pro_entry."\n";
        }
    }
    close $out_fh;
    return @proteins;
}


sub add_functions_to_gff {
    my ($gff_filename, $ftrs) = @_;

    my $fh;
    # Open $gff_filename to read into an array
    my @readin_arr;
    unless (open( $fh, q(<), $gff_filename )) {
        croak "Could not open file '$gff_filename' $!";
    }
    chomp(@readin_arr = <$fh>);
    close($fh);

    ## -- by Seaver: insert the functions into the array @readout_arr
    # Feature Lookup Hash
    my %ftrs_function_lookup = ();
    foreach my $ftr (@$ftrs){
        next if !exists($ftr->{'functions'});
        $ftrs_function_lookup{$ftr->{'id'}}=join(" / ",@{$ftr->{'functions'}});

        # Use these lines if the feature is an old type using singular 'function' field
        #next if !exists($ftr->{'function'});
        #$ftrs_function_lookup{$ftr->{'id'}}=$ftr->{'function'};
    }

    my @readout_arr = ();
    foreach my $current_line (@readin_arr){
        my ($contig_id, $source_id, $feature_type, $start, $end,
	        $score, $strand, $phase, $attributes) = split("\t",$current_line);

        # Some lines in a GFF can be completely empty
        if(!defined($attributes) || $attributes =~ /^s\s*$/){
	        push(@readout_arr, $current_line);
	        next;
        }

        # Populating with attribute key-value pair
        # This is where the feature id is from
        my %ftr_attributes=();
        my @attr_order=();
        my $delimiter="=";

        foreach my $attribute (split(";",$attributes)){
	        chomp $attribute;

            # Sometimes empty string
	        next if $attribute =~ /^\s*$/;

	        # Use of 1 to limit split as '=' character can also be made available later
	        # Sometimes lack of "=", assume spaces instead
	        my ($key,$value)=(undef,undef);
	        $delimiter="=";
	        if($attribute =~ /=/){
	            ($key, $value) = split("=", $attribute, 2);
	        }elsif($attribute =~ /\s/){
	            ($key, $value) = split(" ", $attribute, 2);
	            $delimiter=" ";
	        }

	        if(!defined($key)){
	            print "Warning: $attribute not parsed right\n";
	        }

	        #Force to lowercase in case changes in lookup
	        $ftr_attributes{lc($key)}=$value;
	        push(@attr_order,lc($key));
        }

        # According to the genome loading code, the function must be added to the product attribute (go figure)
        # https://github.com/kbaseapps/GenomeFileUtil/blob/master/lib/GenomeFileUtil/core/FastaGFFToGenome.py#L665-L666
        if(!exists($ftr_attributes{'product'})){
	        push(@attr_order,'product');
        }
 
        # Note that this overwrites anything that was originally in the 'product' field it it previously existed
        # Also note that I'm forcing every feature to have at least an empty product field
        if (!defined($ftr_attributes{'product'})) {
            $ftr_attributes{'product'}="";
        }

        #Look for, and add function
        if(exists($ftrs_function_lookup{$ftr_attributes{'id'}})){
	        $ftr_attributes{'product'}=$ftrs_function_lookup{$ftr_attributes{'id'}};
        }

        # Reform the attributes string
        my @new_attributes=();
        foreach my $attr (@attr_order){
	        $attr.=$delimiter.$ftr_attributes{$attr};
	        push(@new_attributes,$attr);
        }
        my $new_attributes=join(";",@new_attributes);
        my $new_line = join("\t",($contig_id, $source_id, $feature_type, $start, $end,
			                  $score, $strand, $phase, $new_attributes));
        push(@readout_arr,$new_line);
    }

    # Open a new file to write the @readout_arr back to the file
    my $new_gff_filename = catfile($rast_scratch, 'new_genome.gff');
    open( $fh, q(>), $new_gff_filename ) || die "Could not open file '$new_gff_filename' $!";
    # Loop over the array
    foreach (@readout_arr)
    {
        print $fh "$_\n";
    }
    close $fh;

    return $new_gff_filename;
}


sub rast_metagenome {
    my $params = @_;
    my $input_obj_ref = $params->{object_ref};
    my $inputgenome = {
        features => []
    };
    my $out_metag_name = $params -> {output_metagenome_name};
    my $ws = $params -> {output_workspace};

    my $gn_gff_file = catfile($rast_scratch, 'genome.gff');
    my $input_fasta_file = catfile($rast_scratch, 'input.fasta');
    my $trans_file = catfile($rast_scratch, 'protein_translation');
    my $nuc_file = catfile($rast_scratch, 'nucleotide_seq');
    my $output_file = catfile($rast_scratch, 'prodigal_out');
    my $start_file = catfile($rast_scratch, 'start_file');
    my $training_file = '';  # catfile($rast_scratch, 'training_file');

    my $info = $ws_client->get_object_info([{ref=>$input_obj_ref}],0);

    # Check if input is an assembly, if so run Prodigal and parse for proteins
    if ($info->[0]->[2] =~ /Assembly/) {
        my $prodigal_params = build_prodigal_params($input_obj_ref,
                                                    $input_fasta_file,
                                                    $trans_file,
                                                    $nuc_file,
                                                    $output_file,
                                                    $start_file,
                                                    $training_file);

        if (run_prodigal($prodigal_params) == 0) {
            # Prodigal finished run, files are written into $output_file/$trans_file/$nuc_file/$start_file/$training_file
            my $prodigal_result = parse_prodigal_results($trans_file, $nuc_file, $output_file);

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
        }
        else {
            die "Prodigal run failed.";
        }
    }
    else {# input is a (meta)genome, get its protein sequences and gene IDs

        # 1) generating the fasta and gff files for saving the annotated metagenome
        $input_fasta_file = write_fasta_from_metagenome($input_fasta_file, $input_obj_ref);
        $gn_gff_file = write_gff_from_metagenome($gn_gff_file, $input_obj_ref);

        # 2) fetch protein sequences and gene IDs from the above fasta and gff files
        my @proteins = parse_proteins_from_gff_fasta($gn_gff_file,
                                                     $input_fasta_file,
                                                     $trans_file);
        my $i = 1;
        foreach my $prot (@proteins) {
            push(@{$inputgenome->{features}},{
                id => "peg".$i,
                protein_translation => $prot
            });
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
    my $new_gff_file = add_functions_to_gff($gn_gff_file, $ftrs);

    # call $gfu->fasta_gff_to_metagenome() to save the annotated (meta)genome
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

