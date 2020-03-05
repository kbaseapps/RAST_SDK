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
use Config::Simple;
use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catfile splitpath);
use File::Path qw(make_path);
use File::Copy;
use Carp qw(croak);
use Bio::SeqIO;
use Bio::Perl;
use Bio::Tools::CodonTable;
use JSON;
use Encode qw(encode decode);
use File::Basename;

use Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl;
use installed_clients::GenomeAnnotationAPIClient;
use installed_clients::AssemblyUtilClient;
use installed_clients::GenomeFileUtilClient;
use installed_clients::WorkspaceClient;
use installed_clients::KBaseReportClient;

require 'gjoseqlib.pm';


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
sub _build_prodigal_cmd {
    my ($self, $input_file, $trans_file, $nuc_file, $output_file, $output_type,
        $mode, $start_file, $training_file) = @_;

    if (!defined($input_file)) {
       croak "An input FASTA/Genbank file is required for Prodigal to run.\n";
    }
    # setting defaults
    if (!defined($output_type) || $output_type eq '') {
        $output_type = 'gff';
    }
    if (!defined($mode) || $mode eq '') {
        $mode = 'meta';
    }
    my($vol, $f_path, $file) = splitpath($input_file);
    if (!defined($output_file) || $output_file eq '') {
        $output_file = catfile($f_path, 'prodigal_output.'.$output_type);
    }
    if (!defined($trans_file) || $trans_file eq '') {
       $trans_file = catfile($f_path, 'protein_translation');
    }
    if (!defined($nuc_file) || $nuc_file eq '') {
       $nuc_file = catfile($f_path, 'nucleotide_seq');
    }

    # building the Prodigal command
    my @cmd = ('/kb/runtime/bin/prodigal');
    push @cmd, '-i'; push @cmd, "$input_file";
    push @cmd, '-f'; push @cmd, "$output_type";
    push @cmd, '-o'; push @cmd, "$output_file";
    push @cmd, '-a'; push @cmd, "$trans_file";
    push @cmd, '-d'; push @cmd, "$nuc_file";
    push @cmd, '-c';
    push @cmd, '-g'; push @cmd, 11;
    push @cmd, '-q';
    push @cmd, '-p'; push @cmd, "$mode";
    push @cmd, '-m';
    if ($start_file) {
        push @cmd, '-s';
        push @cmd, "$start_file";
    }
    if ($training_file) {
        push @cmd, '-t';
        push @cmd, "$training_file";
    }
    return @cmd;
}


sub _run_prodigal {
    my ($self, @cmd) = @_;
    my $ret = 0;
    eval {
        $ret = system(@cmd);
    };
    if ($@) {
        print Dumper(\@cmd);
        croak "ERROR Prodigal run failed: ".$@."\n";
    }
    print "Prodigal returns: $ret\n";
    return $ret;
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
    my ($self, $trans_file, $output_file, $output_type) = @_;
    my %transH = $self->_parse_translation($trans_file);
    if ($output_type eq 'gff') {
        my ($gff_contents, $attr_del) = $self->_parse_gff($output_file, '=');
        return ($gff_contents, %transH);
    }
    elsif ($output_type eq 'sco') {
        my $sco_tbl = $self->_parse_sco($output_file, %transH);
        return ($sco_tbl, %transH);
    }
    return ([], %transH);
}

sub _parse_translation {
    my ($self, $trans_file) = @_;

    my %transH;
    my ($fh_trans, $trans_id, $comment, $seq);
    $fh_trans = $self->_openRead($trans_file);

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
            croak ("Could not parse record:\n",
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
    # print "Translation table:\n".Dumper(\%transH);
    return %transH;
}

sub _parse_sco {
    my ($self, $sco_file, %transH) = @_;

    my $encoded_tbl = [];
    my $fh_sco = $self->_openRead($sco_file);
    my $contig_id;
    while (defined(my $line = <$fh_sco>)) {
        chomp $line;
        if ($line =~ m/^\# Sequence Data:.*seqhdr=\"([^\"]+)\"/o) {
            my @words = split / /, $1;
            $contig_id = $words[0];
            next;
        }
        
        if ($line =~ m/^\# Model Data/o) {
            next;
        }
        if ($contig_id) {
            if (my ($num, $left, $right, $strand) = ($line =~ m/^\>(\d+)_(\d+)_(\d+)_([+-])/o)) {
                my ($beg, $end, $trunc_flag);
                if ($transH{"$contig_id\t$left\t$right\t$strand"}) {
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
                        # warn "No translation found for \"$sco_file\" line: $line\n";
                    }
                }
                else {
                    # warn "No key \"$contig_id\t$left\t$right\t$strand\" found for \"$sco_file\" line: $line\n";
                }
            }
            else {
                warn "Could not parse calls for \"$sco_file\" line: $line\n";
            }
        }
    }
    close $fh_sco;
    return $encoded_tbl;
}

##----end for prodigal parsing----##

sub _get_fasta_from_assembly {
    my ($self, $assembly_ref) = @_;

    my $au = new installed_clients::AssemblyUtilClient($self->{call_back_url});
    my $output = {};
    eval {
        $output = $au->get_assembly_as_fasta({"ref" => $assembly_ref});
    };
    if ($@) {
        croak "ERROR calling AssemblyUtil.get_assembly_as_fasta: ".$@."\n";
    }
    else {
        return $output->{path};
    }
}

sub _write_fasta_from_metagenome {
    my ($self, $fasta_filename, $input_obj_ref) = @_;

    eval {
        my $genome_obj = $self->_fetch_object_data($input_obj_ref);

        my $fa_file = $self->_get_fasta_from_assembly(
                          $input_obj_ref.";".$genome_obj->{assembly_ref});

        copy($fa_file, $fasta_filename);

        unless (-s $fasta_filename) {print "Fasta file is empty.";}
    };
    if ($@) {
        croak "**_write_fasta_from_metagenome ERROR: ".$@."\n";
    }
    return $fasta_filename;
}

sub _write_gff_from_metagenome {
    my ($self, $gff_filename, $genome_ref) = @_;

    my $obj_info = $self->_fetch_object_info($genome_ref);
    my $in_type = $obj_info->[2];
    my $is_meta_assembly = $in_type =~ /Metagenomes.AnnotatedMetagenomeAssembly/;

    unless ($is_meta_assembly) {
        croak "ValueError: Object is not an AnnotatedMetagenomeAssembly, GFU will throw an error.\n";
    }

    my $gfu = new installed_clients::GenomeFileUtilClient($self->{call_back_url});
    my $gff_result = '';
    eval {
        $gff_result = $gfu->metagenome_to_gff({"metagenome_ref" => $genome_ref});

        copy($gff_result->{file_path}, $gff_filename);

        unless (-s $gff_filename) {print "GFF is empty ";}
    };
    if ($@) {
        croak "**_write_gff_from_metagenome ERROR: ".$@."\n";
    }
    return $gff_filename;
}

sub _save_metagenome {
    my ($self, $ws, $out_metag_name, $obj_ref, $gff_file) = @_;

    my $req_params = "Missing required parameters for saving metagenome.\n";
    my $req1 = "Both 'output_workspace' and 'output_metagenome_name' are required.\n";
    my $req2 = "Both 'obj_ref' and 'gff_file' are required.\n";
    unless (defined($ws) && defined($obj_ref) &&
            defined($out_metag_name) && defined($gff_file)) {
        croak $req_params;
    }

    print "Parameters for saving metegenome-----------------\n";
    print "Workspace name: $_[1]\n";
    print "Metagenome name: $_[2]\n";
    print "Object ref: $_[3]\n";
    print "GFF file: $_[4]\n";

    unless (defined($out_metag_name) && defined($ws)) {
        croak "**In _save_metagenome: $req1";
    }
    unless (defined($obj_ref) && defined($gff_file)) {
        croak "**In _save_metagenome: $req2";
    }
    unless (-e $gff_file) {
        croak "**In _save_metagenome: GFF file not found.\n";
    }
    unless (-s $gff_file) {
        croak "**In _save_metagenome: GFF file is empty.\n";
    }

    eval {
        $self->_fetch_object_info($obj_ref);
    };
    if ($@) {
        croak("**In _save_metagenome ERROR trying to access the input object:\n"
               .$@."\n");
    }

    print "First few 10 lines of the GFF file before call to GFU.ws_obj_gff_to_metagenome-----------\n";
    $self->_print_fasta_gff(0, 10, $gff_file);

    my $gfu = new installed_clients::GenomeFileUtilClient($self->{call_back_url});
    my $annotated_metag = {};
    eval {
        $annotated_metag = $gfu->ws_obj_gff_to_metagenome ({
            "ws_ref" => $obj_ref,
            "gff_file" => {'path' => $gff_file},
            "genome_name" => $out_metag_name,
            "workspace_name" => $ws,
            "generate_missing_genes" => 1});
    };
    if ($@) {
        croak("**In _save_metagenome ERROR calling GenomeFileUtil.ws_obj_gff_to_metagenome:\n"
               .$@."\n");
    }
    else {
        return $annotated_metag;
    }
}

sub _check_annotation_params {
    my ($self, $params) = @_;

    my $missing_params = "Missing required parameters for annotating metagenome.\n";
    unless (defined($params)) {
        print "params is not defined!!!!\n";
        croak $missing_params;
    }
    # print out the content of hash reference
    print "Checking parameters:\n". Dumper($params). "\n";

    if (!keys %$params) {
        print "params is empty!!!!\n";
        croak $missing_params;
    }
    my $req1 = "'output_workspace' is required for running rast_metagenome.\n";
    my $invald1 = "Invalid workspace name:";
    if (!defined($params->{output_workspace}) || $params->{output_workspace} eq '') {
        croak $req1;
    }
    elsif ($params->{output_workspace} !~ m/[^\\w:._-]/) {
        croak $invald1.$params->{output_workspace}.'\n';
    }

    my $req2 = "'object_ref' is required for running rast_metagenome.\n";
    my $invald2 = "Invalid workspace object reference:";
    if (!defined($params->{object_ref}) || $params->{object_ref} eq '') {
        croak $req2;
    }
    elsif ($params->{object_ref} !~ m/[^\\w\\|._-]/) {
        croak $invald2 .$params->{object_ref}.'\n';
    }
    if (!defined($params->{output_metagenome_name})
        || $params->{output_metagenome_name} eq '') {
        $params->{output_metagenome_name} = "rast_annotated_metagenome";
    }
    if (!defined($params->{run_prodigal})
        || $params->{run_prodigal} eq '') {
        $params->{run_prodiagl} = 0;
    }
    return $params;
}


# Call RAST to annotate the proteins/genome
sub _run_rast {
    my ($self, $inputgenome) = @_;
    my $count = scalar @{$inputgenome->{features}};
    print "******Run RAST pipeline on genome with $count features.******\n";
    print "For example, first 3 features: \n".Dumper(@{$inputgenome->{features}}[0..2]);

    my $rasted_gn = {};
    eval {
        my $rast_client = Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl->new();
        $rasted_gn = $rast_client->run_pipeline($inputgenome,
            {stages => [{name => "annotate_proteins_kmer_v2",
                         kmer_v2_parameters => {min_hits => "5",
                                                annotate_hypothetical_only => 0}},
                        {name => "annotate_proteins_kmer_v1",
                         kmer_v1_parameters => {dataset_name => "Release70",
                                                annotate_hypothetical_only => 0}}]}
        );
    };
    if ($@) {
        croak "ERROR calling GenomeAnnotation::GenomeAnnotationImpl->run_pipeline: ".$@."\n";

    }
    return $rasted_gn;
};


sub _generate_stats_from_ama {
    my ($self, $gn_ref) = @_;

    my $gn_info = $self->_fetch_object_info($gn_ref);
    my $is_assembly = $gn_info->[2] =~ m/GenomeAnnotations.Assembly/;
    my $is_meta_assembly = $gn_info->[2] =~ m/Metagenomes.AnnotatedMetagenomeAssembly/;

    my %gn_stats = ();
    $gn_stats{workspace} = $gn_info->[7];
    $gn_stats{id} = $gn_info->[1];
    $gn_stats{genome_type} = $gn_info->[2];

    if ($is_assembly ) {
        if (defined($gn_info->[10])) {
            $gn_stats{gc_content} = $gn_info->[10]->{'GC content'};
            $gn_stats{contig_count} = $gn_info->[10]->{'N Contigs'};
        }
        else {
            $gn_stats{gc_content} = undef;
            $gn_stats{contig_count} = undef;
        }
        $gn_stats{num_features} = 0;
        $gn_stats{feature_counts} = {};
    }
    elsif ($is_meta_assembly) {
        my $gn_data = $self->_fetch_object_data($gn_ref);
        $gn_stats{gc_content} = $gn_data->{gc_content};
        $gn_stats{contig_count} = $gn_data->{num_contigs};
        $gn_stats{num_features} = $gn_data->{num_features};
        $gn_stats{feature_counts} = $gn_data->{feature_counts};
    }
    return %gn_stats;
}

sub _generate_stats_from_gffContents {
    my ($self, $gff_contents) = @_;

    my $gff_count = scalar @{$gff_contents};
    print "INFO: Gathering stats from $gff_count GFFs------\n";

    my %gff_stats = ();
    foreach my $gff_line (@$gff_contents) {
        if(scalar(@$gff_line) < 9){
            #No attributes column
            next;
        }

        my $ftr_attrs = $gff_line->[8];
        my $ftr_attributes;
        if (ref $ftr_attrs ne 'HASH') {
            # assuming it is a string like:
            # 'id=1_2;product=3-oxoacyl-[acyl-carrier protein] reductase (EC 1.1.1.100)'

            # Populating with attribute key-value pair
            # This is where the feature id is from
            foreach my $attr (split(";",$ftr_attrs)) {
                chomp $attr;

                # Sometimes empty string
                next if $attr =~ /^\s*$/;

                # Use of 1 to limit split as '=' character can also be made available later
                # Sometimes lack of "=", assume spaces instead
                my ($key,$value)=(undef,undef);
                if ($attr =~ /=/) {
                    ($key, $value) = split("=", $attr, 2);
                } elsif($attr =~ /\s/) {
                    ($key, $value) = split(" ", $attr, 2);
                }

                if(!defined($key)) {
                    print "Warning: $attr not parsed right\n";
                }

                #Force to lowercase in case changes in lookup
                $ftr_attributes->{lc($key)}=$value;
            }
        }
        else {
            $ftr_attributes = $gff_line->[8];
        }

        if (!defined($ftr_attributes->{'product'})
                || $ftr_attributes->{'product'} eq '') {
            next;
        }
        #Look for function
        my $gene_id = $ftr_attributes->{id};
        my $frole = $ftr_attributes->{product};
        if (!exists($gff_stats{function_roles}{$frole})) {
            $gff_stats{function_roles}{$frole}{gene_count} = 1;
            $gff_stats{function_roles}{$frole}{gene_list} = $gene_id;
        }
        else {
            $gff_stats{function_roles}{$frole}{gene_count} += 1;
            $gff_stats{function_roles}{$frole}{gene_list} = join(";",
                    $gene_id, $gff_stats{function_roles}{$frole}{gene_list});
        }
        $gff_stats{gene_role_map}{$gene_id} = $frole;
    }
    #print "Stats from GFF contents--------\n".Dumper(\%gff_stats);
    return %gff_stats;
}

sub _fetch_subsystem_info {
    my $self = shift;
    my $subsys_json_file = shift;

    print "INFO: Reading subsystem info from json file------\n";

    $subsys_json_file = './templates_and_data/subsystems_March2018.json' unless defined($subsys_json_file);
    my $dirname = dirname(__FILE__);
    my $subsystem = catfile($dirname, $subsys_json_file);
    my $fh = $self->_openRead($subsystem);
    my $json_text = do { local $/; <$fh> };
    close $fh;

    my $json_data = decode_json($json_text);
    my $subsys_data = $json_data->{response}->{docs};
    my $subsys_count = scalar @{$subsys_data};
    #print $subsys_count;  # 920

    my %subsys_info = ();
    for (my $i=0; $i<$subsys_count; $i++) {
        my $subsys_id = $subsys_data->[$i]->{subsystem_id};
        $subsys_info{$subsys_id}{subsys_class} = $subsys_data->[$i]->{class};
        $subsys_info{$subsys_id}{subclass} = $subsys_data->[$i]->{subclass};
        $subsys_info{$subsys_id}{superclass} = $subsys_data->[$i]->{superclass};
        $subsys_info{$subsys_id}{role_names} = $subsys_data->[$i]->{role_name};
        $subsys_info{$subsys_id}{role_ids} = $subsys_data->[$i]->{role_id};
        $subsys_info{$subsys_id}{subsys_name} = $subsys_data->[$i]->{subsystem_name};
        $subsys_info{$subsys_id}{description} = $subsys_data->[$i]->{description};
        $subsys_info{$subsys_id}{notes} = $subsys_data->[$i]->{notes};
    }

    return %subsys_info;
}

sub _find_function_source {
    my $self = shift;
    my $func_tab_ref = shift;
    my $func = shift;

    my $func_src = 'N/A';
    foreach my $ftr_id (sort keys %$func_tab_ref) {
        if (exists($func_tab_ref->{$ftr_id}->{'functions'})) {
            my $funcs = $func_tab_ref->{$ftr_id}->{'functions'};
            if ($funcs =~ /\Q$func\E/) {
                $func_src = $func_tab_ref->{$ftr_id}->{'annotation_src'};
                last;
	    }
        }
    }
    return $func_src;
}

sub _write_html_from_stats {
    my $self = shift;
    my $obj_stats_ref = shift;
    my $gff_stats_ref = shift;
    my $subsys_ref = shift;
    my $func_tab_ref = shift;
    my $template_file = shift;
    if (ref $obj_stats_ref ne "HASH" || ref $gff_stats_ref ne "HASH"
             || ref $subsys_ref ne "HASH" || ref $func_tab_ref ne "HASH") {
        croak "You need to pass in hash references!\n";
    }

    # dereference the hashes
    my %obj_stats = %{ $obj_stats_ref };
    my %gff_stats = %{ $gff_stats_ref };
    my %subsys_info = %{ $subsys_ref };
    $template_file = './templates_and_data/table_report.html' unless defined($template_file);

    my $dirname = dirname(__FILE__);
    my $template = catfile($dirname, $template_file);

    my $fh1 = $self->_openRead($template);
    read $fh1, my $file_content, -s $fh1; # read the whole file into a string
    close $fh1;

    my $report_title = "Feature function report for genome <font color=red>$obj_stats{id}</font>";
    my $rpt_header = "<h3>$report_title:</h3>";
    my $rpt_data = ("data.addColumn('string', 'function role');\n".
                    "data.addColumn('string', 'annotation source');\n".
                    "data.addColumn('number', 'gene count');\n".
                    "data.addColumn('string', 'subsystem name');\n".
                    "data.addColumn('string', 'subsystem class');\n".
                    "data.addRows([\n");

    my $roles = $gff_stats{function_roles};
    my $genes = $gff_stats{gene_role_map};
    foreach my $role_k (sort keys %$roles) {
        # add escape to preserve double quote (")
        (my $new_role_k = $role_k) =~ s/"/\\"/g;
        $rpt_data .= '["<span style=\"white-space:nowrap;\">'."$new_role_k</span>\",";
        my $ann_src = $self->_find_function_source($func_tab_ref, $role_k);
        $rpt_data .= "\"$ann_src\",";
        $rpt_data .= "$roles->{$role_k}->{gene_count},";
        #$rpt_data .= "\"$roles->{$role_k}->{gene_list}\"],\n";

        # search for $role_k in all the $subsys_info items's role_names arrays
        my $subsys_str = '';
        my $class_str = '';
        foreach my $subsys_k (sort keys %subsys_info) {
            #print $subsys_k."---------\n";
            #print Dumper($subsys_info{$subsys_k});
            my $subsys_roles = $subsys_info{$subsys_k}{role_names};
            if ( grep {$_ eq $role_k} @$subsys_roles ) {
                #print "Found $role_k in subsystem $subsys_k\n";
                $subsys_str .= '<br>' unless $subsys_str eq '';
                $class_str .= '<br>' unless $class_str eq '';
                $subsys_str .= $subsys_info{$subsys_k}{subsys_name};
                $class_str .= $subsys_info{$subsys_k}{subsys_class};
            }
        }
        $rpt_data .= "\"$subsys_str\",";
        $rpt_data .= "\"$class_str\"],\n";
    }
    chomp $rpt_data;
    chop $rpt_data;
    $rpt_data .= "\n]);\n";

    my $gene_count = keys %$genes;
    my $rpt_footer = "<p><strong>Contig Count = $obj_stats{contig_count}</strong></p>\n";
    $rpt_footer .= "<p><strong>Feature Count = $gff_stats{num_features}</strong></p>\n";
    $rpt_footer .= "<p><strong>(Functional)Gene Total Count = $gene_count</strong></p>\n";
    $rpt_footer .= "<p><strong>GC Content = $obj_stats{gc_content}</strong></p>\n";

    my $srch1 = "(<replaceHeader>)(.*)(</replaceHeader>)";
    my $srch2 = "(<replaceFooter>)(.*)(</replaceFooter>)";
    my $srch3 = "(<replaceTableData>)(.*)(</replaceTableData>)";

    $file_content =~ s/$srch1/$rpt_header/;
    $file_content =~ s/$srch2/$rpt_footer/;
    $file_content =~ s/$srch3/$rpt_data/;

    my $report_file_path = catfile($self->{metag_dir}, 'genome_report.html');
    my $fh2 = $self->_openWrite($report_file_path);
    print $fh2 $file_content;
    close $fh2;
    #print $file_content;

    my($vol, $f_path, $rfile) = splitpath($report_file_path);
    my @html_report = ({'path'=> $report_file_path,
                        'name'=> $rfile,
                        'label'=> $rfile,
                        'description'=> $report_title});
    return @html_report;
}

#Create a KBaseReport with brief info/stats on a reannotated metagenome
sub _generate_report {
    my ($self, $src_ref, $ama_ref, $src_gff_conts, $ama_gff_conts, $func_tab) = @_;

    my $gn_info = $self->_fetch_object_info($ama_ref);
    my $gn_ws = $gn_info->[7];

    my %src_stats = $self->_generate_stats_from_ama($src_ref);
    my %ama_stats = $self->_generate_stats_from_ama($ama_ref);
    my %src_gff_stats = $self->_generate_stats_from_gffContents($src_gff_conts);
    my %ama_gff_stats = $self->_generate_stats_from_gffContents($ama_gff_conts);
    my %subsys_info = $self->_fetch_subsystem_info();

    my $report_message;
    my ($src_gene_count, $src_role_count, $ama_gene_count, $ama_role_count);
    if (keys %ama_gff_stats) {
        if (keys %src_gff_stats) {
            $src_gene_count = keys $src_gff_stats{gene_role_map};
            $src_role_count = keys $src_gff_stats{function_roles};
        }
        else {
            $src_gene_count = 0;
            $src_role_count = 0;
        }
        $ama_gene_count = keys $ama_gff_stats{gene_role_map};
        $ama_role_count = keys $ama_gff_stats{function_roles};
        $report_message = ("Genome Ref: $ama_ref\n".
                           "Genome type: $ama_stats{genome_type}\n".
                           "Number of contigs: $ama_stats{contig_count}\n".
                           "Number of features before RAST: $src_stats{num_features}\n".
                           "Number of unique function roles before RAST: $src_role_count\n".
                           "Number of features after RAST: $ama_stats{num_features}\n".
                           "Number of unique function roles after RAST: $ama_role_count\n");
        # Note that $ama_stats{feature_counts}
        #           $ama_gff_stats{function_roles}{gene_count}
        #           $ama_gff_stats{function_roles}{gene_list}
        # will be reported in the html file.
    }
    else {
        $report_message = ("Genome Ref: $ama_ref\n".
                           "Genome type: $ama_stats{genome_type}\n".
                           "Number of features before RAST: $src_stats{num_features}\n".
                           "No data on functional roles available\n");
    }

    my @html_files = $self->_write_html_from_stats(\%ama_stats,
                                                   \%ama_gff_stats,
                                                   \%subsys_info,
                                                   $func_tab);

    my $kbr = new installed_clients::KBaseReportClient($self->{call_back_url});
    my $report_info = $kbr->create_extended_report(
        {"message"=>$report_message,
         "objects_created"=>[{"ref"=>$ama_ref, "description"=>"RAST re-annotated metagenome"}],
         "html_links"=> \@html_files,
         "direct_html_link_index"=> 0,
         "html_window_height"=> 366,
         "report_object_name"=>"kb_RAST_metaG_report_".$self->_create_uuid(),
         "workspace_name"=>$gn_ws
        });

    return {"output_genome_ref"=>$ama_ref,
            "report_name"=>$report_info->{name},
            "report_ref"=>$report_info->{ref}};
}

sub _fetch_object_data {
    my ($self, $obj_ref) = @_;
    my $ret_obj_data = {};
    eval {
        $ret_obj_data = $self->{ws_client}->get_objects2(
                             {'objects'=>[{ref=>$obj_ref}]}
                         )->{data}->[0]->{data};
    };
    if ($@) {
        croak "ERROR Workspace.get_objects2 failed: ".$@."\n";
    }
    return $ret_obj_data;
}

sub _fetch_object_info {
    my ($self, $obj_ref) = @_;
    my $obj_info = {};
    eval {
        $obj_info = $self->{ws_client}->get_object_info3(
                                 {objects=>[{ref=>$obj_ref}]}
                        )->{infos}->[0];
    };
    if ($@) {
        croak "ERROR Workspace.get_object_info3 failed: ".$@."\n";
    }
    print "INFO: object info for $obj_ref------\n".Dumper($obj_info);
    return $obj_info;
}

# create a 12 char string unique enough here
sub _create_uuid {
    my $self = shift;

    my $str_uuid = qx(uuidgen);
    chomp $str_uuid;
    $str_uuid = substr($str_uuid, -12);
    return $str_uuid;
}

##----FILE IO ----##
sub _openWrite {
    # Open a file for writing
    my ($self, $fn) = @_;
    # File encoded using UTF-8.
    open my $fh, qw(>:encoding(UTF-8)), $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub _openRead {
    # Open a file for reading
    my ($self, $fn) = @_;
    open my $fh, qw(<), $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

sub _create_metag_dir {
    my ($self, $rast_dir) = @_;
    # create the project directory for metagenome annotation
    my $mg_dir = "metag_annotation_dir_".$self->_create_uuid();
    my $dir = catfile($rast_dir, $mg_dir);

    make_path $dir; # it won't make a directory that already exists
    return $dir;
}

sub _print_fasta_gff {
    my ($self, $start, $num_lines, $f, $k) = @_;

    # Open $f to read into an array
    my $fh = $self->_openRead($f);
    my @read_lines=();
    chomp(@read_lines = <$fh>);
    close($fh);

    my $file_lines = scalar @read_lines;
    print "There are a total of ". $file_lines . " lines from file $f\n";

    $num_lines = $start + $num_lines;
    $num_lines = $num_lines < $file_lines ? $num_lines : $file_lines;
    for (my $i = $start; $i < $num_lines; $i++ ) {
        if (defined($k)) {
            print "$read_lines[$i]\n" if $read_lines[$i] =~ m/$k/;
        }
        else {
            print "$read_lines[$i]\n";
        }
    }
}

##----subs for parsing GFF by Seaver----##
sub _parse_gff {
    my ($self, $gff_filename, $attr_delimiter) = @_;
    print "Parsing GFF contents from file $gff_filename \n";

    # Open $gff_filename to read into an array
    my $fh = $self->_openRead($gff_filename);
    my @gff_lines=();
    chomp(@gff_lines = <$fh>);
    close($fh);

    print "Read in ". scalar @gff_lines . " lines from GFF file $gff_filename\n";

    my @gff_contents=();
    foreach my $current_line (@gff_lines){
        next if $current_line =~ m/^#.*$/;

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


sub _get_feature_function_lookup {
    my ($self, $features) = @_;
    my $ftr_count = scalar @{$features};

    print "INFO: Creating feature function lookup table from $ftr_count RAST features.\n";
    if ($ftr_count > 10) {
        print "INFO: First 10 RAST feature examples:\n".Dumper(@{$features}[0..9]);
    }
    else {
        print "INFO:All $ftr_count RAST features:\n".Dumper(@{$features});
    }

    #Feature Lookup Hash
    my %function_lookup = ();
    foreach my $ftr (@$features){
        next if (!exists($ftr->{'functions'}) && !exists($ftr->{'function'}));

        if (exists($ftr->{'functions'})) {
            $function_lookup{$ftr->{'id'}}{functions}=join(" / ",@{$ftr->{'functions'}});
        }
        elsif (exists($ftr->{'function'})) {
            $function_lookup{$ftr->{'id'}}{functions}=$ftr->{'function'};
        }
        if (exists($ftr->{'annotations'})) {
            #print "\nFound annotation source array:\n".Dumper($ftr->{'annotations'});
            $function_lookup{$ftr->{'id'}}{'annotation_src'}=$ftr->{'annotations'}->[0]->[1];
        }
    }
    my $ksize = keys %function_lookup;
    print "INFO: Feature function look up table contains $ksize entries.\n";
    return %function_lookup;
}

sub _update_gff_functions_from_features {
    my ($self, $gff_contents, $feature_lookup) = @_;

    print "INFO: Updating ".scalar @{$gff_contents}." GFFs with RAST features.\n";

    my %ftrs_function_lookup = %{ $feature_lookup };
    my @new_gff_contents=();
    foreach my $gff_line (@$gff_contents) {
        if(scalar(@$gff_line) < 9){
            #No attributes column
            push(@new_gff_contents,$gff_line);
            next;
        }

        my $ftr_attributes = $gff_line->[8];

        # According to the genome loading code, the function must be added to the product attribute (go figure)
        # https://github.com/kbaseapps/GenomeFileUtil/blob/master/lib/GenomeFileUtil/core/FastaGFFToGenome.py#L665-L666
        # Note that this overwrites anything that was originally in the 'product' field if it previously existed
        # Also note that I'm forcing every feature to have at least an empty product field
        if (!defined($ftr_attributes->{'product'})) {
            $ftr_attributes->{'product'}='';
        }

        #Look for, and add function
        if(exists($ftrs_function_lookup{$ftr_attributes->{'id'}})) {
            $ftr_attributes->{'product'}=$ftrs_function_lookup{$ftr_attributes->{'id'}}->{'functions'};
        }

        $gff_line->[8] = $ftr_attributes;
        push(@new_gff_contents, $gff_line);
    }
    print "INFO: Updated new_gff_contents has ". scalar @new_gff_contents . " entries\n";
    return \@new_gff_contents;
}

sub _write_gff {
    my ($self, $gff_contents, $gff_filename, $attr_delimiter) = @_;

    # Open $gff_filename to write the @$gff_contents array back to the file
    my $fh = $self->_openWrite($gff_filename);

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
    unless (-e $gff_filename) {
        croak "**In _write_gff ERROR: could not find file $gff_filename\n";
    }
    unless (-s $gff_filename) {
        croak "**In _write_gff ERROR: empty file $gff_filename\n";
    }
}


sub _parse_fasta {
    my ($self, $fasta_filename) = @_;
    print "Parsing FASTA contents from file $fasta_filename \n";

    my @fasta_lines = ();
    # Open $fasta_filename to read into an array
    my $fh = $self->_openRead($fasta_filename);
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
    my ($self, $fasta_contents, $gff_contents) = @_;

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
    #print "Gene sequence hash example:\n".Dumper(\%gene_seqs);
    return \%gene_seqs;
}

sub _translate_gene_to_protein_sequences {
    my $self = shift;
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

##----main function----##
sub rast_metagenome {
    my $self = shift;
    my($inparams) = @_;

    print "rast_metagenome input parameter=\n". Dumper($inparams). "\n";
    
    my $params = $self->_check_annotation_params($inparams);
    my $input_obj_ref = $params->{object_ref};

    my $inputgenome = {
        features => []
    };

    my $output_type = 'gff';
    my $gff_filename = catfile($self->{metag_dir}, 'genome.gff');
    my $input_fasta_file = catfile($self->{metag_dir}, 'prodigal_input.fasta');
    my $trans_file = catfile($self->{metag_dir}, 'protein_translation');
    my $nuc_file = catfile($self->{metag_dir}, 'nucleotide_seq');
    my $output_file = catfile($self->{metag_dir}, 'prodigal_output').'.'.$output_type;
    # my $start_file = catfile($self->{metag_dir}, 'start_file');
    # cannot specify metagenomic sequence with a training file
    # my $training_file = catfile($self->{metag_dir}, 'training_file');

    # 1. getting the fasta file from $input_obj_ref according to its type
    my $input_obj_info = $self->_fetch_object_info($input_obj_ref);
    my $in_type = $input_obj_info->[2];
    my $is_assembly = $in_type =~ /GenomeAnnotations.Assembly/;
    my $is_meta_assembly = $in_type =~ /Metagenomes.AnnotatedMetagenomeAssembly/;

    if ($is_assembly) {
        # object is itself an assembly
        my $out_file = $self->_get_fasta_from_assembly($input_obj_ref);
        copy($out_file, $input_fasta_file) || croak "Copy file failed: $!\n";
    }
    elsif ($is_meta_assembly) {
        # input_obj_ref points to a genome
        if (defined($input_obj_info->[10])) {
            my $num_ftrs = $input_obj_info->[10]->{'Number features'};
            print "Input object '$input_obj_ref' is a metagenome and has $num_ftrs features.\n";
        }

        $input_fasta_file = $self->_write_fasta_from_metagenome(
                            $input_fasta_file, $input_obj_ref);
        unless (-e $input_fasta_file) {
            croak "**rast_metagenome ERROR: could not find FASTA file\n";
        }
    }
    else {
        croak ("Only KBaseMetagenomes.AnnotatedMetagenomeAssembly and ".
               "KBaseGenomeAnnotations.Assembly will be annotated by this app.\n");
    }

    # 2. fetching the gff contents, if $input_obj_ref points to an assembly, call Prodigal
    my ($fasta_contents, $gff_contents, $attr_delimiter) = ([], [], "=");

    if ($is_assembly) {
        # object is itself an assembly
        my $mode = 'meta';
        my @prodigal_cmd = $self->_build_prodigal_cmd($input_fasta_file,
                                                      $trans_file,
                                                      $nuc_file,
                                                      $output_file,
                                                      $output_type,
                                                      $mode);

        # 2.1 filing the feature list for rast call
        if ($self->_run_prodigal(@prodigal_cmd) == 0) {
            print "Prodigal finished run, files are written into:\n$output_file\n$trans_file\n$nuc_file\n";

            # print "First 10 lines of the GFF file from Prodigal-----------\n";
            # $self->_print_fasta_gff(0, 10, $output_file);
            my %transH;

            ($gff_contents, %transH) = $self->_parse_prodigal_results(
                                                  $trans_file,
                                                  $output_file,
                                                  $output_type);

            my $count = @$gff_contents;
            print "Prodigal returned $count entries.\n";

            foreach my $entry (@$gff_contents) {
                my ($contig, $source, $ftr_type, $beg, $end, $score, $strand,
                    $phase, $attribs) = @$entry;

                next if $contig =~ m/^#.*/;

                my ($seq, $trunc_left, $trunc_right) = @{$transH{"$contig\t$beg\t$end\t$strand"}};
                my $id = $attribs->{id};
                my $start = ($strand eq q(+)) ? $beg : $end;

                push(@{$inputgenome->{features}}, {
                        id                  => $id,
                        type                => $ftr_type,
                        location            => [[ $contig, $start, $strand, $end ]],
                        annotator           => $source,
                        annotation          => 'Add feature called by PRODIGAL',
                        protein_translation => $seq
                });
            }
            #print "******inputgenome data:******* \n".Dumper($inputgenome);
            copy($output_file, $gff_filename);
        }
    }
    elsif ($is_meta_assembly) {
        # generating the gff file directly from metagenome
        $gff_filename = $self->_write_gff_from_metagenome($gff_filename, $input_obj_ref);
        unless (-e $gff_filename) {
            croak "**rast_metagenome ERROR: could not find GFF file\n";
        }

        # 2.2 filing the feature list for rast call

        # fetch protein sequences and gene IDs from fasta and gff files
        $fasta_contents = $self->_parse_fasta($input_fasta_file);
        ($gff_contents, $attr_delimiter) = $self->_parse_gff($gff_filename, $attr_delimiter);

        my $gene_seqs = $self->_extract_cds_sequences_from_fasta($fasta_contents, $gff_contents);
        my $protein_seqs = $self->_translate_gene_to_protein_sequences($gene_seqs);

        foreach my $gene (sort keys %$protein_seqs){
            push(@{$inputgenome->{features}},{
                id => $gene,
                protein_translation => $protein_seqs->{$gene}
            });
        }
    }

    # 3. call RAST to annotate the proteins/genome
    my $ftr_count = scalar @{$inputgenome->{features}};
    unless ($ftr_count >= 1) {
        print "Empty input genome features, skip rasting, return original genome object.\n";
        return $input_obj_ref;
    }

    #print "--------Print out 1000 (4200~5200) lines in the GFF file before RASTing-------\n";
    #$self->_print_fasta_gff(4200, 1000, $gff_filename);

    my $rasted_genome = $self->_run_rast($inputgenome);
    my $ftrs = $rasted_genome->{features};
    my $rasted_ftr_count = scalar @{$ftrs};
    print "RAST resulted ".$rasted_ftr_count." features.\n";

    print "***********The first 10 or fewer rasted features, for example***************\n";
    my $prnt_lines = ($rasted_ftr_count > 10) ? 10 : $rasted_ftr_count;
    for (my $j=0; $j<$prnt_lines; $j++) {
        my $f_id = $ftrs->[$j]->{id};
        my $f_func = defined($ftrs->[$j]->{function}) ? $ftrs->[$j]->{function} : '';
        my $f_protein = defined($ftrs->[$j]->{protein_translation}) ? $ftrs->[$j]->{protein_translation} : '';
        print "$f_id\t$f_func\t$f_protein\n";
    }

    my %ftr_func_lookup = $self->_get_feature_function_lookup($ftrs);
    my $updated_gff_contents = $self->_update_gff_functions_from_features(
                                   $gff_contents, \%ftr_func_lookup);
    my $new_gff_file = catfile($self->{metag_dir}, 'new_genome.gff');
    $self->_write_gff($updated_gff_contents, $new_gff_file, $attr_delimiter);

    # 4. save rast re-annotated fasta/gff data
    my $out_metag = $self->_save_metagenome($params->{output_workspace},
                                            $params->{output_metagenome_name},
                                            $input_obj_ref, $new_gff_file);
    my $ama_ref = $out_metag->{metagenome_ref};
    my $report_ret = $self->_generate_report(
                         $input_obj_ref, $ama_ref, $gff_contents,
                         $updated_gff_contents, \%ftr_func_lookup);
    return $report_ret;
}

sub doInitialization {
    my $self = shift;

    $self->{_token} = $self->{ctx}->token();
    $self->{_username} = $self->{ctx}->user_id();
    $self->{_method} = $self->{ctx}->method();
    $self->{_provenance} = $self->{ctx}->provenance();

    $self->{ws_url} = $self->{config}->{'workspace-url'};
    $self->{call_back_url} = $ENV{ SDK_CALLBACK_URL };
    $self->{rast_scratch} = $self->{config}->{'scratch'};
    $self->{metag_dir} = $self->_create_metag_dir($self->{rast_scratch});

    die "no workspace-url defined" unless $self->{ws_url};

    $self->{ws_client} = new installed_clients::WorkspaceClient(
                             $self->{ws_url}, token => $self->{_token});

    return 1;
}

sub new {
    my $class = shift;

    my $self = {
        'config' => shift,
        'ctx' => shift
    };

    bless $self, $class;
    $self->doInitialization();

    return $self;             # Return the reference to the hash.
}


1;

