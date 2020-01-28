# -*- perl -*-

# package metag_utils;

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
use Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl;
use Bio::KBase::GenomeAnnotation::Service;
use Workspace::WorkspaceClient;
use AssemblyUtil::AssemblyUtilClient;
use GenomeFileUtil::GenomeFileUtilClient;

use Data::Dumper;
use Carp;
use File::Spec::Functions qw(catfile);
use File::Copy;

use gjoseqlib;


my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = Config::IniFiles->new(-file=>$config_file);
my $ws_url = $config->{'workspace-url'};

my $ws_client = new Workspace::WorkspaceClient($ws_url,token => $token);
my $call_back_url = $ENV{ SDK_CALLBACK_URL };

my $au = new AssemblyUtil::AssemblyUtilClient($call_back_url);
my $gfu = new GenomeFileUtil::GenomeFileUtilClient($call_back_url);
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

    my $output = $au->get_assembly_as_fasta({"ref" => $ref});

    copy($output->{path}, $fasta_file) || die "Could not find file:".$output->{path};
	# print("Genome is smaller than 25000 bp -- using 'metagenome' mode\n");
    my $mode = 'meta';

    $prd_params = {
        input_file => $fasta_file,       # -i (FASTA/Genbank file)
        trans_fle => $trans_file,        # -a (Write protein translations to $trans_file)
        nuc_file => $nuc_file,           # -d (Write nucleotide sequences of genes to $nuc_file)
        output_type => 'gff',          # -f (gbk, gff, or sco, default to gbk)
        output_file => $output_file,     # -o (Specify output file (default writes to stdout).) 
        closed_ends => 1,              # -c (Closed ends.  Do not allow genes to run off edges.)
        N_as_masked_seq => 1,          # -m (Treat runs of N as masked sequence; don't build genes across them.)
        trans_table => 11,             # -g (Specify a translation table to use (default 11).)
        procedure => $mode,            # -p (Select procedure (single or meta).  Default is single.)
        quiet => 0,                    # -q (Run quietly)
        start_file => $start_file,       # -s (Write all potential genes (with scores) to $start_file)
        training_file => $training_file  # -t (Write a training file (if none exists); otherwise, read and use $training_file)
    }
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
    while (defined($line = <$fh_sco>)) {
        chomp $line;
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
                warn "No translation found for \"$sco_file\" line: $line\n";
            }
        }
        else {
            warn "Could not parse calls for \"$sco_file\" line: $line\n";
        }
    }
    close($fh_sco);
    
    return $encoded_tbl;
}

sub rast_metagenome {
    my ($scratch, $params) = @_;
    my $input_obj_ref = $params->{object_ref};
    my $inputgenome = {
  		features => []
    }

    my $fasta_file = catfile($scratch, 'input_contigs.fasta');
    my $trans_file = catfile($scratch, 'protein_translation');
    my $nuc_file = catfile($scratch, 'nucleotide_seq');
    my $output_file = catfile($scratch, 'prodigal_out');
    # my $start_file = catfile($scratch, 'start_file');
    my $training_file = '';  # catfile($scratch, 'training_file');

	my $info = $ws_client->get_object_info([{ref=>$input_obj_ref}],0);

    # Check if input is an assembly, if so run Prodigal and parse for proteins
    my $protein_tbl = [];
    if ($info->[0]->[2] =~ /Assembly/) {
        my $prodigal_params = build_prodigal_params($input_obj_ref,
                                                    $fasta_file,
                                                    $trans_file,
                                                    $nuc_file,
                                                    $output_file,
                                                    $start_file,
                                                    $training_file);

        if (run_prodigal($prodigal_params) == 0) {
            # Prodigal finished run, files are written into $output_file/$trans_file/$nuc_file/$start_file/$training_file
            $prodigal_result = parse_prodigal_results($trans_file, $nuc_file, $output_file);

            my $count = @$prodigal_result;
            my $cur_id_suffix = 1;
            foreach my $entry (@$prodigal_result) {
                # print Data::Dumper->Dump($entry)."\n";
                my ($contig, $beg, undef, $strand, $length, $translation) = @$entry;

                my $id = join(".", “peg“, $cur_id_suffix);
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
    elsif ($info->[0]->[2] =~ /KBaseGenomes.Genome/) { 
        # input is a (meta)genome, get the genome proteins directly from the input
        my $in_genome = $ws_client->get_objects([{ref=>$input_obj_ref}])->[0]->{data};
        if (defined($genome_obj->{features})) {
            push(@{$inputgenome->{features}}, $in_genome->{features};);
        }
    }
    # Call RAST to annotate the proteins/genome
    my $rast_client = Bio::kbase::kbaseenv::ga_client();
    my $genome = $rast_client->run_pipeline($inputgenome,
            {stages => [{name => ‘annotate_proteins_kmer_v2’, kmer_v2_parameters => {}},
		                {name => ‘annotate_proteins_similarity’,
                         similarity_parameters => { annotate_hypothetical_only => 1 }}]}
	);
	my $ftrs = $genome->{features};
	my $return = {};
	$return->{functions} = [];
	for (my $i=0; $i < @{$genome->{features}}; $i++) {
		$return->{functions}->[$i] = [];
		if (defined($genome->{features}->[$i]->{function})) {
			$return->{functions}->[$i] = [split(/\s*;\s+|\s+[\@\/]\s+/,$genome->{features}->[$i]->{function})];
		}
	}
}

