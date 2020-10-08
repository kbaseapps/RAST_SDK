package RAST_SDK::AnnotationUtils;

##########################################################################
# This module is built to handle the annotation of Assembly
# (KBaseGenomeAnnotations.Assembly) and/or Genome
# (KBaseGenomes/KBaseGenomeAnnotations.GenomeAnnotation) object(s)
# and create a KBaseGenomeAnnotations.Assembly/GenomeAnnotation object/set
###########################################################################

our $VERSION = '0.1.9';
use strict;
use warnings;

use Config::Simple;
use Data::Dumper;
use File::Spec::Functions qw(catfile splitpath);
use File::Path qw(make_path);
use File::Copy;
use Carp qw(croak);
use Bio::SeqIO;
use Bio::Perl;
use Bio::Tools::CodonTable;
use JSON;
use Encode qw(encode decode);
use URI::Encode qw(uri_encode uri_decode);
use File::Basename;
use Array::Utils qw(:all);
use Text::Trim qw(trim);
use Data::Structure::Util qw( unbless );
use Data::UUID;
use Try::Tiny;

use GenomeTypeObject;

use Bio::KBase::KBaseEnv;
use Bio::KBase::Utilities;
use Bio::KBase::Exceptions;
use Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl;

use installed_clients::GenomeAnnotationAPIClient;
use installed_clients::GenomeAnnotationClient;
use installed_clients::AssemblyUtilClient;
use installed_clients::GenomeFileUtilClient;
use installed_clients::WorkspaceClient;
use installed_clients::KBaseReportClient;
use installed_clients::kb_SetUtilitiesClient;


#######################Begin refactor annotate_process########################
## helper subs
sub _util_version {
	my ($self) = @_;
    return $VERSION;
}

sub _get_scientific_name_for_NCBI_taxon {
    my ($self, $tax_id, $timestamp) = @_;
    my $url = $self->{_re_url} . "/api/v1/query_results?stored_query=ncbi_fetch_taxon";
    my $content = encode_json({
        'ts' => $timestamp,
        'id'=> $tax_id . '' # make sure tax id is a string
        });
    my $req = HTTP::Request->new(POST => $url);
    $req->header('content-type', 'application/json');
    $req->content($content);
    my $ua = LWP::UserAgent->new();

    my $ret = $ua->request($req);
    if (!$ret->is_success()) {
        print("Error body from Relation Engine on NCBI taxa query:\n" .
            $ret->decoded_content({'raise_error' => 1}));
        # might want to try to parse content to json and extract error, try this for now
        die "Error contacting Relation Engine: " . $ret->status_line();
    }
    my $retjsonref = $ret->decoded_content({'raise_error' => 1, 'ref' => 1});
    my $retjson = decode_json($retjsonref);

    if (!$retjson->{count}) {
        die "No result from Relation Engine for NCBI taxonomy ID " . $tax_id;
    }
    return $retjson->{'results'}[0]{'scientific_name'};
}

#
## get the genome data via GenomeAnnotationAPIClient
#
sub _get_genome {
    my ($self, $gn_ref) = @_;
    my $ga_client = installed_clients::GenomeAnnotationAPIClient->new($self->{call_back_url});
    my $output = $ga_client->get_genome_v1({
        genomes => [{
            "ref" => $gn_ref
        }],
        ignore_errors => 1,
        no_data => 0,
        no_metadata => 0,
        downgrade => 0
    });
    my $genome = $output->{genomes}->[0];
    $genome->{data}->{'_reference'} = $genome->{info}->[6]."/".$genome->{info}->[0]."/".$genome->{info}->[4];
    my $ftr_cnt = @{$genome->{data}->{features}};
    print("Genome $genome->{data}->{'_reference'} downloaded with $ftr_cnt features.\n");
    return $genome->{data};
}

#
## get the genome gff_contents
#
sub _get_genome_gff_contents {
    my ($self, $obj) = @_;

    my $chk = $self->_validate_KB_objref_name($obj);
    return [] unless $chk->{check_passed};

    my $gff_filename = catfile($self->{genome_dir}, 'genome.gff');
    my ($gff_contents, $attr_delimiter) = ([], "=");

    # getting gff_contents from $obj_ref if it is of genome type
    my $obj_info = $self->_fetch_object_info($obj, $chk);
    return [] unless $obj_info;

    my $in_type = $obj_info->[2];
    my $is_genome = ($in_type =~ /KBaseGenomes.Genome/ ||
                     $in_type =~ /KBaseGenomeAnnotatioVns.GenomeAnnotation/);

    $obj = $obj_info->[6].'/'.$obj_info->[0].'/'.$obj_info->[4];
    if ($is_genome) {
        # generating the gff file directly from genome object ref
        $gff_filename = $self->_write_gff_from_genome($obj);
        unless (-s $gff_filename) {
            croak "**rast_genome ERROR: GFF file is empty!\n";
        }
        ($gff_contents, $attr_delimiter) = $self->_parse_gff(
            $gff_filename, $attr_delimiter);
    }
    return $gff_contents;
}


sub _get_contigs_from_fastafile {
    my ($self, $fpath) = @_;

    my @contigs;
    my $contigID_hash = {};
    my $fasta = "";
    my $fh = $self->_openRead($fpath);
    while (my $line = <$fh>) {
        $fasta .= $line;
    }
    close($fh);
    $fasta =~ s/\>([^\n]+)\n/>$1\|\|\|/g;
    $fasta =~ s/\n//g;
    my $array = [split(/\>/, $fasta)];
    for (my $i=0; $i < @{$array}; $i++) {
        if (scalar @contigs > $self->{max_contigs}) {
            Bio::KBase::Exceptions::ArgumentValidationError->throw(
                error => 'too many contigs',
                method_name => 'RAST_SDK::AnnotationUtils._get_contigs_from_fastafile'
            );
        }
        if (length($array->[$i]) > 0) {
            my $subarray = [split(/\|\|\|/,$array->[$i])];
            if (@{$subarray} == 2) {
                my $description = "unknown";
                my $id = $subarray->[0];
                if( $subarray->[0] =~ /^([^\s]+)\s(.+)$/) {
                    $id = $1;
                    $description = $2;
                }
                my $tmp_contigID = 'contigID_'.$i;
                $contigID_hash->{$tmp_contigID} = $id;
                my $contigobject = {
                    id => $tmp_contigID,
                    name => $tmp_contigID,
                    length => length($subarray->[1]),
                    md5 => Digest::MD5::md5_hex($subarray->[1]),
                    sequence => $subarray->[1],
                    description => $description
                };
                push(@contigs, $contigobject);
            }
        }
    }
    return (\@contigs, $contigID_hash);
}


#
## get the contigs via WorkspaceClient and AssemblyUtilClient
## Note: the input $ref has a value which is a Workspace.ref_string
#
sub _get_contigs {
    my ($self, $obj) = @_;

    my $chk = $self->_validate_KB_objref_name($obj);
    return {} unless $chk->{check_passed};

    my $obj_info = $self->_fetch_object_info($obj, $chk);
    return {} unless $obj_info;

    $obj = $obj_info->[6].'/'.$obj_info->[0].'/'.$obj_info->[4];
    my $ret_obj = {
        _reference => $obj,
        id => $obj_info->[0],
        name => $obj_info->[1],
        source_id => $obj_info->[0],
        source => "KBase",
        type => "SingleGenome",
        contigs => []
    };

    my ($contigs, $contigID_hash);
    if ($obj_info->[2] =~ /Assembly/) {
        my $fpath = $self->_get_fasta_from_assembly($obj);
        ($contigs, $contigID_hash) = $self->_get_contigs_from_fastafile($fpath);
        $ret_obj->{contigs} = $contigs;

        my $sortedarray = [sort { $a->{sequence} cmp $b->{sequence} } @{$ret_obj->{contigs}}];
        my $str = "";
        for (my $i=0; $i < @{$sortedarray}; $i++) {
            if (length($str) > 0) {
                $str .= ";";
            }
            $str .= $sortedarray->[$i]->{sequence};
        }
        print("Assembly $ret_obj->{_reference} Downloaded\n");
        $ret_obj->{md5} = Digest::MD5::md5_hex($str);
        $ret_obj->{_kbasetype} = "Assembly";
    } elsif ($obj_info->[2] =~ /ContigSet/) {
        $ret_obj = $self->_fetch_object_data($obj);
        $ret_obj->{_kbasetype} = "ContigSet";
        $ret_obj->{_reference} = $obj;
        print("Contigset $ret_obj->{_reference} Downloaded\n");
    }
    my $totallength = 0;
    my $gclength = 0;
    for (my $i=0; $i < @{$ret_obj->{contigs}}; $i++) {
        my $newseq = $ret_obj->{contigs}->[$i]->{sequence};
        $totallength += length($newseq);
        $newseq =~ s/[atAT]//g;
        $gclength += length($newseq);
    }
    if( $totallength ) {
        $ret_obj->{_gc} = int(1000*$gclength/$totallength+0.5);
        $ret_obj->{_gc} = $ret_obj->{_gc}/1000;
    }
    return ($ret_obj, $contigID_hash);
}

#
#
# Forming the gene calling workflow in the following order:
# Call prodigal
# Call glimmer3
# Call rRNAs
# Call tRNA trnascan
# Call selenoproteins
# Call pyrrolysoproteins
# Call SEED repeat region
# Call strep suis repeats
# Call strep pneumo repeats
# Call crisprs
#
# resolve_overlapping_features
# Call features prophage phispy
# retain_old_anno_for_hypotheticals
#
# return the workflow stages data
#
sub _default_genecall_workflow {
    my $self = shift;
    my @stages = (
        { name => 'call_features_CDS_prodigal' },
        { name => 'call_features_CDS_glimmer3', failure_is_not_fatal => 1,
          glimmer3_parameters => { min_training_len => '2000'} },
        { name => 'call_features_rRNA_SEED' },
        { name => 'call_features_tRNA_trnascan' },
        { name => 'call_selenoproteins', failure_is_not_fatal => 1 },
        { name => 'call_pyrrolysoproteins', failure_is_not_fatal => 1 },
        { name => 'call_features_repeat_region_SEED',
              repeat_region_SEED_parameters => {} },
        { name => 'call_features_strep_suis_repeat',
              condition => '$genome->{scientific_name} =~ /^Streptococcus\s/' },
        { name => 'call_features_strep_pneumo_repeat',
              condition => '$genome->{scientific_name} =~ /^Streptococcus\s/' },
        { name => 'call_features_crispr', failure_is_not_fatal => 1 },
        { name => 'call_features_CDS_genemark' }
      );
      return { stages => \@stages };
}

sub _set_default_parameters {
    my $self = shift;
    return {
        call_features_rRNA_SEED => 1,
        call_features_tRNA_trnascan => 1,
        call_selenoproteins => 1,
        call_pyrrolysoproteins => 1,
        call_features_repeat_region_SEED => 1,
        call_features_strep_suis_repeat => 1,
        call_features_strep_pneumo_repeat => 1,
        call_features_crispr => 1,
        call_features_CDS_glimmer3 => 1,
        call_features_CDS_prodigal => 1,
        call_features_CDS_genemark => 1,
        annotate_proteins_kmer_v2 => 1,
        kmer_v1_parameters => 1,
        annotate_proteins_similarity => 1,
        resolve_overlapping_features => 1,
        call_features_prophage_phispy => 1,
        retain_old_anno_for_hypotheticals => 1
    };
}


sub _get_oldfunchash_oldtype_types {
    my ($self, $in_genome) = @_;

    my $oldfunchash = {};
    my $oldtype     = {};
    my %types = ();

    for my $ftr ( @{ $in_genome->{ features }} ) {
        if (!defined($ftr->{type}) || $ftr->{type} lt '     ') {
            if (defined($ftr->{protein_translation})) {
                $ftr->{type} = 'gene';
            } else {
                $ftr->{type} = 'other';
            }
        }

        # Reset functions in protein features to "hypothetical protein" to make them available
        # for re-annotation in RAST service (otherwise these features will be skipped).
        if (lc($ftr->{type}) eq "cds" || lc($ftr->{type}) eq "peg" ||
                ($ftr->{type} eq "gene" && defined($ftr->{protein_translation}))) {
            if (defined($ftr->{functions})){
                $ftr->{function} = join("; ", @{$ftr->{functions}});
            }
            $oldfunchash->{$ftr->{id}} = $ftr->{function};
            $ftr->{function} = "hypothetical protein";
        }
        elsif ($ftr->{type} eq "gene") {
            $ftr->{type} = 'Non-coding '.$ftr->{type};
        }
        $oldtype->{$ftr->{id}} = $ftr->{type};
        #
        #   Count the input feature types
        #
        if (exists $types{$ftr->{type}}) {
            $types{$ftr->{type}} += 1;
        } else {
            $types{$ftr->{type}} = 1;
        }
    }
    if (exists $in_genome->{non_coding_features}) {
        for my $ftr ( @{ $in_genome->{ non_coding_features } } ) {
            if (!defined($ftr->{type})) {
                $ftr->{type} = "Non-coding";
            }
            $oldtype->{$ftr->{id}} = $ftr->{type};
            #
            #   Count the input feature types
            #
            if (exists $types{"Non-coding ".$ftr->{type}}) {
                $types{"Non-coding ".$ftr->{type}} += 1;
            } else {
                $types{"Non-coding ".$ftr->{type}} = 1;
            }
        }
    } else {
        $in_genome->{non_coding_features} = [];
    }
    return ($oldfunchash, $oldtype, \%types, $in_genome);
}


# When the input object_ref is referencing to a genome object
sub _create_inputgenome_from_genome {
    my ($self, $inputgenome, $input_obj_ref) = @_;

    $inputgenome = $self->_get_genome($input_obj_ref);

    my ($oldfunchash, $oldtype, $types_ref);
    my %types = ();
    ($oldfunchash, $oldtype, $types_ref,
        $inputgenome) = $self->_get_oldfunchash_oldtype_types($inputgenome);
    %types = %{$types_ref};

    my ($contigref, $contigobj, $contigID_hash);
    if($inputgenome->{domain} !~ /Eukaryota|Plant/){
        if (defined($inputgenome->{contigset_ref})) {
            $contigref = $inputgenome->{contigset_ref};
        } elsif (defined($inputgenome->{assembly_ref})) {
            $contigref = $inputgenome->{assembly_ref};
        }
        ($contigobj, $contigID_hash) = $self->_get_contigs($contigref);
    }
    return ($inputgenome, $contigobj, $contigID_hash, $oldfunchash, $oldtype, \%types);
}


# When the input_obj_ref is referencing to an assembly/contigset object
sub _create_inputgenome_from_assembly {
    my ($self, $inputgenome, $input_obj_ref) = @_;

    $inputgenome->{assembly_ref} = $input_obj_ref;  # TOBe confirmed
    my ($contigobj, $contigID_hash) = $self->_get_contigs($input_obj_ref);
    return ($inputgenome, $contigobj, $contigID_hash);
}

## Checking on non_coding_features for ones that have undefined type field
sub _check_NC_features {
    my ($self, $genome) = @_;

    my $ncoding_features = $genome->{non_coding_features};
    my $cnt = 0;
    for my $ncoding_ftr (@{$ncoding_features}) {
        if(exists($ncoding_ftr->{type})) {
            $cnt++;
            #print "Non-conding feature type value= $ncoding_ftr->{type}\n";
        }
    }
    if ($cnt == scalar @{$ncoding_features} and $cnt > 0) {
        print "***INFO***: All $cnt non-coding features have defined field type***\n";
    }
}


## Mapping the contigIDs back to their original (long) names
sub _remap_contigIDs {
    my ($self, $contigID_hash, $gn_contig_ids) = @_;

    return $gn_contig_ids unless $contigID_hash;

    for my $ctg_id (@{$gn_contig_ids}) {
        $ctg_id = $contigID_hash->{$ctg_id};
    }
    return $gn_contig_ids;
}
## end helper subs
#

## begin building the sequence of functions refactored from annotate_process in the RASTImpl
#

#
## According to inputs in $parameters and a pre-created default genome object, $inputgenome
## Returns:
##  ($rast_details = {
##      parameters => $parameters,
##      oldfunchash => $oldfunchash,
##      oldtype => $oldtype,
##      types => \%types,
##      contigobj => $contigobj,
##      fasta_file => $fasta_file,
##      gff_file => $gff_file
##  }, $inputgenomes);
#
sub _set_parameters_by_input {
    my ($self, $parameters, $inputgenome) = @_;

    #print "_set_parameters_by_input's input parameters:\n". Dumper($parameters). "\n";
    my $contigobj;
    my $contigID_hash = {};
    my $oldfunchash = {};
    my $oldtype     = {};
    my %types = ();
    my %rast_details = ();

    # 1. getting the fasta & gff files from $input_obj_ref according to its type
    my $ws = $parameters->{output_workspace};
    my $input_obj_ref = $parameters->{object_ref};
    my $chk = $self->_validate_KB_objref_name($input_obj_ref);
    return ((), {}) unless $chk->{check_passed};
    my $input_obj_info = $self->_fetch_object_info($input_obj_ref, $chk, $ws);
    return ((), {}) unless $input_obj_info;

    my $gn_name = $parameters->{output_genome_name};
    my $in_type = $input_obj_info->[2];
    my $is_assembly = ($in_type =~ /KBaseGenomeAnnotations\.Assembly/ ||
                       $in_type =~ /KBaseGenomes\.ContigSet/);
    my $is_genome = ($in_type =~ /KBaseGenomes\.Genome/ ||
                     $in_type =~ /KBaseGenomeAnnotations\.GenomeAnnotation/);

    if (!$is_assembly && !$is_genome) {
        croak ("Only KBaseGenomes.Genome, KBaseGenomes.ContigSet, ".
               "KBaseGenomeAnnotations.Assembly, and ".
               "KBaseGenomeAnnotations.Assembly will be annotated by this app.\n");
    }

    $input_obj_ref = $input_obj_info->[6].'/'.$input_obj_info->[0].'/'.$input_obj_info->[4];

    my $types_ref;
    if ($is_genome) {
        print "INFO:----$input_obj_ref points to a genome----\n";
        if (defined($input_obj_info->[10])) {
            my $num_ftrs = $input_obj_info->[10]->{'Number features'};
            print "Input object $input_obj_ref is a genome and has $num_ftrs features.\n";
        }

        ($inputgenome, $contigobj, $contigID_hash, $oldfunchash,
            $oldtype, $types_ref) = $self->_create_inputgenome_from_genome(
                                                    $inputgenome, $input_obj_ref);
        $rast_details{oldfunchash} = $oldfunchash;
        $rast_details{oldtype} = $oldtype;
        $rast_details{types} = $types_ref;
    } elsif ($is_assembly) {
        print "INFO:----$input_obj_ref points to an assembly----\n";
        ($inputgenome, $contigobj,
            $contigID_hash) = $self->_create_inputgenome_from_assembly(
                                                    $inputgenome, $input_obj_ref);
    }

    $parameters->{genetic_code} = $inputgenome->{genetic_code};
    $parameters->{domain} = $inputgenome->{domain};
    $parameters->{scientific_name} = $inputgenome->{scientific_name};

    $rast_details{parameters} = $parameters;
    $rast_details{contigobj} = $contigobj;
    $rast_details{contigID_hash} = $contigID_hash;
    $rast_details{is_genome} = $is_genome;
    $rast_details{is_assembly} = $is_assembly;
    return (\%rast_details, $inputgenome);
}


#
## Recording the annotation process by taking a note in $message;
## Inserting contigs into $inputgenome
## Input:
##  ($rast_details = {
##      parameters => $parameters,
##      oldfunchash => $oldfunchash,
##      oldtype => $oldtype,
##      types => \%types,
##      contigobj => $contigobj,
##      fasta_file => $fasta_file,
##      gff_file => $gff_file
##  }, $inputgenome);
## Returns:
##  ($rast_details = {
##      parameters => $parameters,
##      oldfunchash => $oldfunchash,
##      oldtype => $oldtype,
##      types => \%types,
##      contigobj => $contigobj,
##      message => $message,
##      tax_domain => $tax_domain,
##      fasta_file => $fasta_file,
##      gff_file => $gff_file
##  }, $inputgenome);
#
sub _set_messageNcontigs {
    my ($self, $rast_in, $inputgenome) = @_;

    my %rast_details = %{ $rast_in };
    my $parameters = $rast_details{parameters};
    my $contigobj = $rast_details{contigobj};
    my $is_assembly = $rast_details{is_assembly};

    my $tax_domain = 'U';
    my $message = "";
    if (defined($inputgenome->{domain})) {
        $tax_domain = ($inputgenome->{domain} =~ m/^([ABV])/o) ? $inputgenome->{domain} : 'U';
    }
    if ($tax_domain eq 'U' ) {
        $message .= "Some RAST tools will not run unless the taxonomic domain is Archaea, Bacteria, or Virus. \nThese tools include: call selenoproteins, call pyrroysoproteins, call crisprs, and call prophage phispy features.\nYou may not get the results you were expecting with your current domain of $inputgenome->{domain}.\n";
    }
    if (defined($contigobj)) {
        my $count = 0;
        my $size = 0;
        if (defined($contigobj->{contigs})) {
            $inputgenome->{contigs} = $contigobj->{contigs};
            for my $input_contig ( @{$inputgenome->{contigs}} ) {
                $count++;
                $size += length($input_contig->{sequence});
                $input_contig->{dna} = delete $input_contig->{sequence};
            }
        }

        if ($contigobj->{_kbasetype} eq "ContigSet") {
            $inputgenome->{contigset_ref} = $contigobj->{_reference};
        } else {
            $inputgenome->{assembly_ref} = $contigobj->{_reference};
        }

        if ($is_assembly) {
            $message .= "The RAST algorithm was applied to annotating a genome sequence comprised of ".$count." contigs containing ".$size." nucleotides. \nNo initial gene calls were provided.\n";
        } else {
            $message .= "The RAST algorithm was applied to annotating an existing genome: ".$parameters->{scientific_name}.". \nThe sequence for this genome is comprised of ".$count." contigs containing ".$size." nucleotides. \nThe input genome has ".@{$inputgenome->{features}}." existing coding features.\n and ".@{$inputgenome->{non_coding_features}}." existing non-coding features.\n";
            $message .= "NOTE: Older input genomes did not properly separate coding and non-coding features.\n" if (@{$inputgenome->{non_coding_features}} == 0);
        }
    } else {
        if($inputgenome->{domain} !~ /Eukaryota|Plant/){
            $message .= "The RAST algorithm was applied to annotating an existing genome: ".$parameters->{scientific_name}.". \nNo DNA sequence was provided for this genome, therefore new genes cannot be called. \nWe can only functionally annotate the ".@{$inputgenome->{features}}." existing features.\n";
        } else {
            $message .= "The RAST algorithm was applied to functionally annotate ".@{$inputgenome->{features}}." coding features.\n and ".@{$inputgenome->{non_coding_features}}." existing non-coding features in an existing genome: ".$parameters->{scientific_name}.".\n";
            $message .= "NOTE: Older input genomes did not properly separate coding and non-coding features.\n" if (@{$inputgenome->{non_coding_features}} == 0);
        }
    }
    if (defined($rast_details{types})) {
        my %types = %{$rast_details{types}};
        $message .= "Input genome has the following feature types:\n";
        for my $key (sort keys(%types)) {
            $message .= sprintf("\t%-30s %5d \n", $key, $types{$key});
        }
    }

    $rast_details{message} = $message;
    $rast_details{tax_domain} = $tax_domain;
    return (\%rast_details, $inputgenome);
}


#
## Check the gene call specifications in $parameters and set the workflow accordingly
## Input:
##  ($rast_details = {
##      parameters => $parameters,
##      oldfunchash => $oldfunchash,
##      oldtype => $oldtype,
##      types => \%types,
##      contigobj => $contigobj,
##      message => $message,
##      tax_domain => $tax_domain
##  }, $inputgenome);
## Returns:
##  $rast_details = {
##      parameters => $parameters,
##      oldfunchash => $oldfunchash,
##      oldtype => $oldtype,
##      types => \%types,
##      contigobj => $contigobj,
##      message => $message,
##      tax_domain => $tax_domain
##      genecall_workflow => $workflow,
##      extragenecalls => $extragenecalls,
##      genecalls => $genecalls
##  };
#
sub _set_genecall_workflow {
    my ($self, $rast_in, $inputgenome) = @_;

    my %rast_details = %{ $rast_in };
    my $parameters = $rast_details{parameters};
    my $message = $rast_details{message};
    my $tax_domain = $rast_details{tax_domain};
    my $contigobj = $rast_details{contigobj};

    my $workflow = {stages => []};

    my $extragenecalls = "";
    if (defined($parameters->{call_features_rRNA_SEED})
            && $parameters->{call_features_rRNA_SEED} == 1) {
        if (length($extragenecalls) == 0) {
            $extragenecalls = "A scan was conducted for the following additional feature types: ";
        } else {
            $extragenecalls .= "; ";
        }
        $extragenecalls .= "rRNA";
        push(@{$workflow->{stages}},{name => "call_features_rRNA_SEED"});
        if (!defined($contigobj)) {
            Bio::KBase::Utilities::error("Cannot call genes on genome with no contigs!");
        }
    }
    if (defined($parameters->{call_features_tRNA_trnascan})
            && $parameters->{call_features_tRNA_trnascan} == 1) {
        if (length($extragenecalls) == 0) {
            $extragenecalls = "A scan was conducted for the following additional feature types: ";
        } else {
            $extragenecalls .= "; ";
        }
        $extragenecalls .= "tRNA";
        push(@{$workflow->{stages}}, {name => "call_features_tRNA_trnascan"});
        if (!defined($contigobj)) {
            Bio::KBase::Utilities::error("Cannot call genes on genome with no contigs!");
        }
    }
    if (defined($parameters->{call_selenoproteins})
            && $parameters->{call_selenoproteins} == 1) {
        if ($tax_domain ne 'U' ) {
            if (length($extragenecalls) == 0) {
                $extragenecalls = "A scan was conducted for the following additional feature types: ";
            } else {
                $extragenecalls .= "; ";
            }
            $extragenecalls .= "selenoproteins";
            push(@{$workflow->{stages}},{name => "call_selenoproteins"});
            if (!defined($contigobj)) {
                Bio::KBase::Utilities::error("Cannot call genes on genome with no contigs!");
            }
        } else {
            $message .= "Did not call selenoproteins because the tax_domain is 'U'\n\n";
        }
    }
    if (defined($parameters->{call_pyrrolysoproteins})
            && $parameters->{call_pyrrolysoproteins} == 1)   {
        if ($tax_domain ne 'U' ) {
            if (length($extragenecalls) == 0) {
                $extragenecalls = "A scan was conducted for the following additional feature types: ";
            } else {
            $extragenecalls .= "; ";
            }
            $extragenecalls .= "pyrrolysoproteins";
            push(@{$workflow->{stages}},{name => "call_pyrrolysoproteins"});
            if (!defined($contigobj)) {
                Bio::KBase::Utilities::error("Cannot call genes on genome with no contigs!");
            }
        } else {
            $message .= "Did not call pyrrolysoproteins because the tax_domain is 'U'\n\n";
        }
    }
    if (defined($parameters->{call_features_repeat_region_SEED})
            && $parameters->{call_features_repeat_region_SEED} == 1)   {
        if (length($extragenecalls) == 0) {
            $extragenecalls = "A scan was conducted for the following additional feature types: ";
        } else {
            $extragenecalls .= "; ";
        }
        $extragenecalls .= "repeat regions";
        push(@{$workflow->{stages}},{
            name => "call_features_repeat_region_SEED",
            "repeat_region_SEED_parameters" => {
                            "min_identity" => "95",
                            "min_length" => "100"
                     }
        });
        if (!defined($contigobj)) {
            Bio::KBase::Utilities::error("Cannot call genes on genome with no contigs!");
        }
    }
    if (defined($parameters->{call_features_strep_suis_repeat})
			&& defined($parameters->{scientific_name})) {
        if( $parameters->{call_features_strep_suis_repeat} == 1
                && $parameters->{scientific_name} =~ /^Streptococcus\s/) {
            if (length($extragenecalls) == 0) {
                $extragenecalls = "A scan was conducted for the following additional feature types: ";
            } else {
                $extragenecalls .= "; ";
            }
            $extragenecalls .= "strep suis repeats";
            push(@{$workflow->{stages}},{name => "call_features_strep_suis_repeat"});
            if (!defined($contigobj)) {
                Bio::KBase::Utilities::error("Cannot call genes on genome with no contigs!");
            }
        }
    }
    if (defined($parameters->{call_features_strep_pneumo_repeat})
			&& defined($parameters->{scientific_name})) {
        if ($parameters->{call_features_strep_pneumo_repeat} == 1
                && $parameters->{scientific_name} =~ /^Streptococcus\s/) {
            if (length($extragenecalls) == 0) {
                $extragenecalls = "A scan was conducted for the following additional feature types: ";
            } else {
                $extragenecalls .= "; ";
            }
            $extragenecalls .= "strep pneumonia repeats";
            push(@{$workflow->{stages}},{name => "call_features_strep_pneumo_repeat"});
            if (!defined($contigobj)) {
                Bio::KBase::Utilities::error("Cannot call genes on genome with no contigs!");
            }
        }
    }
    if (defined($parameters->{call_features_crispr})
            && $parameters->{call_features_crispr} == 1) {
        if ($tax_domain ne 'U' ) {
            if (length($extragenecalls) == 0) {
                $extragenecalls = "A scan was conducted for the following additional feature types: ";
            } else {
                $extragenecalls .= "; ";
            }
            $extragenecalls .= "crispr";
            push(@{$workflow->{stages}},{name => "call_features_crispr"});
            if (!defined($contigobj)) {
                Bio::KBase::Utilities::error("Cannot call genes on genome with no contigs!");
            }
        } else {
            $message .= "Did not call crisprs because the domain is tax_domain is 'U'\n\n";
        }
    }
    $extragenecalls .= ".\n" if (length($extragenecalls) > 0);

    my $genecalls = "";
    if (defined($parameters->{call_features_CDS_glimmer3})
            && $parameters->{call_features_CDS_glimmer3} == 1)   {
        if (@{$inputgenome->{features}} > 0) {
            $message .= "The existing gene features were cleared due to selection of gene calling with Glimmer3.\n";
#           $inputgenome->{features} = [];
        }
        if (length($genecalls) == 0) {
            $genecalls = "Standard features were called using: ";
        } else {
            $genecalls .= "; ";
        }
        $genecalls .= "glimmer3";
        push(@{$workflow->{stages}},{
            name => "call_features_CDS_glimmer3",
            "glimmer3_parameters" => {
                        "min_training_len" => "2000"}
        });
        if (!defined($contigobj)) {
            Bio::KBase::Utilities::error("Cannot train and call glimmer genes on a genome with no contigs > 2000 nt!\n");
        }
    }
    if (defined($parameters->{call_features_CDS_prodigal})
            && $parameters->{call_features_CDS_prodigal} == 1)   {
        if (@{$inputgenome->{features}} > 0) {
            $message .= "The existing gene features were cleared due to selection of gene calling with Prodigal.\n";
            #$inputgenome->{features} = [];
        }
        if (length($genecalls) == 0) {
            $genecalls = "Standard gene features were called using: ";
        } else {
            $genecalls .= "; ";
        }
        $genecalls .= "prodigal";
        push(@{$workflow->{stages}}, {name => "call_features_CDS_prodigal"});
        if (!defined($contigobj)) {
            Bio::KBase::Utilities::error("Cannot call genes on genome with no contigs!\n");
        }
    }
    $genecalls .= ".\n" if (length($genecalls) > 0);
#	if (defined($parameters->{call_features_CDS_genemark}) && $parameters->{call_features_CDS_genemark} == 1)	{
#		if (@{$inputgenome->{features}} > 0) {
##			$inputgenome->{features} = [];
#			$message .= " Existing gene features were cleared due to selection of gene calling with Glimmer3, Prodigal, or Genmark.";
#		}
#		if (length($genecalls) == 0) {
#			$genecalls = "Standard gene features were called using: ";
#		} else {
#			$genecalls .= "; ";
#		}
#		$genecalls .= "genemark";
#		push(@{$workflow->{stages}},{name => "call_features_CDS_genemark"});
#		if (!defined($contigobj)) {
#			Bio::KBase::Utilities::error("Cannot call genes on genome with no contigs!");
#		}
#	}

    $rast_details{genecall_workflow} = $workflow;
    $rast_details{message} = $message;
    $rast_details{extragenecalls} = $extragenecalls;
    $rast_details{genecalls} = $genecalls;

    return \%rast_details;
}

#
## Create workflow with only the annotation part (i.e., no gene callings)
## Input:
##  $rast_details = {
##      parameters => $parameters,
##      oldfunchash => $oldfunchash,
##      oldtype => $oldtype,
##      types => \%types,
##      contigobj => $contigobj,
##      message => $message,
##      tax_domain => $tax_domain
##      genecall_workflow => $workflow,
##      extragenecalls => $extragenecalls,
##      genecalls => $genecalls
##  };
## Returns:
##  $rast_details = {
##      parameters => $parameters,
##      oldfunchash => $oldfunchash,
##      oldtype => $oldtype,
##      types => \%types,
##      contigobj => $contigobj,
##      message => $message,
##      tax_domain => $tax_domain
##      genecall_workflow => $gc_workflow,
##      annotate_workflow => $workflow,
##      extragenecalls => $extragenecalls,
##      genecalls => $genecalls,
##      annomessage => $annomessage
##  };
#
sub _set_annotation_workflow {
    my ($self, $rast_in) = @_;

    my %rast_details = %{ $rast_in };
    my $parameters = $rast_details{parameters};
    my $message = $rast_details{message};
    my $tax_domain = $rast_details{tax_domain};

    my $annomessage = "";
    my $v1flag = 0;
    my $simflag = 0;
    my $workflow = {stages => []};
    if (defined($parameters->{annotate_proteins_kmer_v2})
            && $parameters->{annotate_proteins_kmer_v2} == 1) {
        if (length($annomessage) == 0) {
            $annomessage = "The genome features were functionally annotated using the following algorithm(s): ";
        }
        $annomessage .= "Kmers V2";
        $v1flag = 1;
        $simflag = 1;
        push(@{$workflow->{stages}},{
            name => "annotate_proteins_kmer_v2",
            "kmer_v2_parameters" => {
                "min_hits" => "5",
                "annotate_hypothetical_only" => 1
            }
        });
    }
    if (defined($parameters->{kmer_v1_parameters})
            && $parameters->{kmer_v1_parameters} == 1) {
        $simflag = 1;
        if (length($annomessage) == 0) {
            $annomessage = "The genome features were functionally annotated using the following algorithm(s): ";
        } else {
            $annomessage .= "; ";
        }
        $annomessage .= "Kmers V1";
        push(@{$workflow->{stages}},{
            name => "annotate_proteins_kmer_v1",
            "kmer_v1_parameters" => {
                "dataset_name" => "Release70",
                "annotate_hypothetical_only" => $v1flag
            }
        });
    }
    if (defined($parameters->{annotate_proteins_similarity})
            && $parameters->{annotate_proteins_similarity} == 1) {
        if (length($annomessage) == 0) {
            $annomessage = "The genome features were functionally annotated using the following algorithm(s): ";
        } else {
            $annomessage .= "; ";
        }
        $annomessage .= "protein similarity";
        push(@{$workflow->{stages}},{
            name => "annotate_proteins_similarity",
            "similarity_parameters" => {
                "annotate_hypothetical_only" => $simflag
            }
        });
    }
    if (defined($parameters->{resolve_overlapping_features})
            && $parameters->{resolve_overlapping_features} == 1) {
        push(@{$workflow->{stages}}, {
            name => "resolve_overlapping_features",
            "resolve_overlapping_features_parameters" => {}
        });
    }
    if (defined($parameters->{call_features_prophage_phispy})
            && $parameters->{call_features_prophage_phispy} == 1) {
        if ($tax_domain ne 'U' ) {
            push(@{$workflow->{stages}}, {name => "call_features_prophage_phispy"});
        } else {
            $message .= "Did not call call features prophage phispy because tax_domain is 'U'\n\n";
        }
    }
    $annomessage .= ".\n" if (length($annomessage) > 0);

    $rast_details{annotate_workflow} = $workflow;
    $rast_details{message} = $message;
    $rast_details{annomessage} = $annomessage;

    return \%rast_details;
}

#
## Merge messages, update $inputgenome and $workflow
## Input:
##  ($rast_details = {
##      parameters => $parameters,
##      oldfunchash => $oldfunchash,
##      oldtype => $oldtype,
##      types => \%types,
##      contigobj => $contigobj,
##      message => $message,
##      tax_domain => $tax_domain
##      genecall_workflow => $gc_workflow,
##      annotate_workflow => $workflow,
##      extragenecalls => $extragenecalls,
##      genecalls => $genecalls,
##      annomessage => $annomessage
##  }, $inputgenome);
## Returns:
##  ($rast_details = {
##      parameters => $parameters,
##      oldfunchash => $oldfunchash,
##      oldtype => $oldtype,
##      types => \%types,
##      contigobj => $contigobj,
##      message => $message,
##      tax_domain => $tax_domain
##      genecall_workflow => $gc_workflow,
##      annotate_workflow => $anno_workflow,
##      extragenecalls => $extragenecalls,
##      renumber_workflow => $workflow,
##      genecalls => $genecalls,
##      annomessage => $annomessage
##  }, $inputgenome);
#
sub _renumber_features {
    my ($self, $rast_ref, $inputgenome) = @_;
    my %rast_details = %{ $rast_ref };
    my $message = $rast_details{message};
    my $annomessage = $rast_details{annomessage};
    my $genecalls = $rast_details{genecalls};
    my $extragenecalls = $rast_details{extragenecalls};

    my $workflow = {stages => []};

    if (length($genecalls) > 0) {
        push(@{$workflow->{stages}},{name => "renumber_features"});

        if (@{$inputgenome->{features}} > 0) {
            my $replace = [];
            for (my $i=0; $i< scalar @{$inputgenome->{features}}; $i++) {
                my $ftr = $inputgenome->{features}->[$i];
                if (!defined($ftr->{protein_translation})
                        || $ftr->{type} =~ /pseudo/) {
                    push(@{$replace}, $ftr);
                }
            }
            $inputgenome->{features} = $replace;
        }
        $message .= $genecalls;
    }
    if (length($extragenecalls) > 0) {
        $message .= $extragenecalls;
    }
    if (length($annomessage) > 0) {
        $message .= $annomessage;
    }
    $rast_details{message} = $message;
    $rast_details{renumber_workflow} = $workflow;

    return (\%rast_details, $inputgenome);
}

#
## Final massage $inputgenome and $workflows
## Input:
##  ($rast_details = {
#      parameters => $parameters,
#      oldfunchash => $oldfunchash,
#      oldtype => $oldtype,
#      types => \%types,
#      contigobj => $contigobj,
#      message => $message,
#      tax_domain => $tax_domain
#      genecall_workflow => $gc_workflow,
#      annotate_workflow => $anno_workflow,
#      extragenecalls => $extragenecalls,
#      renumber_workflow => $workflow,
#      genecalls => $genecalls,
#      annomessage => $annomessage
#  }, $inputgenome);
## Returns:
##  $rast_details = {
#      parameters => $parameters,
#      oldfunchash => $oldfunchash,
#      oldtype => $oldtype,
#      types => \%types,
#      contigobj => $contigobj,
#      message => $message,
#      tax_domain => $tax_domain
#      genecall_workflow => $gc_workflow,
#      annotate_workflow => $anno_workflow,
#      extragenecalls => $extragenecalls,
#      renumber_workflow => $workflow,
#      genecalls => $genecalls,
#      annomessage => $annomessage,
#      genehash => $genehash
#  }, $inputgenome);
#
sub _pre_rast_call {
    my ($self, $rast_ref, $inputgenome) = @_;

    my %rast_details = %{ $rast_ref };
    my $genome = $inputgenome;
    my $genehash = {};
    if (defined($genome->{features})) {
        for (my $i=0; $i < @{$genome->{features}}; $i++) {
            # Caching feature functions for future comparison against new functions
            # defined by RAST service in order to calculate number of updated features.
            # If function is not set we treat it as empty string to avoid perl warning.
            my $func = $genome->{features}->[$i]->{function};
            if (not defined($func)) {
                $func = "";
            }
            $genehash->{$genome->{features}->[$i]->{id}}->{$func} = 1;
        }
    }
    if (defined($genome->{non_coding_features})) {
        for (my $i=0; $i < @{$genome->{non_coding_features}}; $i++) {
            $genehash->{$genome->{non_coding_features}->[$i]->{id}} = 1;
        }
    }
    if (defined($inputgenome->{features})) {
        ## Checking if it's recent old genome containing both CDSs and genes in "features" array
        if ((not defined($inputgenome->{cdss})) && (not defined($inputgenome->{mrnas}))) {
            my $genes = [];
            my $non_genes = [];
            for (my $i=0; $i < @{$inputgenome->{features}}; $i++) {
                my $feature = $inputgenome->{features}->[$i];
                if ($feature->{type} eq "gene"
                        && (not defined($feature->{protein_translation}))) {
                    # gene without protein translation
                    push(@{$genes}, $feature);
                } else {
                    push(@{$non_genes}, $feature);
                }
            }
            if ((scalar @{$genes} > 0) && (scalar @{$non_genes} > 0)) {
                # removing genes from features
                $inputgenome->{features} = $non_genes;
            }
        }
        ## Temporary switching feature type from 'gene' to 'CDS':
        for (my $i=0; $i < @{$inputgenome->{features}}; $i++) {
            my $feature = $inputgenome->{features}->[$i];
            if ($feature->{type} eq "gene" and defined($feature->{protein_translation})) {
                $feature->{type} = "CDS";
            }
        }
        # When RAST annotates the genome it becomes incompatible with the new
        # spec. Removing this attribute triggers an "upgrade" to the genome
        # that fixes these problems when saving with GFU
        delete $inputgenome->{feature_counts};
        if (defined($inputgenome->{ontology_events})) {
            my $ont_event = {
                "id" => "SSO",
                "method" => Bio::KBase::Utilities::method(),
                "method_version" => $self->_util_version(),
                "ontology_ref" => "KBaseOntology/seed_subsystem_ontology",
                "timestamp" => Bio::KBase::Utilities::timestamp()
            };
            push(@{$inputgenome->{ontology_events}}, $ont_event);
        }
    }

    $rast_details{genehash} = $genehash;
    return (\%rast_details, $inputgenome);
}
#
## end building the sequence of functions refactored from annotate_process in the RASTImpl
#
#
#--------------------------------------------------------------------------
## Runs RAST with input workflow for gene calls and/or annotation on genome
#--------------------------------------------------------------------------
#
sub _run_rast_workflow {
    my ($self, $in_genome, $workflow) = @_;

    my $count = scalar @{$in_genome->{features}};
    print "******INFO: Run RAST pipeline on $in_genome->{id} with $count features.******\n";

    my $rasted_gn = $in_genome;
    my $g_data_type = (ref($rasted_gn) eq 'HASH') ? 'ref2Hash' : ref($rasted_gn);
    if ($g_data_type eq 'GenomeTypeObject') {
        print "**********Genome input passed to the run_rast_workflow with:\n".
              Dumper($workflow)." is of type of $g_data_type, prepare it**********.\n";
        $rasted_gn = $rasted_gn->prepare_for_return();
    }
    try {
        my $rast_client = Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl->new();
        $rasted_gn = $rast_client->run_pipeline($rasted_gn, $workflow);
        print "********SUCCEEDED: calling rast run_pipeline with\n".Dumper($workflow).
              "\non $in_genome->{id}.\n";
    } catch {
        print "********ERROR calling rast run_pipeline with\n".Dumper($workflow).
              "\non $in_genome->{id}:\n$_\n";
        $rasted_gn = $in_genome;
    };
    return $rasted_gn;
}

#
## process the rast genecall & annotation result
#
sub _post_rast_ann_call {
    my ($self, $genome, $inputgenome, $parameters, $contigobj) = @_;

    delete $genome->{contigs};
    delete $genome->{feature_creation_event};
    delete $genome->{analysis_events};
    if (defined($genome->{genetic_code})) {
        $genome->{genetic_code} = $genome->{genetic_code}+0;
    }
    $genome->{id} = $parameters->{output_genome_name};
    if (!defined($genome->{source})) {
        $genome->{source} = "KBase";
        $genome->{source_id} = $parameters->{output_genome_name};
    }
    if (defined($inputgenome->{gc_content})) {
        $genome->{gc_content} = $inputgenome->{gc_content};
    }
    if (defined($genome->{gc})) {
        $genome->{gc_content} = $genome->{gc}+0;
        delete $genome->{gc};
    }
    if (!defined($genome->{gc_content})) {
        $genome->{gc_content} = 0.5;
    }
    if (not defined($genome->{non_coding_features})) {
        $genome->{non_coding_features} = [];
    }
    $genome = $self->_move_non_coding_features($genome, $contigobj);
    return $genome;
}

sub _move_non_coding_features {
    my ($self, $genome, $contigobj) = @_;

    my @splice_list = ();
    if (defined($genome->{features})) {
        for (my $i=0; $i < scalar @{$genome->{features}}; $i++) {
            my $ftr = $genome->{features}->[$i];
            if (defined($ftr->{aliases}) && scalar @{$ftr->{aliases}} > 0)  {
                if (ref($ftr->{aliases}->[0]) !~ /ARRAY/) {
                # Found some pseudogenes that have wrong structure for aliases
                    my $tmp = [];
                    for my $key (@{$ftr->{aliases}})  {
                        my @ary = ('alias', $key);
                        push(@{$tmp}, \@ary);
                    }
                    $ftr->{aliases} = $tmp;
                }
            }
            if (defined ($ftr->{type}) && $ftr->{type} ne 'gene' && $ftr->{type} ne 'CDS') {
                push(@splice_list, $i);
            }
        }
    }

    my $count = 0;
    #
    #   Move non-coding features from features to non_coding_features
    #   They can have more than one location and they need md5 and dna_sequence_length
    #
    #   print "Array size =".scalar @{$genome->{features}}."\n";
    #   print "Number to splice out = ".scalar @splice_list."\n";
    for my $key (reverse @splice_list) {
        if ($key =~ /\D/ ) {
            print "INVALID:size=".scalar @{$genome->{features}}."\n";
            # print Dumper $key;
        }
        else {
            my $ftr = $genome->{features}->[$key];
            if (defined($ftr->{location})) {
                $ftr->{dna_sequence_length} = 0;
                $ftr->{md5} = '';
                for (my $i=0; $i < scalar @{$ftr->{location}}; $i++) {
                    $ftr->{location}->[$i]->[1]  = $ftr->{location}->[$i]->[1]+0;
                    $ftr->{location}->[$i]->[3]  = $ftr->{location}->[$i]->[3]+0;
                    $ftr->{dna_sequence_length} += $ftr->{location}->[$i]->[3]+0;
                }
            }
            delete $ftr->{feature_creation_event};
            my $non = splice(@{$genome->{features}}, $key, 1);
            push(@{$genome->{non_coding_features}}, $non);

            $count++;
        }
    }

    if ($count > 0) {
        print "\n***INFO:***$count features have been moved to non_coding_features.\n";
    }
    my $nc_ftr_count = @{$genome->{non_coding_features}};
    print "\n***INFO:***There are $nc_ftr_count non_coding_features.\n";

    if (defined($contigobj) && defined($contigobj->{contigs})
            && scalar(@{$contigobj->{contigs}})>0 ) {
        $genome->{num_contigs} = @{$contigobj->{contigs}};
        $genome->{md5} = $contigobj->{md5};
    }
    return $genome;
}

## Getting the seed ontology dictionary
sub _build_seed_ontology {
    my ($self, $rast_ref, $genome, $inputgenome) = @_;

    print "**************_build_seed_ontology*************\n";
    my %rast_details  = %{ $rast_ref };
    my $genehash = $rast_details{genehash};
    my $parameters = $rast_details{parameters};
    my $oldfunchash = $rast_details{oldfunchash};
    my $oldtype = $rast_details{oldtype};

    my $output = $self->{ws_client}->get_objects2({
                            'objects'=>[{
                                workspace => "KBaseOntology",
                                name => "seed_subsystem_ontology"
                            }]
                        })->{data}->[0]->{data};

    #Building a hash of standardized seed function strings
    my $num_coding = 0;
    my $funchash = {};
    for my $term (keys(%{$output->{term_hash}})) {
        my $rolename = lc($output->{term_hash}->{$term}->{name});
        $rolename =~ s/[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+//g;
        $rolename =~ s/\s//g;
        $rolename =~ s/\#.*$//g;
        $funchash->{$rolename} = $output->{term_hash}->{$term};
    }
    my $newftrs = 0;
    my $newncfs = 0;
    my $update_cdss = 'N';
    my $proteins = 0;
    my $others = 0;
    my $seedfunctions = 0;
    my $genomefunchash;
    my $seedfunchash;
    my $advancedmessage = '';
    my %types = ();
    if (defined($genome->{features})) {
        for (my $i=0; $i < @{$genome->{features}}; $i++) {
            my $ftr = $genome->{features}->[$i];
            if (defined($genehash) && !defined($genehash->{$ftr->{id}})) {
                # Let's count number of features with functions updated by RAST service.
                # If function is not set we treat it as empty string to avoid perl warning.
                $newftrs++;
                my $func = $ftr->{function};
                if (not defined($func)) {
                    $func = "";
                }
                if ($ftr->{type} eq 'CDS' && !defined($genome->{cdss}->[$i])) {
                    #Some of the optional gene finders don't respect new genome format
                    #They may not have the cdss
                    $update_cdss = 'Y';
                }
            } else {
                $num_coding++;
            }
            if (!defined($ftr->{type}) && $ftr->{id} =~ m/(\w+)\.\d+$/) {
                $ftr->{type} = $1;
            }
            if (defined($ftr->{protein_translation})) {
                $ftr->{protein_translation_length} = length($ftr->{protein_translation})+0;
                $ftr->{md5} = Digest::MD5::md5_hex($ftr->{protein_translation});
                $proteins++;
            } else {
                $others++;
            }
            if (defined($ftr->{dna_sequence})) {
                $ftr->{dna_sequence_length} = length($ftr->{dna_sequence})+0;
            }
            if (defined($ftr->{quality}->{weighted_hit_count})) {
                $ftr->{quality}->{weighted_hit_count} = $ftr->{quality}->{weighted_hit_count}+0;
            }
            if (defined($ftr->{quality}->{hit_count})) {
                $ftr->{quality}->{hit_count} = $ftr->{quality}->{hit_count}+0;
            }
            if (defined($ftr->{annotations})) {
                delete $ftr->{annotations};
            }
            if (defined($ftr->{location})) {
                $ftr->{location}->[0]->[1] = $ftr->{location}->[0]->[1]+0;
                $ftr->{location}->[0]->[3] = $ftr->{location}->[0]->[3]+0;
            }
            delete $ftr->{feature_creation_event};
            if (defined($ftr->{function}) && length($ftr->{function}) > 0 && defined($ftr->{protein_translation})) {
                for (my $j=0; $j < @{$genome->{features}}; $j++) {
                    if ($i ne $j) {
                        if ($ftr->{location}->[0]->[0] eq $genome->{features}->[$j]->{location}->[0]->[0]) {
                            if ($ftr->{location}->[0]->[1] == $genome->{features}->[$j]->{location}->[0]->[1]) {
                                if ($ftr->{location}->[0]->[2] eq $genome->{features}->[$j]->{location}->[0]->[2]) {
                                    if ($ftr->{location}->[0]->[3] == $genome->{features}->[$j]->{location}->[0]->[3]) {
                                        if ($ftr->{type} ne $genome->{features}->[$j]->{type}) {
                                            $genome->{features}->[$j]->{function} = $ftr->{function};
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        for (my $i=0; $i < @{$genome->{features}}; $i++) {
            my $ftr = $genome->{features}->[$i];
            if (defined($oldfunchash->{$ftr->{id}}) && (!defined($ftr->{function}) || $ftr->{function} =~ /hypothetical\sprotein/)) {
                if (defined($parameters->{retain_old_anno_for_hypotheticals}) && $parameters->{retain_old_anno_for_hypotheticals} == 1) {
                    $ftr->{function} = $oldfunchash->{$ftr->{id}};
                }
            }
            if (defined($ftr->{function}) && length($ftr->{function}) > 0) {
                my $function = $ftr->{function};
                my $array = [split(/\#/,$function)];
                $function = shift(@{$array});
                $function =~ s/\s+$//;
                $array = [split(/\s*;\s+|\s+[\@\/]\s+/,$function)];
                my $marked = 0;
                for (my $j=0;$j < @{$array}; $j++) {
                    my $rolename = lc($array->[$j]);
                    $rolename =~ s/[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+//g;
                    $rolename =~ s/\s//g;
                    $rolename =~ s/\#.*$//g;
                    $genomefunchash->{$rolename} = 1;
                    if (defined($funchash->{$rolename})) {
                        $seedfunchash->{$rolename} = 1;
                        if ($marked == 0) {
                            $seedfunctions++;
                            $marked = 1;
                        }
                        my $ont_term = $ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}};
                        if (!defined($ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}})) {
                            $ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}} = {
                                 evidence => [],
                                 id => $funchash->{$rolename}->{id},
                                 term_name => $funchash->{$rolename}->{name},
                                 ontology_ref => "KBaseOntology/seed_subsystem_ontology",
                                 term_lineage => [],
                            };
                        }
                        my $found = 0;
                        if (ref($ont_term) eq 'ARRAY'){
                            push(@{$ont_term}, $#{$inputgenome->{ontology_events}});
                        } else {
                            for (my $k = 0; $k < @{$ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}}; $k++) {
                                if ($ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}->[$k]->{method} eq Bio::KBase::Utilities::method()) {
                                    $ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}->[$k]->{timestamp} = Bio::KBase::Utilities::timestamp();
                                    $ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}->[$k]->{method_version} = $self->_util_version();
                                    $found = 1;
                                    last;
                                }
                            }
                            if ($found == 0) {
                                push(
                                    @{$ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}},
                                    {
                                        method         => Bio::KBase::Utilities::method(),
                                        method_version => $self->_util_version(),
                                        timestamp      => Bio::KBase::Utilities::timestamp()
                                    });
                            }
                        }
                        if (exists ($ftr->{ontology_terms}->{SSO})) {
                            for my $sso (keys($ftr->{ontology_terms}->{SSO})) {
                                if ($sso =~ /SSO:000009304/) {
                                    $advancedmessage .= "Found selenocysteine-containing gene $ftr->{ontology_terms}->{SSO}->{$sso}\n";
                                }
                                elsif ($sso =~ /SSO:000009291/) {
                                    $advancedmessage .= "Found pyrrolysine-containing gene $ftr->{ontology_terms}->{SSO}->{$sso}\n";
                                }
                            }
                        }
                    }
                }
            }
        }
        ## Rolling protein features back from 'CDS' to 'gene':
        for (my $i=0; $i < @{$genome->{features}}; $i++) {
            my $ftr = $genome->{features}->[$i];
            if ($ftr->{type} eq "CDS") {
                $ftr->{type} = "gene"
            }
            if (exists $oldtype->{$ftr->{id}} && $oldtype->{$ftr->{id}} =~ /gene/) {
                $ftr->{type} = $oldtype->{$ftr->{id}};
            }
            if (defined($ftr->{location})) {
                $ftr->{location}->[0]->[1] = $ftr->{location}->[0]->[1]+0;
                $ftr->{location}->[0]->[3] = $ftr->{location}->[0]->[3]+0;
            }
#           $ftr->{type} = 'gene';
            my $type = '';
            if ( defined$ftr->{type}) {
                $type = 'Coding '.$ftr->{type};
            } else {
                $type = 'Coding ';
            }
            #   Count the output feature types
            if (exists $types{$type}) {
                $types{$type} += 1;
            } else {
                $types{$type} = 1;
            }
            delete $genome->{features}->[$i]->{type} if (exists $ftr->{type});
            delete $genome->{features}->[$i]->{protein_md5} if (exists $ftr->{protein_md5});
        }
        if ((not defined($genome->{cdss})) || (not defined($genome->{mrnas})) || $update_cdss eq 'Y') {
            ## Reconstructing new feature arrays ('cdss' and 'mrnas') if they are not present:
            my $cdss = [];
            $genome->{cdss} = $cdss;
            my $mrnas = [];
            $genome->{mrnas} = $mrnas;
            for (my $i=0; $i < @{$genome->{features}}; $i++) {
                my $feature = $genome->{features}->[$i];
                if (defined($feature->{protein_translation})) {
                    my $gene_id = $feature->{id};
                    my $cds_id = $gene_id . "_CDS";
                    my $mrna_id = $gene_id . "_mRNA";
                    my $location = $feature->{location};
                    my $md5 = $feature->{md5};
                    my $function = $feature->{function};
                    if (not defined($function)) {
                        $function = "";
                    }
                    my $ontology_terms = $feature->{ontology_terms};
                    if (not defined($ontology_terms)) {
                        $ontology_terms = {};
                    }
                    my $protein_translation = $feature->{protein_translation};
                    my $protein_translation_length = length($protein_translation);
                    my $aliases = [];
                    if (defined($feature->{aliases})) {
                        $aliases = $feature->{aliases};
                        # Found some pseudogenes that have wrong structure for aliases
                        if (ref($aliases) !~ /ARRAY/) {
                            $aliases = [$feature->{aliases}];
                        }
                    }
                    push(@{$cdss}, {
                        id => $cds_id,
                        location => $location,
                        md5 => $md5,
                        parent_gene => $gene_id,
                        parent_mrna => $mrna_id,
                        function => $function,
                        ontology_terms => $ontology_terms,
                        protein_translation => $protein_translation,
                        protein_translation_length => $protein_translation_length,
                        aliases => $aliases
                    });
                    push(@{$mrnas}, {
                        id => $mrna_id,
                        location => $location,
                        md5 => $md5,
                        parent_gene => $gene_id,
                        cds => $cds_id
                    });
                    $feature->{cdss} = [$cds_id];
                    $feature->{mrnas} = [$mrna_id];
                }
            }
        }
    }

    $rast_details{types} = \%types;
    $rast_details{num_coding} = $num_coding;
    $rast_details{newncfs} = $newncfs;
    $rast_details{newftrs} = $newftrs;
    $rast_details{genomefunchash} = $genomefunchash;
    $rast_details{seedfunctions} = $seedfunctions;
    $rast_details{seedfunchash} = $seedfunchash;

    return ($genome, \%rast_details);
}

sub _summarize_annotation {
    my ($self, $rast_ref, $genome, $inputgenome) = @_;

    print "**************In _summarize_annotation*************\n";
    my %rast_details = %{ $rast_ref };
    my $message = $rast_details{message};
    my $contigobj = $rast_details{contigobj};
    my %types = ();
    if (exists($rast_details{types})) {
        %types = %{$rast_details{types}};
    }
    my $num_coding = $rast_details{num_coding};
    my $genehash = $rast_details{genehash};
    my $newncfs = $rast_details{newncfs};
    my $newftrs = $rast_details{newftrs};
    my $genomefunchash = $rast_details{genomefunchash};
    my $seedfunctions = $rast_details{seedfunctions};
    my $seedfunchash = $rast_details{seedfunchash};

    my $num_non_coding = 0;
    if (defined($genome->{non_coding_features})) {
        for (my $i=0; $i < @{$genome->{non_coding_features}}; $i++) {
            my $ftr = $genome->{non_coding_features}->[$i];
            if (defined($genehash) && !defined($genehash->{$ftr->{id}})) {
                # Let's count number of non_coding_features with functions updated by RAST service.
                # If function is not set we treat it as empty string to avoid perl warning.
                $newncfs++;
                $newftrs++;
                my $func = $ftr->{function};
                if (not defined($func)) {
                    $func = "";
                }

            } else {
                $num_non_coding++;
            }
            my $type = '';
            if ( defined$ftr->{type} && $ftr->{type} =~ /coding/) {
                $type = $ftr->{type} ;
            } elsif (defined$ftr->{type}) {
                $type = 'Non-coding '.$ftr->{type} ;
            } else {
                $type = 'Non-coding ';
            }
            #   Count the output feature types
            if (exists $types{$type}) {
                $types{$type} += 1;
            } else {
                $types{$type} = 1;
            }
        }
    }
    if (defined($inputgenome)) {
        if (defined($num_coding) && defined($num_non_coding) && defined($newftrs) && defined($newncfs)) {
        $message .= "In addition to the remaining original $num_coding coding features and $num_non_coding non-coding features, ".$newftrs." new features were called, of which $newncfs are non-coding.\n";
        }
        if (%types) {
            $message .= "Output genome has the following feature types:\n";
            for my $key (sort keys(%types)) {
                $message .= sprintf("\t%-30s %5d \n", $key, $types{$key});
            }
        }
    }

    if (defined($genomefunchash)) {
        $message .= "Overall, the genes have ".keys(%{$genomefunchash})." distinct functions";
    }
    if (defined($seedfunctions) && defined($seedfunchash)) {
        $message .= "\nThe genes include ".$seedfunctions." genes with a SEED annotation ontology across ".keys(%{$seedfunchash})." distinct SEED functions.\n";
    }
    $message .= "The number of distinct functions can exceed the number of genes because some genes have multiple functions.\n";
    print($message);

    if (!defined($genome->{assembly_ref})) {
        delete $genome->{assembly_ref};
    }
    if (defined($contigobj)) {
        $genome->{gc_content} = $contigobj->{_gc};
        if ($contigobj->{_kbasetype} eq "ContigSet") {
            $genome->{contigset_ref} = $contigobj->{_reference};
        } else {
            $genome->{assembly_ref} = $contigobj->{_reference};
        }
    }
    $rast_details{message} = $message;
    $rast_details{types} = \%types;
    return ($genome, \%rast_details);
}

sub _save_annotation_results {
    my ($self, $genome, $rast_ref) = @_;

    print "SEND OFF FOR SAVING\n";
#   print "***** Domain       = $genome->{domain}\n";
#   print "***** Genitic_code = $genome->{genetic_code}\n";
#   print "***** Scientific_namee = $genome->{scientific_name}\n";
#   print "***** Number of features=".scalar  @{$genome->{features}}."\n";
#   print "***** Number of non_coding_features=".scalar  @{$genome->{non_coding_features}}."\n";
#   print "***** Number of cdss=    ".scalar  @{$genome->{cdss}}."\n";
#   print "***** Number of mrnas=   ".scalar  @{$genome->{mrnas}}."\n";

    my %gc_rast = %{ $rast_ref };
    my $parameters = $gc_rast{parameters},
    my $message = $gc_rast{message};
    my $rasted_gn = $genome;

    my $g_data_type = (ref($rasted_gn) eq 'HASH') ? 'ref2Hash' : ref($rasted_gn);
    if ($g_data_type eq 'GenomeTypeObject') {
        print "INFO***Genome input passed to _save_annotation_results is of type of $g_data_type, prepare it before saving***.\n";
        $rasted_gn = $rasted_gn->prepare_for_return();
    }

    $self->_check_NC_features($rasted_gn);
    if( $rasted_gn->{contig_ids} ) {
        my $remapped_ctg_ids = $self->_remap_contigIDs($gc_rast{contigID_hash},
                                                       $rasted_gn->{contig_ids});
        $rasted_gn->{contig_ids} = $remapped_ctg_ids;
    }

    my $gfu_client = installed_clients::GenomeFileUtilClient->new($self->{call_back_url});
    my ($gaout, $gaout_info);

    try {
        $gaout = $gfu_client->save_one_genome({
            workspace => $parameters->{output_workspace},
            name => $parameters->{output_genome_name},
            data => $rasted_gn,
            provenance => [{
                "time" => DateTime->now()->datetime()."+0000",
                service_ver => $self->_util_version(),
                service => "RAST_SDK",
                method => Bio::KBase::Utilities::method(),
                method_params => [$parameters],
                input_ws_objects => [],
                resolved_ws_objects => [],
                intermediate_incoming => [],
                intermediate_outgoing => []
            }],
            hidden => 0
        });

        $gaout_info = $gaout->{info};
        my $ref = $gaout_info->[6]."/".$gaout_info->[0]."/".$gaout_info->[4];
        Bio::KBase::KBaseEnv::add_object_created({
            "ref" => $ref,
            "description" => "Annotated genome"
        });
        Bio::KBase::Utilities::print_report_message({
            message => "<pre>".$message."</pre>",
            append => 0,
            html => 0
        });
        print "One genome has been saved\n";
        return ({"ref" => $ref}, $message);
    } catch {
        my $err_msg = "ERROR: Calling GFU.save_one_genome failed with error message:$_\n";
        return ({}, $err_msg);
    };
}

## _build_workflows
#  input: $parameters
#
## return:(
#  {genecall_workflow=>$workflow1,
#   annotation_workflow=>$workflow2,
#   renumber_workflow=>$workflow3
#   oldfunchash=>$oldfunchash,
#   oldtype=>$oldtype,
#   message=>$message,
#   types=>%types,
#   contigobj=>$contigobj,
#   parameters=>$parameters,
#   ...
#  }, $inputgenome)
#
sub _build_workflows {
    my ($self, $parameters) = @_;

    ## refactor 1 -- set the parameter values from $parameters and initiate $inputgenome
    # 1. creating default genome object
    my $inputgenome = {
        id => $parameters->{output_genome_name},
        genetic_code => $parameters->{genetic_code},
        scientific_name => $parameters->{scientific_name},
        domain => $parameters->{domain},
        contigs => [],
        features => []
    };
    if ($parameters->{ncbi_taxon_id}) {
        $inputgenome->{taxon_assignments} = {'ncbi' => '' . $parameters->{ncbi_taxon_id}};
    }

    my (%rast_details, $rast_ref);
    ($rast_ref, $inputgenome) = $self->_set_parameters_by_input(
                                    $parameters, $inputgenome);

    # 2. merge with the default gene call settings
    %rast_details = %{ $rast_ref };
    $parameters = $rast_details{parameters};
    my $default_params = $self->_set_default_parameters();
    $parameters = { %$default_params, %$parameters };
    $rast_details{parameters} = $parameters;

    ## refactor 2 -- taking notes in $message
    ($rast_ref, $inputgenome) = $self->_set_messageNcontigs(
                                    \%rast_details, $inputgenome);

    ## refactor 3 -- set gene call workflow
    $rast_ref = $self->_set_genecall_workflow($rast_ref, $inputgenome);

    ## refactor 4 -- set annotation workflow
    $rast_ref = $self->_set_annotation_workflow($rast_ref);

    ## refactor 5 -- combine messages and renumber features
    ($rast_ref, $inputgenome) = $self->_renumber_features($rast_ref, $inputgenome);

    ## refactor 6 -- pre_rast_call
    ($rast_ref, $inputgenome) = $self->_pre_rast_call($rast_ref, $inputgenome);

    return ($rast_ref, $inputgenome);
}

#
# Makes gene calls on the input genome with the gene call stages defined in $wf_gcs.
# returns a genes called new genome.
#
sub _run_rast_genecalls {
    my ($self, $in_genome, $wf_gcs) = @_;

    return $self->_run_rast_workflow($in_genome, $wf_gcs);
}

##----end gene call subs----##

sub _get_fasta_from_assembly {
    my ($self, $assembly_ref) = @_;

    my $au = installed_clients::AssemblyUtilClient->new($self->{call_back_url});
    my $outfa = {};
    try {
        $outfa = $au->get_assembly_as_fasta({"ref" => $assembly_ref});
        return $outfa->{path};
    } catch {
        croak "ERROR calling AssemblyUtil.get_assembly_as_fasta:\n$_\n";
    };
}

sub _write_fasta_from_genome {
    my ($self, $input_obj_ref) = @_;

    my $fa_file = '';
    try {
        my $genome_obj = $self->_fetch_object_data($input_obj_ref);
        if ($genome_obj->{assembly_ref}) {
            print "*******Input genome's assembly ref is: $genome_obj->{assembly_ref}**********\n";

            $fa_file = $self->_get_fasta_from_assembly(
                          $input_obj_ref.";".$genome_obj->{assembly_ref});
            unless (-e $fa_file && -s $fa_file) {print "Fasta file is empty!!!!";}
        }
        else {
            warn ("**_write_fasta_from_genome ERROR:\n".
            "No assembly can be found for genome $input_obj_ref\n");
        }
    } catch {
        warn "**_write_fasta_from_genome ERROR: ".$_."\n";
    } finally {
        return $fa_file;
    };
}

sub _write_gff_from_genome {
    my ($self, $genome_ref) = @_;

    my $chk = $self->_validate_KB_objref_name($genome_ref);
    return '' unless $chk->{check_passed};

    my $obj_info = $self->_fetch_object_info($genome_ref, $chk);
    return '' unless $obj_info;

    $genome_ref = $obj_info->[6].'/'.$obj_info->[0].'/'.$obj_info->[4];
    my $in_type = $obj_info->[2];
    my $is_assembly = ($in_type =~ /KBaseGenomeAnnotations\.Assembly/ ||
                       $in_type =~ /KBaseGenomes\.ContigSet/);
    my $is_genome = ($in_type =~ /KBaseGenomes\.Genome/ ||
                     $in_type =~ /KBaseGenomeAnnotations\.GenomeAnnotation/);

    unless ($is_genome) {
        croak "ValueError: Object is not an KBaseAnnotations.GenomeAnnotation or Genome, will throw an error.\n";
    }
    my $gfu = installed_clients::GenomeFileUtilClient->new($self->{call_back_url});
    my $gff_result;
    try {
        $gff_result = $gfu->genome_to_gff({"genome_ref" => $genome_ref});

        unless (-s $gff_result->{file_path}) {print "GFF is empty ";}
        return $gff_result->{file_path};
    } catch {
        croak "**_write_gff_from_genome ERROR:\n$_\n";
    };
}

sub _validate_KB_objref_name {
    my ($self, $obj_str) = @_;
    my $ret = {
            "is_ref" => 0,
            "is_name" => 0,
            "check_passed" => 1
       };

    my @str_arr = split ('/', $obj_str);

    if( @str_arr > 3 ) {
        $ret->{check_passed} = 0;
        return $ret;
    }

    if( @str_arr == 1 ) {
        ## assuming it is a string as an object name
        if( $str_arr[0] =~ m/[^\\w\\|._-]/ ) {
            $ret->{is_name} = 1;
            return $ret;
        }
    }

    unless( $str_arr[0] =~ m/[^\\w\\|._-]/ ) {
        $ret->{check_passed} = 0;
        return $ret;
    }

    unless( $str_arr[1] =~ m/[^\\w\\|._-]/ ) {
        $ret->{check_passed} = 0;
        return $ret;
    }

    if( exists($str_arr[2]) ) {
        unless( $str_arr[2] =~ m/^\d+$/ ) {
            $ret->{check_passed} = 0;
            return $ret;
        }
    }
    $ret->{is_ref} = 1;
    return $ret;
}

sub _check_annotation_params {
    my ($self, $params) = @_;

    my $missing_params = "Missing required parameters for annotating genome.\n";
    unless (defined($params)) {
        print "params is not defined!!!!\n";
        croak $missing_params;
    }

    if (!keys %$params) {
        print "params is empty!!!!\n";
        croak $missing_params;
    }
    my $req1 = "'output_workspace' is required for running rast_genome.\n";
    my $invald1 = "Invalid workspace name:";
    if (!defined($params->{output_workspace}) || $params->{output_workspace} eq '') {
        croak $req1;
    }
    elsif ($params->{output_workspace} !~ m/[^\\w:._-]/) {
        croak $invald1.$params->{output_workspace}.'\n';
    }

    my $req2 = "'object_ref' is required for running rast_genome.\n";
    if (!defined($params->{object_ref}) || $params->{object_ref} eq '') {
        croak $req2;
    }
    if (!defined($params->{output_genome_name})
        || $params->{output_genome_name} eq '') {
        $params->{output_genome_name} = "rast_annotated_genome";
    }
    if (!defined($params->{scientific_name})
        || $params->{scientific_name} eq '') {
        $params->{scientific_name} = "Unknown species";
    }
    if (!defined($params->{genetic_code})) {
        $params->{genetic_code} = 11;
    }
    if (!defined($params->{domain})
        || $params->{domain} eq '') {
        $params->{domain} = "Bacteria";
    }
    if (!defined($params->{create_report})) {
        $params->{create_report} = 0;
    }
    return $params;
}

sub _check_bulk_annotation_params {
    my ($self, $params) = @_;

    my $missing_params = "ERROR: Missing required parameters for annotating genomes/assemblies.\n";
    unless (defined($params)) {
        print "params is not defined!!!!\n";
        croak $missing_params;
    }

    if (!keys %$params) {
        print "params is empty!!!!\n";
        croak $missing_params;
    }
    my $req1 = "'output_workspace' is required.\n";
    my $invald1 = "Invalid workspace name:";
    if (!defined($params->{output_workspace}) || $params->{output_workspace} eq '') {
        croak $req1;
    }
    elsif ($params->{output_workspace} !~ m/[^\\w:._-]/) {
        croak $invald1.$params->{output_workspace}.'\n';
    }
    if (!defined($params->{output_GenomeSet_name})
        || $params->{output_GenomeSet_name} eq '') {
        $params->{output_GenomeSet_name} = "rasted_GenomeSet_name";
    }
    if (!defined($params->{output_AMASet_name})
        || $params->{output_AMASet_name} eq '') {
        $params->{output_AMASet_name} = "rasted_AMASet_name";
    }
    if (!defined($params->{create_report})) {
        $params->{create_report} = 0;
    }
    if (!defined($params->{input_genomes})) {
        $params->{input_genomes} = [];
    }
    elsif (ref $params->{input_genomes} ne 'ARRAY') {
        $params->{input_genomes} = [$params->{input_genomes}];
    }
    $params->{input_genomes} = $self->_uniq_ref($params->{input_genomes});
    if (!defined($params->{input_AMAs})) {
        $params->{input_AMAs} = [];
    }
    elsif (ref $params->{input_AMAs} ne 'ARRAY') {
        $params->{input_AMAs} = [$params->{input_AMAs}];
    }
    $params->{input_AMAs} = $self->_uniq_ref($params->{input_AMAs});
    if (!defined($params->{input_assemblies})) {
        $params->{input_assemblies} = [];
    }
    elsif (ref $params->{input_assemblies} ne 'ARRAY') {
        $params->{input_assemblies} = [$params->{input_assemblies}];
    }
    $params->{input_assemblies} = $self->_uniq_ref($params->{input_assemblies});
    unless ($params->{input_text}) {
        $params->{input_text} = '';
    }
    if (!defined($params->{genetic_code})) {
        $params->{genetic_code} = 11;
    }
    if (!defined($params->{domain})
        || $params->{domain} eq '') {
        $params->{domain} = "Bacteria";
    }
    return $params;
}


sub _generate_stats_from_aa {
    my ($self, $gn_ref) = @_;

    my %gn_stats = ();
    my $chk = $self->_validate_KB_objref_name($gn_ref);
    return %gn_stats unless $chk->{check_passed};

    print "++++++++++++++_generate_stats_from_aa on $gn_ref++++++++++++\n";
    my $gn_info = $self->_fetch_object_info($gn_ref, $chk);
    return %gn_stats unless $gn_info;

    $gn_ref = $gn_info->[6].'/'.$gn_info->[0].'/'.$gn_info->[4];
    my $in_type = $gn_info->[2];
    my $is_assembly = ($in_type =~ /KBaseGenomeAnnotations\.Assembly/ ||
                       $in_type =~ /KBaseGenomes\.ContigSet/);
    my $is_genome = ($in_type =~ /KBaseGenomes\.Genome/ ||
                     $in_type =~ /KBaseGenomeAnnotations\.GenomeAnnotation/);
    my $is_meta_assembly = $in_type =~ m/Metagenomes\.AnnotatedMetagenomeAssembly/;

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
    elsif ($is_meta_assembly || $is_genome) {
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
    for my $gff_line (@$gff_contents) {
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
            for my $attr (split(";",$ftr_attrs)) {
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
    for my $ftr_id (sort keys %$func_tab_ref) {
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

    my $report_title = "Feature function report for genome <font color=green>$obj_stats{id}</font>";
    my $rpt_header = "<h3>$report_title:</h3>";
    my $rpt_data = ("data.addColumn('string', 'function role');\n".
                    #"data.addColumn('string', 'annotation source');\n".
                    "data.addColumn('number', 'gene count');\n".
                    "data.addColumn('string', 'subsystem name');\n".
                    "data.addColumn('string', 'subsystem class');\n".
                    "data.addRows([\n");

    my $roles = $gff_stats{function_roles};
    my $genes = $gff_stats{gene_role_map};
    for my $role_k (sort keys %$roles) {
        # add escape to preserve double quote (")
        (my $new_role_k = $role_k) =~ s/"/\\"/g;
        $new_role_k = uri_decode($new_role_k);
        $rpt_data .= '["<span style=\"white-space:nowrap;\">'."$new_role_k</span>\",";
        #my $ann_src = $self->_find_function_source($func_tab_ref, $role_k);
        #$rpt_data .= "\"$ann_src\",";
        $rpt_data .= "$roles->{$role_k}->{gene_count},";
        #$rpt_data .= "\"$roles->{$role_k}->{gene_list}\"],\n";

        # search for $role_k in all the $subsys_info items's role_names arrays
        my $subsys_str = '';
        my $class_str = '';
        for my $subsys_k (sort keys %subsys_info) {
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
    my $rpt_footer .= "<p><strong>Gene Total Count = $gene_count</strong></p>\n";
    if (defined($obj_stats{contig_count}) && defined($obj_stats{num_features})) {
        $rpt_footer .= "<p><strong>Contig Count = $obj_stats{contig_count}</strong></p>\n";
        $rpt_footer .= "<p><strong>Feature Count = $obj_stats{num_features}</strong></p>\n";
        # $rpt_footer .= "<p><strong>GC Content = $obj_stats{gc_content}</strong></p>\n";
    }
    my $srch1 = "(<replaceHeader>)(.*)(</replaceHeader>)";
    my $srch2 = "(<replaceFooter>)(.*)(</replaceFooter>)";
    my $srch3 = "(<replaceTableData>)(.*)(</replaceTableData>)";

    $file_content =~ s/$srch1/$rpt_header/;
    #$file_content =~ s/$srch2/$rpt_footer/;
    $file_content =~ s/$srch3/$rpt_data/;

    my $report_file_path = catfile($self->{genome_dir}, 'genome_report.html');
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

#Create a KBaseReport with brief info/stats on a reannotated genome
sub _generate_genome_report {
    my ($self, $aa_ref, $aa_gff_conts, $func_tab, $ftr_cnt, $msg) = @_;

    my $chk = $self->_validate_KB_objref_name($aa_ref);
    unless( $chk->{check_passed} ) {
		return {"output_genome_ref"=>$aa_ref,
                "workspace_name"=>undef,
                "report_name"=>undef,
                "report_ref"=>undef};
    }

    my $gn_info = $self->_fetch_object_info($aa_ref, $chk);
    unless( $gn_info ) {
        return {"output_genome_ref"=>$aa_ref,
                "workspace_name"=>undef,
                "report_name"=>undef,
                "report_ref"=>undef};
    }
    my $gn_ws = $gn_info->[7];
    $aa_ref = $gn_info->[6].'/'.$gn_info->[0].'/'.$gn_info->[4];
    my %aa_stats = $self->_generate_stats_from_aa($aa_ref);
    my %aa_gff_stats = $self->_generate_stats_from_gffContents($aa_gff_conts);
    my %subsys_info = $self->_fetch_subsystem_info();

    my @html_files = $self->_write_html_from_stats(\%aa_stats,
                                                   \%aa_gff_stats,
                                                   \%subsys_info,
                                                   $func_tab);

    my $kbr = installed_clients::KBaseReportClient->new($self->{call_back_url});
    my $report_message = $msg;
    my $report_info = $kbr->create_extended_report(
            {"message"=>$report_message,
             "objects_created"=>[{"ref"=>$aa_ref, "description"=>"RAST re-annotated genome"}],
             "html_links"=> \@html_files,
             "direct_html_link_index"=> 0,
             "html_window_height"=> 366,
             "report_object_name"=>"kb_RAST_genome_report_".$self->_create_uuid(),
             "workspace_name"=>$gn_info->[7]
        });

    return {"output_genome_ref"=>$aa_ref,
            "workspace_name"=>$report_info->{workspace_name},
            "report_name"=>$report_info->{name},
            "report_ref"=>$report_info->{ref}};
}

#
## Check the given ref string to make sure it conforms to the expected
## object reference format and actually exists in KBase.
## Return the object data if it passes, die otherwise.
#
sub _fetch_object_data {
    my ($self, $obj_ref) = @_;
    my $ret_obj_data = {};

    unless ($obj_ref =~ m/^\d+\/\d+\/\d+$/) {
        print "invalid object reference: $obj_ref\n";
        return {};
    }
    try {
        $ret_obj_data = $self->{ws_client}->get_objects2(
                            {'objects'=>[{ref=>$obj_ref}]}
                        )->{data}->[0]->{data};
    } catch {
        warn "ERROR Workspace.get_objects2 failed to access $obj_ref:\n$_\n";
    } finally {
        return $ret_obj_data;
    };
}

#
## Check the given ref string to make sure it conforms to the expected
## object reference format and actually exists in KBase.
## Return the object info if it passes, undef otherwise.
#
sub _fetch_object_info {
    my ($self, $obj, $chk, $ws) = @_;

    my $ret_obj_info = undef;
    my $objs = [];
    if( $chk->{is_ref} ) {
        $objs = [{ref=>$obj}];
    } elsif( $chk->{is_name} && defined($ws) ) {
        my $ref = Bio::KBase::KBaseEnv::buildref($ws, $obj);
        $objs = [{ref=>$ref}];
    }

    try {
        $ret_obj_info = $self->{ws_client}->get_object_info3(
                                 {objects=>$objs}
                        )->{infos}->[0];
        print "INFO: object info for $obj------\n".Dumper($ret_obj_info);
        return $ret_obj_info;
    } catch {
        warn "INFO: Workspace.get_object_info3 failed to access $obj.\n";
        warn "ERROR message:$_\n";
        return undef;
    };
}


# create a unique ID
sub _create_uuid {
    return Data::UUID->new()->create_str();
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

# create the project directory for rast annotation
sub _create_rast_subdir {
    my ($self, $rast_dir, $sub_dir) = @_;
    my $subdir = $sub_dir.$self->_create_uuid();
    my $dir = catfile($rast_dir, $subdir);

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

sub _get_file_lines {
    my ($self, $filename) = @_;
    # Open $filename to read into an array
    my $fh = $self->_openRead($filename);
    my @file_lines=();
    chomp(@file_lines = <$fh>);
    close($fh);
    return \@file_lines;
}


##----subs for parsing GFF by Seaver----##
sub _parse_gff {
    my ($self, $gff_filename, $attr_delimiter) = @_;
    print "Parsing GFF contents from file $gff_filename \n";

    my $gff_lines = $self->_get_file_lines($gff_filename);
    print "Read in ".scalar @{$gff_lines}." lines from GFF file $gff_filename\n";

    my @gff_contents=();
    for my $current_line (@{$gff_lines}){
        chomp $current_line;
        next if (!$current_line || $current_line =~ m/^#.*$/);
        next if $current_line =~ m/^##gff-version\s*\d$/;

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
        for my $attribute (split(";",$attributes)) {
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
            #Added the trim/chomp function to avoid cases like extra space or newline:
            # ' parent' => 'mRNA_2',
            # ' product' => 'malate/citrate symporter
            #'
            # 'id' => 'EPWB_RS00005
            #'
            #
            trim $key;
            trim $value;
            chomp $value;

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
    #if ($ftr_count > 10) {
    #    print "INFO: First 10 RAST feature examples:\n".Dumper(@{$features}[0..9]);
    #}
    #else {
    #    print "INFO:All $ftr_count RAST features:\n".Dumper(@{$features});
    #}

    #Feature Lookup Hash
    my %function_lookup = ();
    for my $ftr (@$features){
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

## Combine all the workflows into one
sub _combine_workflows {
    my ($self, $rast_ref) = @_;

    my %rast_details = %{ $rast_ref };
    my $comp_workflow = {stages => []};
    my $wf_genecall = $rast_details{genecall_workflow};
    my $wf_annotate = $rast_details{annotate_workflow};
    my $wf_renum = $rast_details{renumber_workflow};

    if ($wf_genecall && @{$wf_genecall->{stages}}) {
        # print "There are genecall workflows:\n".Dumper($wf_genecall);
        push @{$comp_workflow->{stages}}, @{$wf_genecall->{stages}};
    }
    if ($wf_annotate && @{$wf_annotate->{stages}}) {
        # print "There are annotation workflows:\n".Dumper($wf_annotate);
        push @{$comp_workflow->{stages}}, @{$wf_annotate->{stages}};
    }
    ## could be skipped: running the renumber_features workflow
    if ($wf_renum && @{$wf_renum->{stages}}) {
        # print "There is a renumber feature workflow:\n".Dumper($wf_renum);
        push @{$comp_workflow->{stages}}, @{$wf_renum->{stages}};
    }

    return $comp_workflow;
}


##----main function----##
#
# Runs through the workflows on the input genome with genecall and annotation workflow
# stages defined in $rast_ref.
# returns a RAST-ed new genome, the original input genome together with the updated
# rast-ing details.
#
sub rast_genome {
    my $self = shift;
    my($inparams) = @_;

    print "rast_genome input parameter=\n". Dumper($inparams). "\n";
    ## 0. Doing nothing on an invalid input
    my $chk = $self->_validate_KB_objref_name($inparams->{object_ref});
    unless( $chk->{check_passed} ) {
        print "INFO: Invalid workspace object reference: $inparams->{object_ref}";
        return {};
    }

    my $params = $self->_check_annotation_params($inparams);
    my $input_obj_ref = $params->{object_ref};

    ## 1. build the workflows for RAST-ing
    my ($rast_ref, $inputgenome) = $self->_build_workflows($params);
    my %rast_details = %{ $rast_ref };

    ## 2. combine all workflows and run it at once
    my $rasted_genome = $self->_run_rast_workflow(
                                    $inputgenome,
                                    $self->_combine_workflows($rast_ref));
    my $rast_ftrs = $rasted_genome->{features};
    my $rast_ftr_count = scalar @{$rast_ftrs};
    unless ($rast_ftr_count > 0) {
        print( "Empty input genome features after full workflow, "
               ."return an empty object.\n" );
        return {};
    }

    print "\n***********RAST calling resulted in ".$rast_ftr_count." features.\n";

    my $genome_final = $self->_post_rast_ann_call($rasted_genome,
                                                  $inputgenome,
                                                  $rast_details{parameters},
                                                  $rast_details{contigobj});

    my $cd_ftr_count = @{$genome_final->{features}};
    my $nc_ftr_count = @{$genome_final->{non_coding_features}};
    print( "\n***********Post rast processing resulted in $cd_ftr_count coding features "
           ."and $nc_ftr_count non_coding features.\n");

    ## 3. build seed ontology
    ($genome_final, $rast_ref) = $self->_build_seed_ontology(
                                         $rast_ref, $genome_final, $inputgenome);

    ## 4. summarize annotation
    ($genome_final, $rast_ref) = $self->_summarize_annotation(
                                         $rast_ref, $genome_final, $inputgenome);

    my $ftrs = $genome_final->{features};
    my $rasted_ftr_count = scalar @{$ftrs};
    print "\n***********Finally RAST resulted ".$rasted_ftr_count." features.\n";

    my $prnt_num = 10;
    if ($rasted_ftr_count) {
        print "***********The first $prnt_num or fewer rasted features, for example***************\n";
        my $prnt_lines = ($rasted_ftr_count > $prnt_num) ? $prnt_num : $rasted_ftr_count;
        for my $ftr (@{$ftrs}[0, $prnt_lines - 1]) {
            my $f_id = $ftr->{id};
            # if $ftr->{ function } is defined, set $f_func to it; otherwise, set it to ''
            my $f_func = $ftr->{ function } // '';
            my $f_protein = $ftr->{protein_translation} // '';
            print "$f_id\t$f_func\t$f_protein\n";
        }
    }

    ## 5. save the annotated genome
    ## If saving failed, return {} so that downstream can skip it
    my ($aa_out, $out_msg) = $self->_save_annotation_results($genome_final, $rast_ref);
    unless ( %$aa_out ) {
        print( $out_msg."No annotation is saved on $input_obj_ref\n" );
        return {};
    }

    print "_save_annotation_results returned:\n".Dumper($aa_out);
    my $aa_ref = $aa_out->{ref};
    my $rast_ret = {
        output_genome_ref => $aa_ref,
        output_workspace => $params->{output_workspace},
        report_name => undef,
        report_ref => undef
    };

    my %ftr_func_lookup = $self->_get_feature_function_lookup($ftrs);
    my $gff_contents = $self->_get_genome_gff_contents($aa_ref);
    if (defined($aa_ref) && defined($params->{create_report}) &&
            $params->{create_report} == 1) {
        $rast_ret = $self->_generate_genome_report(
                          $aa_ref, $gff_contents, \%ftr_func_lookup,
                          $rasted_ftr_count, $out_msg);
    }
    return $rast_ret;
}

#
## Helper function 1 for bulk_rast_genomes
sub _build_param_from_obj {
    my ($self, $obj, $chk, $ws, $out_genomeSet) = @_;

    my $obj_info = $self->_fetch_object_info($obj, $chk, $ws);
    return undef unless $obj_info;

    $obj = $obj_info->[6].'/'.$obj_info->[0].'/'.$obj_info->[4];
    my $obj_name = $obj_info->[1];
    return {
        object_ref => $obj,
        output_workspace => $ws,
        output_genome_name => $out_genomeSet . '_' .$obj_name,
        create_report => 0
    };
}

#
## Helper function 3 for bulk_rast_genomes
## Parse the inputs into an array $bulk_inparams of the following object structure:
 # {
 #    object_ref => $asmb,
 #    output_workspace => $ws,
 #    output_genome_name => $out_genomeSet . '_' .$obj_name,
 #    scientific_name => $scientific_name,
 #    genetic_code => $genetic_code,
 #    domain => $domain,
 #    create_report => 0
 # }
##
#
sub _get_bulk_rast_parameters {
    my $self = shift;
    my ($params) = @_;

    my $ws = $params->{output_workspace};
    my $out_genomeSet = $params->{output_GenomeSet_name};
    my $in_assemblies = $params->{input_assemblies};
    my $in_genomes = $params->{input_genomes};
    my $in_text = $params->{input_text};

    my $chk = {};
    my $bulk_inparams = ();
    for my $asmb (@$in_assemblies) {
        $chk = $self->_validate_KB_objref_name($asmb);
        next unless $chk->{check_passed};
        next if $self->_value_in_array($asmb, $bulk_inparams);
        my $tmp_parm = $self->_build_param_from_obj($asmb, $chk, $ws, $out_genomeSet);
        push @$bulk_inparams, $tmp_parm if defined $tmp_parm;
    }

    for my $gn (@$in_genomes) {
        $chk = $self->_validate_KB_objref_name($gn);
        next unless $chk->{check_passed};
        next if $self->_value_in_array($gn, $bulk_inparams);
        my $tmp_parm = $self->_build_param_from_obj($gn, $chk, $ws, $out_genomeSet);
        push @{$bulk_inparams}, $tmp_parm if defined $tmp_parm;
    }

    if ($in_text) {
        my $input_list = [split(/[\n;\|]+/, $in_text)];
        $input_list = $self->_uniq_ref($input_list);
        for my $gn (@{$input_list}) {         
            $chk = $self->_validate_KB_objref_name($gn);
            next unless $chk->{check_passed};
            next if $self->_value_in_array($gn, $bulk_inparams);
            my $tmp_parm = $self->_build_param_from_obj($gn, $chk, $ws, $out_genomeSet);
            push @{$bulk_inparams}, $tmp_parm if defined $tmp_parm;
        }
    }
    return $bulk_inparams;
}

##
 # Loop through the array to call the above rast_genome function on each.
 # After creating a genomeSet in the workspace, return the following object:
 # {
 #     "output_genomeSet_ref"=>$ws."/".$out_genomeSet,
 #     "output_workspace"=>$ws,
 #     "report_name"=>$kbutil_output->{report_name},
 #     "report_ref"=>$kbutil_output->{report_ref}
 # }
##
sub bulk_rast_genomes {
    my ( $self, $inparams ) = @_;
    print "*********bulk_rast_genomes input parameter=\n". Dumper($inparams). "\n";

    my $params = $self->_check_bulk_annotation_params($inparams);

    my $ws = $params->{output_workspace};
    my $out_genomeSet = $params->{output_GenomeSet_name};
    my $bulk_inparams = $self->_get_bulk_rast_parameters($params);
    #
    # Throw an error IF $bulks_inparams is NOT a ref to a non-empty ARRAY
    #
    my $empty_input_msg = ("ERROR:Missing required inputs--must specify at least one genome \n".
            "and/or a string of genome names separated by ';', '\n' or '|' (without quotes).\n");
    Bio::KBase::Exceptions::ArgumentValidationError->throw(
        error        => $empty_input_msg,
        method_name  => 'bulk_rast_genomes'
    ) unless ref $bulk_inparams eq 'ARRAY' && @{$bulk_inparams};

    my $anngns = [];
    for my $parm (@{$bulk_inparams}) {
        my $rast_out = $self->rast_genome($parm);
        if (keys %{ $rast_out }) {
            push @$anngns, $rast_out->{output_genome_ref};
        }
    }

    # if $anngns is empty, return empty hash
    return {} unless @$anngns;
    print "Bulk rasted ".scalar @$anngns." genomes.\n"; # .Dumper($anngns);

    # create, save and then return that GenomeSet object's ref
    my $kbutil = installed_clients::kb_SetUtilitiesClient->new($self->{call_back_url});
    my $kbutil_output = $kbutil->KButil_Build_GenomeSet({
        workspace_name => $ws,
        input_refs => $anngns,
        output_name => $out_genomeSet,
        desc => 'GenmeSet generated from RAST annotated genomes/assemblies'
    });

    return {"output_genomeSet_ref"=>$ws."/".$out_genomeSet,
            "output_workspace"=>$ws,
            "report_name"=>$kbutil_output->{report_name},
            "report_ref"=>$kbutil_output->{report_ref}};
}


#
## for use with arrays, e.g.,
## my @unique = uniq(@data);
#
sub _uniq {
    my $self = shift;
    keys { map { $_ => 1 } @_ };
}

#
## for use with array refs, e.g.,
## my $unique = uniq($data);
#
sub _uniq_ref {
    my $self = shift;
    [ keys { map { $_ => 1 } @{$_[0]} } ];
}

#
## for use to check if a key value is in an array of hash
#
sub _value_in_array {
     my ($self, $val, $arr) = @_;
     for my $href (@{$arr}) {
        if ($href->{object_ref} eq $val) {
            return 1;
        }
     }
     return 0;
}

#
## Initialization
#
sub doInitialization {
    my $self = shift;

    $self->{_token} = $self->{ctx}->token();
    $self->{_username} = $self->{ctx}->user_id();
    $self->{_method} = $self->{ctx}->method();
    $self->{_provenance} = $self->{ctx}->provenance();

    $self->{ws_url} = $self->{config}->{'workspace-url'};
    $self->{call_back_url} = $ENV{ SDK_CALLBACK_URL };
    $self->{rast_scratch} = $self->{config}->{'scratch'};
    $self->{genome_dir} = $self->_create_rast_subdir($self->{rast_scratch},
                                                     "genome_annotation_dir_");

    $self->{max_contigs} = 10000;

    die "no workspace-url defined" unless $self->{ws_url};

    $self->{ws_client} = installed_clients::WorkspaceClient->new(
                             $self->{ws_url}, token => $self->{_token});

    return 1;
}

#
## constructor
#
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

