use strict;
use Data::Dumper;
use Test::More;
use Test::Exception;
use Config::Simple;
use Time::HiRes qw(time);
use Workspace::WorkspaceClient;
use JSON;
use File::Copy;
use AssemblyUtil::AssemblyUtilClient;
use Storable qw(dclone);
use RAST_SDK::RAST_SDKImpl

local $| = 1;
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('RAST_SDK');
my $ws_url = $config->{"workspace-url"};
my $ws_name = undef;
my $ws_client = new Workspace::WorkspaceClient($ws_url,token => $token);
my $auth_token = Bio::KBase::AuthToken->new(token => $token, ignore_authrc => 1);
my $ctx = LocalCallContext->new($token, $auth_token->user_id);
$RAST_SDK::RAST_SDKServer::CallContext = $ctx;
my $impl = new RAST_SDK::RAST_SDKImpl();

sub get_ws_name {
    if (!defined($ws_name)) {
        my $suffix = int(time * 1000);
        $ws_name = 'test_RAST_SDK_' . $suffix;
        $ws_client->create_workspace({workspace => $ws_name});
    }
    return $ws_name;
}

sub make_impl_call {
    my ($method, $params) = @_;
    # Prepare json file input
    my $input = {
        method => $method,
        params => $params,
        version => "1.1",
        id => "1"
    };
    open my $fh, ">", "/kb/module/work/input.json";
    print $fh encode_json($input);
    close $fh;
    my $output_path = "/kb/module/work/output.json";
    if (-e $output_path) {
        unlink($output_path);
    }
    # Run run_async.sh
    system("sh", "/kb/module/scripts/run_async.sh");
    # Load json file with output
    open my $fh2, "<", $output_path;
    my $output_json = <$fh2>;
    close $fh2;
    my $json = JSON->new;
    my $output = $json->decode($output_json);
    if (defined($output->{error})) {
        die encode_json($output->{error});
    }
    my $ret_obj = $output->{result}->[0];
    return $ret_obj;
}

sub prepare_assembly {
    my($assembly_obj_name) = @_;
    my $fasta_data_path = "/kb/module/test/data/bogus.fna";
    my $fasta_temp_path = "/kb/module/work/tmp/bogus.fna";
    copy $fasta_data_path, $fasta_temp_path;
    my $call_back_url = $ENV{ SDK_CALLBACK_URL };
    my $au = new AssemblyUtil::AssemblyUtilClient($call_back_url);
    my $ret = $au->save_assembly_from_fasta({
        file => {path => $fasta_temp_path},
        workspace_name => get_ws_name(),
        assembly_name => $assembly_obj_name
    });
    unlink($fasta_temp_path);
    return $ret;
}

sub check_genome_obj {
    my($genome_obj) = @_;
    ok(defined($genome_obj->{features}), "Features array is present");
    ok(scalar @{ $genome_obj->{features} } eq 1, "Number of features");
    ok(defined($genome_obj->{cdss}), "CDSs array is present");
    ok(scalar @{ $genome_obj->{cdss} } eq 1, "Number of CDSs");
    ok(defined($genome_obj->{mrnas}), "mRNAs array is present");
    ok(scalar @{ $genome_obj->{mrnas} } eq 1, "Number of mRNAs");
}

sub test_annotate_assembly {
    my($assembly_obj_name) = @_;
    my $genome_obj_name = "genome.1";
    my $params={"input_contigset"=>$assembly_obj_name,
             "scientific_name"=>'bogus',
             "domain"=>'B',
             "genetic_code"=>'11',
             "call_features_rRNA_SEED"=>'1',
             "call_features_tRNA_trnascan"=>'1',
             "call_selenoproteins"=>'1',
             "call_pyrrolysoproteins"=>'1',
             "call_features_repeat_region_SEED"=>'1',
             "call_features_insertion_sequences"=>'0',
             "call_features_strep_suis_repeat"=>'1',
             "call_features_strep_pneumo_repeat"=>'1',
             "call_features_crispr"=>'1',
             "call_features_CDS_glimmer3"=>'0',
             "call_features_CDS_prodigal"=>'1',
             "annotate_proteins_kmer_v2"=>'1',
             "kmer_v1_parameters"=>'1',
             "annotate_proteins_similarity"=>'1',
             "resolve_overlapping_features"=>'1',
             "find_close_neighbors"=>'1',
             "call_features_prophage_phispy"=>'0',
             "output_genome"=>$genome_obj_name,
             "workspace"=>get_ws_name()
           };
    $impl->annotate_genome($params);
    #my $ret = make_impl_call("RAST_SDK.annotate_genome", $params);
    my $genome_ref = get_ws_name() . "/" . $genome_obj_name;
    my $genome_obj = $ws_client->get_objects([{ref=>$genome_ref}])->[0]->{data};
    open my $fh, ">", "/kb/module/work/tmp/bogus_genome.json";
    print $fh encode_json($genome_obj);
    close $fh;
    # no detailed checking right now
    #check_genome_obj($genome_obj);
}

sub load_genome_from_json {
    my($assembly_ref) = @_;
    my $genome_obj_dir = "/kb/module/work/tmp/";
    my $genome_file_name = "bogus_genome.json";
    my $genome_obj_path = $genome_obj_dir . $genome_file_name;
    unless (-e $genome_obj_path) {
        $genome_obj_dir = "/kb/module/test/data/";
        $genome_obj_path = $genome_obj_dir . $genome_file_name;
    }
    open my $fh2, "<", $genome_obj_path;
    my $genome_json = <$fh2>;
    close $fh2;
    my $json = JSON->new;
    my $genome_obj = $json->decode($genome_json);
    $genome_obj->{assembly_ref} = $assembly_ref;
    # TODO: we need to set $genome_obj->{features}->[*]->{ontology_terms}->{SSO}->{"*"}->{ontology_ref} = "KBaseOntology/seed_subsystem_ontology";
    # And the same with $genome_obj->{cdss}
    return $genome_obj;
}

sub save_genome_to_ws {
    my($genome_obj,$genome_obj_name) = @_;
    my $ret = $ws_client->save_objects({workspace=>get_ws_name(),objects=>[{data=>$genome_obj,
            type=>"KBaseGenomes.Genome", name=>$genome_obj_name}]})->[0];
    return $ret->[6]."/".$ret->[0]."/".$ret->[4];
}

sub prepare_new_genome {
    my($assembly_obj_name,$genome_obj_name) = @_;
    my $genome_obj = load_genome_from_json($assembly_obj_name);
    return save_genome_to_ws($genome_obj, $genome_obj_name);
}

sub test_reannotate_genome {
    my($genome_obj_name, $genome_ref) = @_;
    my $params={"input_genome"=>$genome_obj_name,
             "call_features_rRNA_SEED"=>'0',
             "call_features_tRNA_trnascan"=>'0',
             "call_selenoproteins"=>'0',
             "call_pyrrolysoproteins"=>'0',
             "call_features_repeat_region_SEED"=>'0',
             "call_features_insertion_sequences"=>'0',
             "call_features_strep_suis_repeat"=>'0',
             "call_features_strep_pneumo_repeat"=>'0',
             "call_features_crispr"=>'0',
             "call_features_CDS_glimmer3"=>'0',
             "call_features_CDS_prodigal"=>'0',
             "annotate_proteins_kmer_v2"=>'0',
             "kmer_v1_parameters"=>'0',
             "annotate_proteins_similarity"=>'1',
             "resolve_overlapping_features"=>'0',
             "find_close_neighbors"=>'1',
             "call_features_prophage_phispy"=>'0',
             "output_genome"=>$genome_obj_name,
             "workspace"=>get_ws_name()
           };
    if (defined $genome_ref){
        $params->{input_genome} = $genome_ref;
    }
    return $impl->annotate_genome($params);
    #my $ret = make_impl_call("RAST_SDK.annotate_genome", $params);
    my $genome_ref = get_ws_name() . "/" . $genome_obj_name;
    my $genome_obj = $ws_client->get_objects([{ref=>$genome_ref}])->[0]->{data};
    # no detailed checking right now
    #check_genome_obj($genome_obj);
}

sub prepare_old_genome {
    my($assembly_obj_name,$genome_obj_name) = @_;
    my $genome_obj = load_genome_from_json($assembly_obj_name);
    delete $genome_obj->{cdss};
    delete $genome_obj->{mrnas};
    for (my $i=0; $i < @{$genome_obj->{features}}; $i++) {
        my $ftr = $genome_obj->{features}->[$i];
        $ftr->{type} = "CDS";
        delete $ftr->{cdss};
        delete $ftr->{mrnas};
    }
    return save_genome_to_ws($genome_obj, $genome_obj_name);
}

sub prepare_recent_old_genome {
    my($assembly_obj_name,$genome_obj_name) = @_;
    my $genome_obj = load_genome_from_json($assembly_obj_name);
    delete $genome_obj->{cdss};
    delete $genome_obj->{mrnas};
    my $genes = [];
    for (my $i=0; $i < @{$genome_obj->{features}}; $i++) {
        my $ftr = $genome_obj->{features}->[$i];
        $ftr->{type} = "CDS";
        delete $ftr->{cdss};
        delete $ftr->{mrnas};
	my $ftr2 = dclone $ftr;
        $ftr2->{type} = "gene";
        delete $ftr2->{protein_translation};
        $ftr2->{location}->[0]->[3] += 1;
        push(@{$genes}, $ftr2);
    }
    for (my $i=0; $i < @{$genes}; $i++) {
        push(@{$genome_obj->{features}}, $genes->[$i]);
    }
    return save_genome_to_ws($genome_obj, $genome_obj_name);
}

my $assembly_obj_name = "contigset.1";
my $assembly_ref = prepare_assembly($assembly_obj_name);
lives_ok {
        test_annotate_assembly($assembly_obj_name);
    }, "test_annotate_assembly";
    my $genome_obj_name = "genome.2";
lives_ok{
        my $genome_ref = prepare_new_genome($assembly_ref, $genome_obj_name);
        test_reannotate_genome($genome_obj_name, $genome_ref);
    }, "test_reannotate_genome";
lives_ok{
        $genome_obj_name = "genome.3";
        my $genome_ref = prepare_old_genome($assembly_ref, $genome_obj_name);
        test_reannotate_genome($genome_obj_name, $genome_ref);
    }, "test_reannotate_genome";
lives_ok{
        $genome_obj_name = "genome.4";
        my $genome_ref = prepare_recent_old_genome($assembly_ref, $genome_obj_name);
        test_reannotate_genome($genome_obj_name, $genome_ref);
    }, 'test_reannotate_genome';
done_testing(4);

my $err = undef;
if ($@) {
    $err = $@;
}
eval {
    if (defined($ws_name)) {
        $ws_client->delete_workspace({workspace => $ws_name});
        print("Test workspace was deleted\n");
    }
};
if (defined($err)) {
    if(ref($err) eq "Bio::KBase::Exceptions::KBaseException") {
        die("Error while running tests: " . $err->trace->as_string);
    } else {
        die $err;
    }
}

