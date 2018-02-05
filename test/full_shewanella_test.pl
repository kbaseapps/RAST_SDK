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
use GenomeAnnotationAPI::GenomeAnnotationAPIClient;
use Storable qw(dclone);

local $| = 1;
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('RAST_SDK');
my $ws_url = $config->{"workspace-url"};
my $ws_name = undef;
my $ws_client = new Workspace::WorkspaceClient($ws_url,token => $token);
my $call_back_url = $ENV{ SDK_CALLBACK_URL };
my $gaa = new GenomeAnnotationAPI::GenomeAnnotationAPIClient($call_back_url);

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
    unless (-e $output_path) {
        die "Output file of JSON-RPC call (in CLI mode) wasn't created";
    }
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
    my $fasta_data_path = "/kb/module/test/data/shewy.fna";
    my $fasta_temp_path = "/kb/module/work/tmp/shewy.fna";
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

sub load_genome_from_json {
    my($assembly_obj_name,$genome_obj_path) = @_;
    my $genome_json;
    {
        local $/;  #Enable 'slurp' mode
        open my $fh2, "<", $genome_obj_path;
        $genome_json = <$fh2>;
        close $fh2;
    }
    my $json = JSON->new;
    my $genome_obj = $json->decode($genome_json);
    $genome_obj->{assembly_ref} = get_ws_name() . "/" . $assembly_obj_name;
    if (defined($genome_obj->{taxon_ref})) {
        delete $genome_obj->{taxon_ref};
    }
    if (defined($genome_obj->{genbank_handle_ref})) {
        delete $genome_obj->{genbank_handle_ref};
    }
    # TODO: we need to set $genome_obj->{features}->[*]->{ontology_terms}->{SSO}->{"*"}->{ontology_ref} = "KBaseOntology/seed_subsystem_ontology";
    # And the same with $genome_obj->{cdss}
    return $genome_obj;
}

sub save_genome_to_ws {
    my($genome_obj,$genome_obj_name) = @_;
    #my $info = $ws_client->save_objects({workspace=>get_ws_name(),objects=>[{data=>$genome_obj,
    #        type=>"KBaseGenomes.Genome-12.1", name=>$genome_obj_name}]})->[0];
    my $info = $gaa->save_one_genome_v1({workspace=>get_ws_name(), data=>$genome_obj, name=>$genome_obj_name})->{info};
    return $info->[6]."/".$info->[0]."/".$info->[4];
}

sub prepare_genome {
    my($assembly_obj_name,$genome_obj_name,$input_json_file) = @_;
    my $genome_obj = load_genome_from_json($assembly_obj_name,$input_json_file);
    my $genome_upa = save_genome_to_ws($genome_obj, $genome_obj_name);
    return ($genome_obj, $genome_upa);
}

sub reannotate_genome {
    my($genome_obj_name, $genome_upa) = @_;
    my $params={"input_genome"=>$genome_upa,
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
             "annotate_proteins_kmer_v2"=>'1',
             "kmer_v1_parameters"=>'1',
             "annotate_proteins_similarity"=>'0',
             "resolve_overlapping_features"=>'0',
             "find_close_neighbors"=>'0',
             "call_features_prophage_phispy"=>'0',
             "output_genome"=>$genome_obj_name,
             "workspace"=>get_ws_name()
           };
    return make_impl_call("RAST_SDK.annotate_genome", $params);
}

my $diff_count = 0;
lives_ok {
    print("Loading assembly to WS\n");
    my $assembly_obj_name = "assembly.1";
    prepare_assembly($assembly_obj_name);
    my $genome_obj_name = "genome.1";
    print("Loading genome to WS\n");
    my $genome_json_path = "/kb/module/test/data/shewy_genome.json";
    my ($orig_genome, $genome_upa) = prepare_genome($assembly_obj_name, $genome_obj_name, $genome_json_path);
    my $orig_funcs = {};
    for (my $i=0; $i < @{$orig_genome->{features}}; $i++) {
        my $ftr = $orig_genome->{features}->[$i];
        my $func = $ftr->{function};
        if (not defined($func)) {
            $func = "";
        }
        $orig_funcs->{$ftr->{id}} = $func;
    }
    open my $fh, ">", "/kb/module/work/tmp/shewy_genome_funcs.json";
    print $fh encode_json($orig_funcs);
    close $fh;
    print("Running RAST annotation\n");
    my $ret = reannotate_genome($genome_obj_name, $genome_upa);
    my $report_ref = $ret->{report_ref};
    my $report_obj = $ws_client->get_objects([{ref=>$report_ref}])->[0]->{data};
    my $report_json = encode_json($report_obj);
    my $json = JSON->new;
    my $report_text = $json->decode($report_json)->{text_message};
    print "Report: " . $report_text . "\n\n";
    my $genome_ref = get_ws_name() . "/" . $genome_obj_name;
    my $genome_obj = $ws_client->get_objects([{ref=>$genome_ref}])->[0]->{data};
    open my $fh, ">", "/kb/module/work/tmp/shewy_genome_reannot.json";
    print $fh encode_json($genome_obj);
    close $fh;
    for (my $i=0; $i < @{$genome_obj->{features}}; $i++) {
        my $ftr = $genome_obj->{features}->[$i];
        my $func = $ftr->{function};
        if (not defined($func)) {
            $func = "";
        }
        my $orig_func = $orig_funcs->{$ftr->{id}};
        if (not($func eq $orig_func)) {
            $diff_count++;
        }
    }
    } "Pipeline Runs";
print "Number of features with changed function: " . $diff_count . "\n";
ok($diff_count > 1000, "Number of updated features (" . $diff_count . ")");
done_testing(2);

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

