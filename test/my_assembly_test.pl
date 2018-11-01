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
use File::Slurp;

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
    open my $fh2, "<", $output_path;
    my $output_json = <$fh2>;
    #print $output_json;
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
    my $fasta_data_path = "/kb/module/test/data/Clostridium_thermocellum_ATCC27405.fa";
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
    ok(scalar @{ $genome_obj->{features} } gt 0, "Number of features");
    ok(defined($genome_obj->{cdss}), "CDSs array is present");
    ok(scalar @{ $genome_obj->{cdss} } gt 0, "Number of CDSs");
    ok(defined($genome_obj->{mrnas}), "mRNAs array is present");
    ok(scalar @{ $genome_obj->{mrnas} } gt 0, "Number of mRNAs");
}

sub test_annotate_assembly {
    my($assembly_obj_name) = @_;
    my $genome_obj_name = "genome.1";
    my $params={"input_contigset"=>$assembly_obj_name,
             "scientific_name"=>'bogus',
             "domain"=>'B',
             "genetic_code"=>'11',
             "call_features_rRNA_SEED"=>'0',
             "call_features_tRNA_trnascan"=>'1',
             "call_selenoproteins"=>'0',
             "call_pyrrolysoproteins"=>'0',
             "call_features_repeat_region_SEED"=>'0',
             "call_features_insertion_sequences"=>'0',
             "call_features_strep_suis_repeat"=>'0',
             "call_features_strep_pneumo_repeat"=>'0',
             "call_features_crispr"=>'0',
             "call_features_CDS_glimmer3"=>'0',
             "call_features_CDS_prodigal"=>'1',
             "annotate_proteins_kmer_v2"=>'1',
             "kmer_v1_parameters"=>'0',
             "annotate_proteins_similarity"=>'0',
             "resolve_overlapping_features"=>'0',
             "call_features_prophage_phispy"=>'0',
             "output_genome"=>$genome_obj_name,
             "workspace"=>get_ws_name()
           };
    #$impl->annotate_genome($params);
    my $ret = make_impl_call("RAST_SDK.annotate_genome", $params);
    my $genome_ref = get_ws_name() . "/" . $genome_obj_name;
    my $genome_obj = $ws_client->get_objects([{ref=>$genome_ref}])->[0]->{data};
    open my $fh, ">", "/kb/module/work/tmp/bogus_genome.json";
    print $fh encode_json($genome_obj);
    close $fh;
    # no detailed checking right now
    check_genome_obj($genome_obj);
}


sub save_genome_to_ws {
    my($genome_obj,$genome_obj_name) = @_;
    #my $ret = $ws_client->save_objects({workspace=>get_ws_name(),objects=>[{data=>$genome_obj,
    #        type=>"KBaseGenomes.Genome-11.0", name=>$genome_obj_name}]})->[0];
    my $ret = $gaa->save_one_genome_v1({workspace=>get_ws_name(), data=>$genome_obj, name=>$genome_obj_name})->{info};
    return $ret->[6]."/".$ret->[0]."/".$ret->[4];
}




my $assembly_obj_name = "contigset.1";
my $assembly_ref = prepare_assembly($assembly_obj_name);
my $genome_ref_new;
my $genome_ref_old;
lives_ok {
        test_annotate_assembly($assembly_obj_name);
    }, "test_annotate_assembly";
done_testing(7);

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

