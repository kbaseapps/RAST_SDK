use strict;
use Data::Dumper;
use Test::More;
use Config::Simple;
use Time::HiRes qw(time);
use Bio::KBase::AuthToken;
use Bio::KBase::workspace::Client;
use JSON;
use File::Copy;
use AssemblyUtil::AssemblyUtilClient;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday };
};

local $| = 1;
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('RAST_SDK');
my $ws_url = $config->{"workspace-url"};
my $ws_name = undef;
my $ws_client = new Bio::KBase::workspace::Client($ws_url,token => $token);
my $auth_token = Bio::KBase::AuthToken->new(token => $token, ignore_authrc => 1);

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
    #chdir('/kb/module');
    system("sh", "/kb/module/scripts/run_async.sh");
    # Load json file with output
    open my $fh2, "<", $output_path;
    my $output_json = <$fh2>;
    close $fh2;
    my $json = JSON->new;
    my $output = $json->decode($output_json);
    if (defined($output->{error})) {
        die $output->{error};
    }
    my $ret_obj = $output->{result}->[0];
    return $ret_obj;
}

eval {
    my $obj_name = "contigset.1";
    my $fasta_data_path = "/kb/module/test/data/bogus.fna";
    my $fasta_temp_path = "/kb/module/work/tmp/bogus.fna";
    copy $fasta_data_path, $fasta_temp_path;
    my $call_back_url = $ENV{ SDK_CALLBACK_URL };
    my $au = new AssemblyUtil::AssemblyUtilClient($call_back_url);
    $au->save_assembly_from_fasta({
        file => {path => $fasta_temp_path},
        workspace_name => get_ws_name(),
        assembly_name => $obj_name
    });
    my $genome_obj_name = "genome.1";
    my $params={"input_contigset"=>$obj_name,
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
             "find_close_neighbors"=>'0',
             "call_features_prophage_phispy"=>'0',
             "output_genome"=>$genome_obj_name,
             "workspace"=>get_ws_name()
           };
    my $ret = make_impl_call("RAST_SDK.annotate_genome", $params);
    print encode_json($ret);
    my $genome_ref = get_ws_name() . "/" . $genome_obj_name;
    my $genome_obj = $ws_client->get_objects([{ref=>$genome_ref}])->[0]->{data};
    open my $fh, ">", "/kb/module/work/tmp/bogus_genome.json";
    print $fh encode_json($genome_obj);
    close $fh;
    done_testing(1);
};
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

