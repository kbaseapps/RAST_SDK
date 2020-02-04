use strict;
use Data::Dumper;
use Test::More;
use Test::Exception;
use Config::Simple;
use Time::HiRes qw(time);
use Workspace::WorkspaceClient;
use JSON;
use File::Copy;
use GenomeAnnotationAPI::GenomeAnnotationAPIClient;
use Storable qw(dclone);
use File::Slurp;

use lib "/kb/module/test";
use testRASTutil;

print "PURPOSE:\n";
print "    1.  Test annotation of genetic code 25 for Mircea Podar. \n";
print "    2.  Test the Prodigal is using the -meta when the domain is 'Unknown'\n\n";

local $| = 1;
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('RAST_SDK');
my $ws_url = $config->{"workspace-url"};
my $ws_name = undef;
my $ws_client = new Workspace::WorkspaceClient($ws_url,token => $token);
my $call_back_url = $ENV{ SDK_CALLBACK_URL };
my $gaa = new GenomeAnnotationAPI::GenomeAnnotationAPIClient($call_back_url);

my $assembly_obj_name = "GCA_000350285.1_OR1_genomic.fna";
my $assembly_ref = prepare_assembly($assembly_obj_name);
my $genome_obj_name = 'SR1_bacterium_MGEHA';
my $genome_set_name = "New_GenomeSet";

my $params={"input_contigset"=>$assembly_obj_name,
             "scientific_name"=>'candidate division SR1 bacterium MGEHA',
             "domain"=>'U',
             "genetic_code"=>'25',
             "call_features_CDS_prodigal"=>'1',
           };

$params = &set_params($genome_obj_name,$params);

lives_ok {
	print("######## Running RAST annotation ########\n");
	my $ret = &make_impl_call("RAST_SDK.annotate_genome", $params);
	my $genome_ref = get_ws_name() . "/" . $genome_obj_name;
	my $genome_obj = $ws_client->get_objects([{ref=>$genome_ref}])->[0]->{data};
    
	print "\n\nOUTPUT OBJECT DOMAIN = $genome_obj->{domain}\n";
	print "OUTPUT OBJECT G_CODE = $genome_obj->{genetic_code}\n";
	print "OUTPUT SCIENTIFIC NAME = $genome_obj->{scientific_name}\n";

    ok(defined($genome_obj->{features}), "Features array is present");
    ok(scalar @{ $genome_obj->{features} } gt 0, "Number of features");
}, "test_annotate_assembly";
print "Summary for $assembly_obj_name\n";

done_testing(10);

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

