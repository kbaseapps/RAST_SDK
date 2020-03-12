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
use installed_clients::GenomeFileUtilClient;
use installed_clients::kb_SetUtilitiesClient;
use Storable qw(dclone);
use Bio::KBase::kbaseenv;

use lib "/kb/module/test";
use testRASTutil;

print "PURPOSE:\n";
print "    1.  Test annotate Multiple Assemblies. \n";
print "        Test two assemblies, an assemblySet, an asseblySet plus singleton, and redundant assemblies\n";
print "    2.  Minimum test for tRNA\n";
print "    3.  In debug mode, using a genome reference in CI/prod.\n";
print "         Otherwise, load an assembly in data dir. This takes more time but every user has access.\n";
print "         For this reason, the tests aren't for a specific number of changes or names.\n\n";

local $| = 1;
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('RAST_SDK');
my $ws_url = $config->{"workspace-url"};
my $ws_name = undef;
my $ws_client = new Workspace::WorkspaceClient($ws_url,token => $token);
my $call_back_url = $ENV{ SDK_CALLBACK_URL };
my $au = new AssemblyUtil::AssemblyUtilClient($call_back_url);
my $gfu = new installed_clients::GenomeFileUtilClient($call_back_url);
my $su = new installed_clients::kb_SetUtilitiesClient($call_back_url);


my $DEBUG = 'N';
my $assembly_obj_name1 = "Dactylopius coccus";
my $assembly_ref1 = "40046/5/1";
   $assembly_ref1 = "40619/11/1";
my $assembly_obj_name2 = "Drosophila melanogaster";
my $assembly_ref2 = "40046/6/1";
   $assembly_ref2 = "40619/39/1";
my $assembly_obj_name3 = "Nomada ferruginata";
my $assembly_ref3 = "40046/7/1";
my $genome_set_name = "New_GenomeSet";


if ($DEBUG ne 'Y') {
	$assembly_obj_name1 = "bogus.fna";
	$assembly_ref1 = prepare_assembly($assembly_obj_name1);
		
	$assembly_obj_name2 = "bogus2.fna";
	$assembly_ref2 = prepare_assembly($assembly_obj_name2);
	
	$assembly_obj_name3 = "bogus2.fna";
	$assembly_ref3 = prepare_assembly($assembly_obj_name2);
}
my $params={"input_genomes"=>'',
             "call_features_tRNA_trnascan"=>'1',
			"output_genome"=>$genome_set_name
           };

#
#	Set up the needed AssemblySet
#
my $assembly_set_name = 'new_assembly_set';
my $assembly_set = $su->KButil_Build_AssemblySet({
     	workspace_name => get_ws_name(),
      	input_refs => [$assembly_ref1,$assembly_ref2],
      	output_name => $assembly_set_name,
      	desc => 'Test AssemblySet Description'
  	});

#
#	TEST TWO ASSEMBLIES -
#
print "ASSEMBLYREF1 = $assembly_ref1 and ASSEMBLYREF2 = $assembly_ref2\n";
lives_ok {
	if ($DEBUG ne 'Y') {
		my $genome_refs  = [$assembly_ref1,$assembly_ref2];
		$params->{input_genomes} = $genome_refs;
		my ($genome_set_obj) = &submit_set_annotation($genome_set_name, $genome_refs, $params);

		my $data = $ws_client->get_objects([{ref=>$genome_set_obj}])->[0]->{refs};
		my $number_genomes = scalar @{ $data};
    	ok($number_genomes == 2, "Input: Two Assemblies. Output: $number_genomes in output GenomeSet");
	} else {
		1;
	}
} "Two Assemblies";

#
#	BUILD AND TEST A ASSEMBLY SET
#

lives_ok {
	if ($DEBUG ne 'Y') {
		my $gs = get_ws_name() . "/" . $assembly_set_name ;
		my $info = $ws_client->get_objects([{ref=>$gs}])->[0]->{info};
		my $newref = $info->[6]."/".$info->[0]."/".$info->[4];
		$params->{input_genomes} = [$newref];

		my ($genome_set_obj) = &submit_set_annotation($genome_set_name, $params->{input_genomes}, $params);
		my $data = $ws_client->get_objects([{ref=>$genome_set_obj}])->[0]->{refs};
		my $number_genomes = scalar @{ $data};
    	ok($number_genomes == 2, "Input: AssemblySet with two. Output: $number_genomes in output GenomeSet");
	} else {
		1;
	}
} "AssemblySet with two Assemblies";

#
#	BUILD AND TEST AN ASSEMBLY SET PLUS ANOTHER ASSEMBLY
#

lives_ok {
	if ($DEBUG ne 'Y') {
		my $gs = get_ws_name() . "/" . $assembly_set_name ;
		my $info = $ws_client->get_objects([{ref=>$gs}])->[0]->{info};
		my $newref = $info->[6]."/".$info->[0]."/".$info->[4];
		$params->{input_genomes} = [$newref,$assembly_ref3];

		my ($genome_set_obj) = &submit_set_annotation($genome_set_name, $params->{input_genomes}, $params);
		my $data = $ws_client->get_objects([{ref=>$genome_set_obj}])->[0]->{refs};
		my $number_genomes = scalar @{ $data};
    	ok($number_genomes == 3, "Input: AssemblySet plus one. Output: $number_genomes in output GenomeSet");

	} else {
		1;
	}
} "AssemblySet plus single assembly";

#
#	TEST REDUNDANT ASSEMBLIES -
#
lives_ok {
	if ($DEBUG ne 'Y') {
		$params->{input_genomes} = [$assembly_ref1,$assembly_ref1];
		my ($genome_set_obj) = &submit_set_annotation($genome_set_name, $params->{input_genomes}, $params);
		my $data = $ws_client->get_objects([{ref=>$genome_set_obj}])->[0]->{refs};
		my $number_genomes = scalar @{ $data};
    	ok($number_genomes = 1, "Input: Two redundant. Output: $number_genomes in output GenomeSet");
		
	} else {
		1;
	}
} "Redundant Assemblies";

done_testing(8);

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


