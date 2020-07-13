use strict;
use Data::Dumper;
use Test::More;
use Test::Exception;
use Config::Simple;
use Time::HiRes qw(time);
use JSON;
use File::Copy;
use installed_clients::WorkspaceClient;
use installed_clients::AssemblyUtilClient;
use installed_clients::GenomeFileUtilClient;
use installed_clients::kb_SetUtilitiesClient;
use Storable qw(dclone);
use Bio::KBase::KBaseEnv;

use RASTTestUtils;

print "Running $0\n";
print "PURPOSE:\n";
print "    1.  Test of Multiple Genomes and GenomeSets. \n";
print "        a. Two individual genomes and annotate.\n";
print "        b. Take the two individual genomes, make a set and annotate.\n";
print "        c. Take the genomeSet, add an individual genomes and annotate.\n";
print "    2.  The optional parameters are minimal because the sets are the test.\n";
print "    3.  In debug mode, using a genome reference in CI/prod.\n";
print "         Otherwise, load a genome in data dir. This takes more time but every user has access.\n";
print "         For this reason, the test isn't for a specific number of changes\n\n";

local $| = 1;
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = Config::Simple->new($config_file)->get_block('RAST_SDK');
my $ws_url = $config->{"workspace-url"};
my $ws_name = undef;
my $ws_client = installed_clients::WorkspaceClient->new($ws_url,token => $token);
my $call_back_url = $ENV{ SDK_CALLBACK_URL };
my $au = installed_clients::AssemblyUtilClient->new($call_back_url);
my $gfu = installed_clients::GenomeFileUtilClient->new($call_back_url);
my $su = installed_clients::kb_SetUtilitiesClient->new($call_back_url);

my $DEBUG = 'N';
my $genome_obj_name1 = "Dactylopius_coccus";
my $genome_ref1 = "40046/2/1";
my $genome_obj_name2 = "Drosophila_melanogaster";
my $genome_ref2 = "40046/3/1";
my $genome_obj_name3 = "Nomada_ferruginata";
my $genome_ref3 = "40046/4/1";

# PROD
#my $genome_ref1 = "19217/183555/4";
#my $genome_ref2 = "19217/191347/4";
#my $genome_ref3 = "19217/165153/4";

if ($DEBUG ne 'Y') {
	my $tmp_obj;
	$genome_obj_name1 = "Carsonella";
	my $genome_gbff_name1 = "Carsonella.gbk";
	($tmp_obj, $genome_ref1) = prepare_gbff($genome_gbff_name1,$genome_obj_name1);
	$genome_obj_name2 = "Methanosarcina";
	my $genome_gbff_name2 = "Methanosarcina_acetivorans_C2A.gbff";
	($tmp_obj, $genome_ref2) = prepare_gbff($genome_gbff_name2,$genome_obj_name2);
	my $genome_gbff_name3 = "GCF_000287295.1_ASM28729v1_genomic.gbff";
	($tmp_obj, $genome_ref3) = prepare_gbff($genome_gbff_name3,$genome_obj_name3);

}

my $params={"input_genomes"=>[$genome_ref1,$genome_ref2],
             "call_features_tRNA_trnascan"=>'1',
           };

my ($orig_genome1,$orig_funcs1) = &get_genome($genome_ref1);
my ($orig_genome2,$orig_funcs2) = &get_genome($genome_ref2);
my ($orig_genome3,$orig_funcs3) = &get_genome($genome_ref3);

#
#	TEST TWO GENOMES -
#
lives_ok {

	if ($DEBUG ne 'Y') {
		$params->{input_genomes} = [$genome_ref1,$genome_ref2];
		my ($genome_set_obj) = &submit_set_annotation('multi_genomes', $params->{input_genomes}, $params);

		my $data = $ws_client->get_objects([{ref=>$genome_set_obj}])->[0]->{refs};
		my $number_genomes = scalar @{ $data};
    	ok($number_genomes == 2, "Input: Two Genomes. Output: $number_genomes in output GenomeSet");
	}
	1;
} "Two individual Genome";

#
#	BUILD AND TEST A GENOME SET
#

lives_ok {
	if ($DEBUG ne 'Y') {
		my $genome_set_name = 'new_genome_set';
		my $genome_set = $su->KButil_Build_GenomeSet({
        	workspace_name => get_ws_name(),
        	input_refs => [$genome_ref1,$genome_ref2],
        	output_name => $genome_set_name,
        	desc => 'GenomeSet Description'
    	});

		my $gs = get_ws_name() . "/" . $genome_set_name ;
		my $info = $ws_client->get_objects([{ref=>$gs}])->[0]->{info};
		my $newref = $info->[6]."/".$info->[0]."/".$info->[4];
		$params->{input_genomes} = [$newref];

		my ($genome_set_obj) = &submit_set_annotation($genome_set_name, $params->{input_genomes}, $params);

		my $data = $ws_client->get_objects([{ref=>$genome_set_obj}])->[0]->{refs};
		my $number_genomes = scalar @{ $data};
    	ok($number_genomes == 2, "Input: One GenomeSet. Output: $number_genomes in output GenomeSet");
	} else {
		1;
	}
} "Build and Test a Set";

#
#	BUILD AND TEST A GENOME SET PLUS ANOTHER GENOME
#

lives_ok {
	if ($DEBUG ne 'Y') {
		my $genome_set_name = 'genome_plus_set';
		my $genome_set = $su->KButil_Build_GenomeSet({
        	workspace_name => get_ws_name(),
        	input_refs => [$genome_ref1,$genome_ref2],
        	output_name => $genome_set_name,
        	desc => 'GenomeSet Description'
    	});

		my $gs = get_ws_name() . "/" . $genome_set_name ;
		my $info = $ws_client->get_objects([{ref=>$gs}])->[0]->{info};
		my $newref = $info->[6]."/".$info->[0]."/".$info->[4];
		$params->{input_genomes} = [$newref,$genome_ref3];

		my ($genome_set_obj) = &submit_set_annotation($genome_set_name, $params->{input_genomes}, $params);

		my $data = $ws_client->get_objects([{ref=>$genome_set_obj}])->[0]->{refs};
		my $number_genomes = scalar @{ $data};
                ok($number_genomes == 2, "Input: GenomeSet plus one. Output: $number_genomes in output GenomeSet");
	} else {
		1;
	}
} "Test a Set Plus one";
done_testing(6);

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


