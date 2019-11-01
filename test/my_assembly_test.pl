use strict;
use Data::Dumper;
use Test::Deep;
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

use lib "/kb/module/test";
use testRASTutil;

print "PURPOSE:\n";
print "    1.  Test annotation of one small assembly. \n";
print "    2.  Test that the saved genome isn't using defaults. Must be Archaea and genetic code 4\n";
print "    3.  Test passing a taxon ID vs. a scientific name.\n";
print "    4.  Test 2 failure modes for contacting the RE.\n";
print "\n";

local $| = 1;
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('RAST_SDK');
my $ws_url = $config->{"workspace-url"};
my $ws_name = undef;
my $ws_client = new Workspace::WorkspaceClient($ws_url,token => $token);
my $call_back_url = $ENV{ SDK_CALLBACK_URL };
my $gaa = new GenomeAnnotationAPI::GenomeAnnotationAPIClient($call_back_url);

my $assembly_obj_name = "Acidilobus_sp._CIS.fna";
my $assembly_ref = prepare_assembly($assembly_obj_name);
my $genome_obj_name = 'Acidilobus_sp_CIS';
my $genome_set_name = "New_GenomeSet";

my $params={"input_contigset"=>$assembly_obj_name,
             "scientific_name"=>'Acidilobus sp 7',
             "domain"=>'A',
             "genetic_code"=>'4',
             "call_features_CDS_prodigal"=>'1',
           };

$params = &set_params($genome_obj_name,$params);

# Test processing a single assembly and saving as a genome.
lives_ok {
	print("######## Running RAST annotation ########\n");
	my $ret = &make_impl_call("RAST_SDK.annotate_genome", $params);
	my $genome_ref = get_ws_name() . "/" . $genome_obj_name;
	my $genome_obj = $ws_client->get_objects([{ref=>$genome_ref}])->[0]->{data};
    
	print "\n\nOUTPUT OBJECT DOMAIN = $genome_obj->{domain}\n";
	print "OUTPUT OBJECT G_CODE = $genome_obj->{genetic_code}\n";

    ok(defined($genome_obj->{features}), "Features array is present");
    ok(scalar @{ $genome_obj->{features} } gt 0, "Number of features");
    ok(defined($genome_obj->{cdss}), "CDSs array is present");
    ok(scalar @{ $genome_obj->{cdss} } gt 0, "Number of CDSs");
    ok(defined($genome_obj->{mrnas}), "mRNAs array is present");
    ok(scalar @{ $genome_obj->{mrnas} } gt 0, "Number of mRNAs");
    ok($genome_obj->{scientific_name} eq "Acidilobus sp 7", "Sci name is correct");
    ok(!defined($genome_obj->{taxon_assignments}), "Taxon assignments is undefined");
}, "test_annotate_assembly";
print "Summary for $assembly_obj_name\n";

# Test processing a single assembly and saving as a genome while specifing a NCBI taxon ID.
lives_ok {
	print("######## Running RAST annotation with Taxon ID ########\n");
    my $params_copy = { %$params };
    # this tax ID's species name changed in the 2018-12 NCBI dump and again in the 2019-02 dump
    # so it is a good test case for making sure the timestamp is passed to the RE correctly.
    # It depends on the RE containing 2018 NCBI data, which it currently does
    $params_copy->{ncbi_taxon_id} = 2448083;  # spec says this must be a string, we're sloppy tho
    $params_copy->{relation_engine_timestamp_ms} = 1545000000000;  # epoch ms

	my $ret = &make_impl_call("RAST_SDK.annotate_genome", $params_copy);
	my $genome_ref = get_ws_name() . "/" . $genome_obj_name;
	my $genome_obj = $ws_client->get_objects([{ref=>$genome_ref}])->[0]->{data};
    
	print "\n\nOUTPUT OBJECT DOMAIN = $genome_obj->{domain}\n";
	print "OUTPUT OBJECT G_CODE = $genome_obj->{genetic_code}\n";

    ok(defined($genome_obj->{features}), "Features array is present");
    ok(scalar @{ $genome_obj->{features} } gt 0, "Number of features");
    ok(defined($genome_obj->{cdss}), "CDSs array is present");
    ok(scalar @{ $genome_obj->{cdss} } gt 0, "Number of CDSs");
    ok(defined($genome_obj->{mrnas}), "mRNAs array is present");
    ok(scalar @{ $genome_obj->{mrnas} } gt 0, "Number of mRNAs");
    ok($genome_obj->{scientific_name} eq "Metarhizium sp. MJH 2018c", "Sci name is correct");
    cmp_deeply($genome_obj->{taxon_assignments}, {'ncbi' => '2448083'},
        "Taxon assignments is correct");
}, "test_annotate_assembly";
print "Summary for $assembly_obj_name\n";

# Test processing a single assembly and saving as a genome set.
print "ASSEMBLYREF = $assembly_ref\n";
lives_ok {

	my $params={"input_genomes"=>[$assembly_ref],
             "call_features_tRNA_trnascan"=>'1',
			"output_genome"=>$genome_set_name
           };

	my ($genome_set_obj,$params) = &submit_set_annotation($genome_set_name, $params->{input_genomes}, $params);
	my $data = $ws_client->get_objects([{ref=>$genome_set_obj}])->[0]->{refs};
	my $number_genomes = scalar @{ $data};
    ok($number_genomes == 1, "Input: One Assembly. Output: $number_genomes in output GenomeSet");
    
    my $genome_obj = $ws_client->get_objects([{ref=>$data->[0]}])->[0]->{data};
    ok($genome_obj->{scientific_name} eq "unknown taxon", "Sci name is correct");
    ok(!defined($genome_obj->{taxon_assignments}), "Taxon assignments is undefined");

	my $report = "/kb/module/work/tmp/annotation_report.$genome_set_name";
	my $local_path = "/kb/module/test/report_output/annotation_report.$genome_set_name";
    copy $report, $local_path;

	ok(-e $local_path,'File found');
} "Create a Report";
done_testing(23);

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

