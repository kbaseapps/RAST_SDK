use strict;
use Data::Dumper;
use Test::More;
use Test::Exception;
use Config::Simple;
use Time::HiRes qw(time);
use Workspace::WorkspaceClient;
use JSON;
use File::Copy;
use installed_clients::AssemblyUtilClient;
use installed_clients::GenomeFileUtilClient;
use Storable qw(dclone);
use Bio::KBase::KBaseEnv;

use RASTTestUtils;

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

print "    In debug mode, use a genome reference in CI/prod.\n";
print "    Otherwise, load a genome in data dir. This takes more time but every user has access.\n";
print "    For this reason, the tests aren't for a specific number of changes\n\n";

my $DEBUG = 'N';

my %list = ('carson' => 1, 'clos' => 1, 'meth' => 1, 'shew' => 1, 'strepp' => 1, 'streps' => 1);

my $num_tests = 0;

if (exists $list{carson} && $list{carson} > 0) {
	print "PURPOSE:\n";
	print "    1.  Test annotation of one small genome. (Carsonella)\n";
	print "    2.  Test for rRNA, tRNA, and protein_similarity\n";

	my $genome_ref = "15792/210698/1";
  	my $genome_gbff_name = "Carsonella.gbk";
	my $genome_obj_name = "Carsonella";

	my $params={
			"input_genome"=>$genome_ref,
             "call_features_rRNA_SEED"=>'1',
             "call_features_tRNA_trnascan"=>'1',
             "annotate_proteins_similarity"=>'1',
         };
	my %counts = &annotate($genome_obj_name,$genome_gbff_name, $genome_ref, $params);
	ok($counts{num_diff} > 0, "Number of new annotation (" . $counts{num_diff} . ")");
	ok($counts{num_new} > 0, "Number of new features (" . $counts{num_new} . ")");
	$num_tests += 3;
}

if (exists $list{clos} && $list{clos} > 0) {
	print "PURPOSE:\n";
	print "    1.  Clostridium botulinum is used as a test for selenocysteine-containing and for CRISPRS.\n";
	print "    2.  Both of these and one prophage should exist in this genome\n";

	my $genome_obj_name = "Clostridium_botulinum";
    my $genome_gbff_name = "Clostridium_botulinum_310.gbff";
	my $genome_ref = '15792/147849/1';

	my $params={
             "input_genome"=>$genome_ref,
             "call_selenoproteins"=>'1',
             "call_features_crispr"=>'1',
             "retain_old_anno_for_hypotheticals"=>'1',
             "call_features_prophage_phispy"=>'1',
           };
	my %counts = &annotate($genome_obj_name,$genome_gbff_name, $genome_ref, $params);
	ok($counts{seleno} == 0, "Number of Selenocysteine-containing genes = " . $counts{seleno} . "\n");
	ok($counts{crispr} == 0, "Number of features with CRISPR = " . $counts{crispr} . "\n");
	ok($counts{prophage} == 0, "Number of features with prophage = " . $counts{prophage} . "\n");
	$num_tests += 4;
}

if (exists $list{meth} && $list{meth} > 0) {
	print "PURPOSE:\n";
	print "    1.  Methanosarcina_acetivorans_C2A.gbff is used as a test for pyrrolysine-containing genes.\n";
	print "    2.  This should exist in this genome\n";

	my $genome_obj_name = "Methanosarcina_acetivorans_C2A";
    my $genome_gbff_name = "Methanosarcina_acetivorans_C2A.gbff";
	my $genome_ref = "15792/114497/2";

	my $params={
            "input_genome"=>$genome_ref,
             "call_pyrrolysoproteins"=>'1',
             "call_features_prophage_phispy"=>'1',
           };
	my %counts = &annotate($genome_obj_name,$genome_gbff_name, $genome_ref, $params);
	ok($counts{prophage} == 0, "Number of features with prophage = " . $counts{prophage} . "\n");
	ok($counts{pyrro} == 0, "Number of Pyrrolysine-containing genes = " . $counts{pyrro} . "\n");
	$num_tests += 3;

}

if (exists $list{shew} && $list{shew} > 0) {
	print "PURPOSE:\n";
	print "    1.  Shewanella_sp._ANA-3 is used as a test for RAST repeats, kmer v2 and kmer v2.\n";
	print "    2.  Both of these should exist in this genome\n";
	print "    3.  In debug mode, using a genome reference in CI/prod.\n";

	my $genome_obj_name = "Shewanella_sp._ANA-3";
	my $genome_gbff_name = "Shewanella_sp._ANA-3.gbff";
	my $genome_ref = "15792/209997/1";

	my $params={
            "input_genome"=>$genome_ref,
             "call_features_repeat_region_SEED"=>'1',
             "annotate_proteins_kmer_v2"=>'1',
             "kmer_v1_parameters"=>'1',
           };
	my %counts = &annotate($genome_obj_name,$genome_gbff_name, $genome_ref, $params);
	ok($counts{num_diff} > 0, "Number of new annotation (" . $counts{num_diff} . ")");
	ok($counts{num_new} > 0, "Number of new features (" . $counts{num_new} . ")");
	$num_tests += 3;
}

if (exists $list{strepp} && $list{strepp} > 0) {

	print "PURPOSE:\n";
	print "    1.  Streptococcus pneumoniae is a test for the strep-pneumoniae specific repeat\n";
	print "    2.  Streptococcus pneumoniae should also have SEED repeats\n";

	my $genome_obj_name = "Streptococcus_pneumoniae_D39";
	my $genome_gbff_name = "Streptococcus_pneumoniae_D39.gbff";
	my $genome_ref = "15792/103507/1";

	my $params={
            "input_genome"=>$genome_ref,
             "call_features_repeat_region_SEED"=>'1',
             "call_features_strep_pneumo_repeat"=>'1',
             "retain_old_anno_for_hypotheticals"=>'1',
           };
	my %counts = &annotate($genome_obj_name,$genome_gbff_name, $genome_ref, $params);
	ok($counts{repeat} > 0, "Number of repeats (" . $counts{repeat} . ")");
	$num_tests += 2;
}

if (exists $list{streps} && $list{streps} > 0) {
	print "PURPOSE:\n";
	print "    1.  Streptococcus suis is a test for the strep-suis specific repeat\n";
	print "    2.  This should exist in this genome\n";
	print "    3.  In debug mode, using a genome reference in CI/prod.\n";
	print "         Otherwise, load a genome in data dir. This takes more time but every user has access.\n";
	print "         For this reason, the test isn't for a specific number of changes\n\n";

	my $genome_obj_name = "Streptococcus_suis";
	my $genome_gbff_name = "Streptococcus_suis.gbff";
	my $genome_ref = '15792/60950/1';

	my $params={
	     "input_genome"=>$genome_ref,
	     "call_features_strep_suis_repeat"=>'1',
	     };
	my %counts = &annotate($genome_obj_name,$genome_gbff_name, $genome_ref, $params);
	ok($counts{repeat} > 0, "Number of repeats (" . $counts{repeat} . ")");
	$num_tests += 2;
}

sub annotate {
	my ($genome_obj_name,$genome_gbff_name, $genome_ref, $params) = @_;

	unless ($DEBUG eq 'Y')
	{
		my $tmp_genome_obj;
		($tmp_genome_obj, $genome_ref) = prepare_gbff($genome_gbff_name,$genome_obj_name);
        $params->{input_genome} = $genome_ref;
	}

	my $num_func_in  = 0;
	my $num_func_out = 0;
	my %counts = ('num_new' => 0, 'num_diff' => 0,
		'seleno' => 0, 'crispr' => 0, 'prophage' => 0, 'pyrro' => 0,'repeat' => 0);

	lives_ok {
		my ($orig_genome) = &get_genome($genome_ref);
		my ($orig_funcs) = &summarize_input($orig_genome);

		print "number of input features = ".scalar  @{$orig_genome->{features}}."\n";
		print "number of input non-coding features = ".scalar  @{$orig_genome->{non_coding_features}}."\n";
		$num_func_in  = scalar  @{$orig_genome->{features}} + @{$orig_genome->{non_coding_features}};

		print("######## Running RAST annotation ########\n");
		my ($genome_obj,$params) = &submit_annotation($genome_obj_name, $genome_ref, $params);

		print "number of returned features = ".scalar  @{$genome_obj->{features}}."\n";
		print "number of returned non-coding features = ".scalar  @{$genome_obj->{non_coding_features}}."\n";
		$num_func_out = scalar  @{$genome_obj->{features}} + @{$genome_obj->{non_coding_features}};

		for (my $i=0; $i < @{$genome_obj->{features}}; $i++) {
			my $ftr = $genome_obj->{features}->[$i];
			my $func      = defined($ftr->{function}) ? $ftr->{function} : "";
			my $orig_func = defined($orig_funcs->{$ftr->{id}}) ? $orig_funcs->{$ftr->{id}} : "";
			if ($func ne $orig_func) {
				$counts{num_diff}++;
			}

			if ($ftr->{ontology_terms}) {
				#print Dumper $ftr->{ontology_terms}->{SSO};
				my @roles = keys(%{$ftr->{ontology_terms}->{SSO}});
				foreach (@roles) {
					$counts{seleno}++ if ($_ =~ /SSO:000009304/);
					$counts{pyrro}++  if ($_ =~ /SSO:000009291/);
				}
			}
	    }

		for (my $i=0; $i < @{$genome_obj->{non_coding_features}}; $i++) {
			my $ftr = $genome_obj->{non_coding_features}->[$i];
			my $func      = defined($ftr->{function}) ? $ftr->{function} : "";
			my $orig_func = defined($orig_funcs->{$ftr->{id}}) ? $orig_funcs->{$ftr->{id}} : "";
			if ($func ne $orig_func) {
				$counts{num_diff}++;
			}
			if ($ftr->{function} &&  $ftr->{type} =~ /crispr/) {
				$counts{crispr}++;
			} elsif ($ftr->{function} &&  $ftr->{function} =~ /phiSpy/) {
				$counts{prophage}++;
			}
			if ($ftr->{type} =~ /repeat/) {
				$counts{repeat}++;
			}
		}

		print "**** Number of features post-annotation = $num_func_out\n";
		$counts{num_new} = $num_func_out - $num_func_in;
	} "Annotate $genome_obj_name";

	print "Summary for $genome_obj_name\n";
	print "Number of features with changed function: " . $counts{num_diff} . "\n";
	return (%counts);
}


done_testing($num_tests);

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


