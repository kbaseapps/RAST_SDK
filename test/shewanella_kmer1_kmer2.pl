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
use GenomeFileUtil::GenomeFileUtilClient;
use Storable qw(dclone);
use Bio::KBase::kbaseenv;

use lib "/kb/module/test";
use testRASTutil;

my $purpose = "Shewanella_sp._ANA-3 is used as a test for RAST repeats, kmer v2 and kmer v2.\n";
$purpose   .= "Both of these should exist in this genome\n";


local $| = 1;
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('RAST_SDK');
my $ws_url = $config->{"workspace-url"};
my $ws_name = undef;
my $ws_client = new Workspace::WorkspaceClient($ws_url,token => $token);
my $call_back_url = $ENV{ SDK_CALLBACK_URL };
my $au = new AssemblyUtil::AssemblyUtilClient($call_back_url);
my $gfu = new GenomeFileUtil::GenomeFileUtilClient($call_back_url);

sub reannotate_genome {
    my($genome_obj_name, $genome_upa) = @_;
    my $params={"input_genome"=>$genome_upa,
             "call_features_rRNA_SEED"=>'0',
             "call_features_tRNA_trnascan"=>'0',
             "call_selenoproteins"=>'0',
             "call_pyrrolysoproteins"=>'0',
             "call_features_repeat_region_SEED"=>'1',
             "call_features_insertion_sequences"=>'0',
             "call_features_strep_suis_repeat"=>'0',
             "call_features_strep_pneumo_repeat"=>'0',
             "call_features_crispr"=>'0',
             "call_features_CDS_glimmer3"=>'0',
             "call_features_CDS_prodigal"=>'0',
             "annotate_proteins_kmer_v2"=>'1',
             "kmer_v1_parameters"=>'1',
             "annotate_proteins_similarity"=>'0',
             "retain_old_anno_for_hypotheticals"=>'0',
             "resolve_overlapping_features"=>'0',
             "call_features_prophage_phispy"=>'0',
             "output_genome"=>$genome_obj_name,
             "workspace"=>get_ws_name()
           };
    return make_impl_call("RAST_SDK.annotate_genome", $params);
}

my $diff_count   = 0;
my $num_func_in  = 0;
my $num_func_out = 0;
my $genome_obj_name = "Shewanella_sp._ANA-3";
lives_ok {
    my $genome_gbff_name = "Shewanella_sp._ANA-3.gbff";
    my ($tmp_genome_obj, $genome_ref) = prepare_gbff($genome_gbff_name,$genome_obj_name);
#	my $genome_ref = "15792/209997/1";

	my ($orig_genome,$orig_funcs) = &get_and_prep($genome_ref);

	print "number of input features = ".scalar  @{$orig_genome->{features}}."\n";
	print "number of input non-coding features = ".scalar  @{$orig_genome->{non_coding_features}}."\n";
	$num_func_in  = scalar  @{$orig_genome->{features}} + @{$orig_genome->{non_coding_features}};

    print("######## Running RAST annotation ########\n");
	my ($genome_obj,$params) = &submit_annotation($genome_obj_name, $genome_ref);

	print "number of returned features = ".scalar  @{$genome_obj->{features}}."\n";
	print "number of returned non-coding features = ".scalar  @{$genome_obj->{non_coding_features}}."\n";
	$num_func_out = scalar  @{$genome_obj->{features}} + @{$genome_obj->{non_coding_features}};

	my %types = ();
    for (my $i=0; $i < @{$genome_obj->{features}}; $i++) {
        my $ftr = $genome_obj->{features}->[$i];
		if (defined($ftr->{functions}) && scalar @{$ftr->{functions}} > 0){
					$ftr->{function} = join("; ", @{$ftr->{functions}});
		}
        my $func      = defined($ftr->{function}) ? $ftr->{function} : "Not defined";
		my $orig_func = defined($orig_funcs->{$ftr->{id}}) ? $orig_funcs->{$ftr->{id}} : "";
        if ($func ne $orig_func) {
            $diff_count++;
        }
        %types = &count_types($ftr->{type},%types);
    }
    for (my $i=0; $i < @{$genome_obj->{non_coding_features}}; $i++) {
        my $ftr = $genome_obj->{non_coding_features}->[$i];
		if (defined($ftr->{functions}) && scalar @{$ftr->{functions}} > 0){
					$ftr->{function} = join("; ", @{$ftr->{functions}});
		}
        my $func      = defined($ftr->{function}) ? $ftr->{function} : "";
		my $orig_func = defined($orig_funcs->{$ftr->{id}}) ? $orig_funcs->{$ftr->{id}} : "";
        if ($func ne $orig_func) {
            $diff_count++;
        }
        %types = &count_types($ftr->{type},%types);
    }
	print "**** Number of features post-annotation = $num_func_out\n";
    } "Pipeline Runs";
print "Summary for $genome_obj_name\n";
ok($diff_count > 10,"Some of the features changed function ($diff_count)\n");
ok($num_func_out > $num_func_in, "More features out than in.  OUT=$num_func_out IN=$num_func_in\n");
done_testing(3);

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


