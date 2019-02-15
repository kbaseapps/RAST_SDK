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
use installed_clients::kb_SetUtilitiesClient;
use Storable qw(dclone);
use Bio::KBase::kbaseenv;

use lib "/kb/module/test";
use testRASTutil;

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
my $su = new installed_clients::kb_SetUtilitiesClient($call_back_url);

sub reannotate_genomes {
    my($genome_obj_name, $genome_upa) = @_;
    my $params={"input_genomes"=>$genome_upa,
             "call_features_rRNA_SEED"=>'0',
             "call_features_tRNA_trnascan"=>'0',
             "call_selenoproteins"=>'1',
             "call_pyrrolysoproteins"=>'0',
             "call_features_repeat_region_SEED"=>'0',
             "call_features_insertion_sequences"=>'0',
             "call_features_strep_suis_repeat"=>'0',
             "call_features_strep_pneumo_repeat"=>'0',
             "call_features_crispr"=>'0',
             "call_features_CDS_glimmer3"=>'0',
             "call_features_CDS_prodigal"=>'1',
             "annotate_proteins_kmer_v2"=>'0',
             "kmer_v1_parameters"=>'0',
             "annotate_proteins_similarity"=>'0',
             "resolve_overlapping_features"=>'0',
             "call_features_prophage_phispy"=>'0',
             "output_genome"=>$genome_obj_name,
             "workspace"=>get_ws_name()
           };
    return make_impl_call("RAST_SDK.annotate_genomes", $params);
}

sub prepare_assembly {
    my($assembly_obj_name) = @_;
    #my $fasta_data_path = "/kb/module/test/data/Clostridium_thermocellum_ATCC27405.fa";
    my $fasta_data_path = "/kb/module/test/data/$assembly_obj_name";
    my $fasta_temp_path = "/kb/module/work/tmp/$assembly_obj_name";
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


my $DEBUG = 'Y';
my $assembly_obj_name1 = "Dactylopius coccus";
my $assembly_ref1 = "40046/5/1";
   $assembly_ref1 = "40619/11/1";

if ($DEBUG ne 'Y') {
	$assembly_obj_name1 = "bogus.fna";
	$assembly_ref1 = prepare_assembly($assembly_obj_name1);
}

#
#	TEST ONE ASSEMBLY -
#
print "ASSEMBLYREF1 = $assembly_ref1\n";
lives_ok {
		my ($genome_obj,$params) = &submit_multi_annotation('multi_genomes', [$assembly_ref1]);
		my $gs = get_ws_name() . "/" . $genome_obj ;
		my $info = $ws_client->get_objects([{ref=>$gs}])->[0]->{info};
		my $newref = $info->[6]."/".$info->[0]."/".$info->[4];
		# If you got this far, the new object was created.

		my $report = "/kb/module/work/tmp/microbial_genome_report.multi_genomes";
		my $local_path = "/kb/module/test/report_output/microbial_genome_report.multi_genomes";
	    copy $report, $local_path;

		ok(-e $local_path,'File found');
} "Pipeline Runs";

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


