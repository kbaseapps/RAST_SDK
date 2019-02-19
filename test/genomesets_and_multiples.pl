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
             "call_features_tRNA_trnascan"=>'1',
             "call_selenoproteins"=>'0',
             "call_pyrrolysoproteins"=>'0',
             "call_features_repeat_region_SEED"=>'0',
             "call_features_insertion_sequences"=>'0',
             "call_features_strep_suis_repeat"=>'0',
             "call_features_strep_pneumo_repeat"=>'0',
             "call_features_crispr"=>'0',
             "call_features_CDS_glimmer3"=>'0',
             "call_features_CDS_prodigal"=>'0',
             "annotate_proteins_kmer_v2"=>'0',
             "kmer_v1_parameters"=>'0',
             "annotate_proteins_similarity"=>'0',
             "retain_old_anno_for_hypotheticals"=>'0',
             "resolve_overlapping_features"=>'0',
             "call_features_prophage_phispy"=>'0',
             "output_genome"=>$genome_obj_name,
             "workspace"=>get_ws_name()
           };
    return make_impl_call("RAST_SDK.annotate_genomes", $params);
}

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

my ($orig_genome1,$orig_funcs1) = &get_and_prep($genome_ref1);
my ($orig_genome2,$orig_funcs2) = &get_and_prep($genome_ref2);
my ($orig_genome3,$orig_funcs3) = &get_and_prep($genome_ref3);

#
#	TEST TWO GENOMES -
#
lives_ok {
	if ($DEBUG ne 'Y') {
		my ($genome_obj,$params) = &submit_multi_annotation('multi_genomes', [$genome_ref1,$genome_ref2]);
		my $gs = get_ws_name() . "/" . $genome_obj ;
		my $info = $ws_client->get_objects([{ref=>$gs}])->[0]->{info};
		my $newref = $info->[6]."/".$info->[0]."/".$info->[4];
		# If you got this far, the new object was created.
	}
	1;
} "Pipeline Runs";

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

		my ($genome_obj,$params) = &submit_set_annotation($genome_set_name, [$newref]);
		$info = $ws_client->get_objects([{ref=>$gs}])->[0]->{info};
	# 	If you got this far, the new object was created.
	} else {
		1;
	}
} "Pipeline Runs";

#
#	BUILD AND TEST A GENOME SET PLUS ANOTHER GENOME
#

lives_ok {
	if ($DEBUG eq 'Y') {
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

		my ($genome_obj,$params) = &submit_set_annotation($genome_set_name, [$newref,$genome_ref3]);
		$info = $ws_client->get_objects([{ref=>$gs}])->[0]->{info};
	# 	If you got this far, the new object was created.
	} else {
		1;
	}
} "Pipeline Runs";
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


