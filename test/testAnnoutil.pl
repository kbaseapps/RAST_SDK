use Test::Most;
use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catfile splitpath);
use File::Copy;
use Carp qw(croak);
use File::Compare;
use Config::Simple;
use Bio::KBase::AuthToken;

use installed_clients::WorkspaceClient;
use installed_clients::GenomeFileUtilClient;
use RAST_SDK::RAST_SDKImpl;

use strict;
use_ok "anno_utils";
use testRASTutil;


## global variables
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('RAST_SDK');
my $auth_token = Bio::KBase::AuthToken->new(
        token => $token, ignore_authrc => 1, auth_svc=>$config->{'auth-service-url'});
my $ws_url = $config->{"workspace-url"};
my $ws = get_ws_name();
my $ws_client = new installed_clients::WorkspaceClient($ws_url,token => $token);
my $call_back_url = $ENV{ SDK_CALLBACK_URL };
my $ctx = LocalCallContext->new($token, $auth_token->user_id);
$RAST_SDK::RAST_SDKServer::CallContext = $ctx;


my $tmp_write_dir = 'data/write_tmp';  ## For saving temporary test files

my $gbff_file = 'data/Clostridium_botulinum.gbff';

my $out_name = 'annotated_metag';
my $outgn_name = 'rast_annotated_genome';
my $fasta1 = 'data/short_one.fa';
my $fasta3 = 'data/GCA_000350285.1_OR1_genomic.fna';
my $fasta4 = 'data/metag_test/Test_v1.0.fa';
my $gff1 = 'data/short_one.gff';
my $fasta2 = 'data/metag_test/59111.assembled.fna';
my $fa_LOng = 'data/LOng_contig_names.fa';
my $gff2 = 'data/metag_test/59111.assembled.gff';
my $fasta_scrt = 'fasta_file.fa';
my $gff_scrt = 'gff_file.gff';

my $rast_impl = new RAST_SDK::RAST_SDKImpl();
my $annoutil = new anno_utils($config, $ctx);

my $scratch = $config->{'scratch'}; #'/kb/module/work/tmp';
my $rast_genome_dir = $annoutil->_create_rast_subdir($scratch, "genome_annotation_dir_");


##-----------------Test Blocks--------------------##

my $ecoli_gff = 'data/metag_test/ecoli_out.gff';
my $ecoli_sco = 'data/metag_test/ecoli_out.sco';
my $trans_file = 'data/metag_test/translationfile';
my %trans_tab;
my $sco_tab = [];

my $obj_Echinacea = "55141/242/1";  #prod genome
my $obj_Echinacea_ann = "55141/247/1";  #prod genome
my $obj_Ecoli = "55141/212/1";  # prod genome
my $obj_Ecoli_ann = "55141/252/1";  # prod genome
my $obj_asmb = "55141/243/1";  # prod assembly
my $obj_asmb_ann = "55141/244/1";  # prod assembly
my $obj_asmb_refseq = "55141/266/3";  # prod assembly
my $obj1 = "37798/14/1";  # appdev
my $obj2 = "37798/15/1";  # appdev
my $obj3 = "55141/77/1";  # prod KBaseGenomeAnnotations.Assembly
my $obj_65386_1 = '65386/2/1';  # same as 63171/436/1, i.e., GCF_003058445.1
my $obj_65386_2 = '65386/12/1';

my $asmb_fasta = $annoutil->_get_fasta_from_assembly($obj_asmb);

#####RAST_SDK Module test objects #####
my $obj2_1 = "63171/315/1";

# used for function lookup in report testing
my $test_ftrs = [{
 'id' => '10000_1',
 'protein_translation' => 'VVVLLHGGCCEDMQRGRRESAPDLTLVVYPHALHALDMRLPDRTVLGMRLGFDAHAAADARRQVLDFLTARGVAPPDR*'
 },
 {
 'function' => 'L-carnitine dehydratase/bile acid-inducible protein F (EC 2.8.3.16)',
 'quality' => {
 'hit_count' => 3
 },
 'annotations' => [
     [
     'Function updated to L-carnitine dehydratase/bile acid-inducible protein F (EC 2.8.3.16)',
     'annotate_proteins_kmer_v1',
     '1583389302.95393',
     '6c670c83-2a11-49ff-97bf-b1c3e2121f30'
     ]
 ],
 'id' => '10000_2'
 },
 {
 'id' => '10000_3',
 'protein_translation' => 'MLAVNEPTVVLASAETKSLGPVVTGDRLETEAEVERTDGRKRWVKVTVRRAGAPVMEGQFLAVVPDRHILDAKDARR*'
 },
 {
 'id' => '10000_madeup',
 'function' => 'completely fake function',
 'annotations' => [
     [
     'completely fake function',
     'annotate_madeup_source',
     '1583389302.95394',
     '6c670c83-2a11-49ff-97bf-b1c3e2121f33'
     ]
 ],
 }];


## Re-mapping the contigIDs back to their original (long) names
subtest '_remap_contigIDs' => sub {
    my $contigID_hash = {
          'contigID_1' => 'NZ_CP028859.1',
          'contigID_2' => 'NZ_CP028860.1'
    };
    my $gn = {
        contigs => [
            {id => 'contigID_1',
             name => 'contigID_1'},
            {id => 'contigID_2',
             name => 'contigID_2'}
        ]
    };
    $gn = $annoutil->_remap_contigIDs($contigID_hash, $gn);
    my $contig1 = $gn->{contigs}[0];
    my $contig2 = $gn->{contigs}[1];
    ok ($contig1->{id} eq $contigID_hash->{contigID_1}, "mapped contigID correctly");
    ok ($contig2->{id} eq $contigID_hash->{contigID_2}, "mapped contigID correctly");
};

subtest '_get_contigs_from_fastafile' => sub {
    my $fa_file = $fa_LOng;
    my $contigs;
    my $contigID_hash;
    lives_ok {
        ($contigs, $contigID_hash) = $annoutil->_get_contigs_from_fastafile($fa_file);
    } '_get_contigs_from_fastafile runs successfully';
    is ($contigID_hash->{contigID_1},
        'CP0035411Candidatus_Carsonella_ruddii_CE_isolate_Thao2000_complete_genome',
        'contigID hash map generated correctly');
    ok ($contigs->[0]->{id} eq 'contigID_1', 'contig id correctly shortened for downstream');
    ok ($contigs->[0]->{name} eq 'contigID_1', 'contig name correctly shortened for downstream');
};

subtest '_get_genome' => sub {
    my $parameters = {
         output_genome_name => 'test_out_gn_name',
         output_workspace => $ws,
         object_ref => $obj_Ecoli
    };

    my $obj_ref = $parameters->{object_ref};
    my $ret_gn;
    lives_ok { 
        $ret_gn = $annoutil->_get_genome($obj_ref);
        # print "genome object returned on $obj_ref:\n".Dumper(keys %$ret_gn);
    } '_get_genome runs successfully';
    ok (@{$ret_gn->{features}} > 0, 'Genome has features!');
    is ($ret_gn->{assembly_ref}, '2901/78/1', 'found genome assembly ref');
};

subtest '_get_contigs' => sub {
    my $parameters = {
         output_genome_name => 'test_out_gn_name',
         output_workspace => $ws,
         object_ref => $obj_Ecoli
    };

    my $obj_ref = $parameters->{object_ref};
    my $obj = $annoutil->_fetch_object_data($obj_ref);
    
    lives_ok { 
        my $contig_obj1 = $annoutil->_get_contigs($obj->{assembly_ref});
        print "anno_tuils _get_contigs returns:\n".Dumper(keys %$contig_obj1);
    } '_get_contigs runs successfully on genome';

    lives_ok { 
        my $contig_obj2 = $annoutil->_get_contigs($obj_asmb);
        print "anno_tuils _get_contigs returns:\n".Dumper(keys %$contig_obj2);
    } '_get_contigs runs successfully on assembly';

};


## Testing the feature_function_lookup related functions
subtest '_get_feature_function_lookup' => sub {
    my %ffunc_lookup = ();
    lives_ok { 
        %ffunc_lookup = $annoutil->_get_feature_function_lookup($test_ftrs);
    } '_get_feature_function_lookup runs successfully on assembly';
    ok (exists($ffunc_lookup{'10000_2'}), 'found one key'); 
    ok (exists($ffunc_lookup{'10000_madeup'}), 'found another key'); 
    ok (exists($ffunc_lookup{'10000_2'}->{functions}), 'found functions in one'); 
    ok (exists($ffunc_lookup{'10000_madeup'}->{functions}), 'found functions in another'); 
    ok (exists($ffunc_lookup{'10000_2'}->{annotation_src}), 'found annotation_src in one'); 
    ok (exists($ffunc_lookup{'10000_madeup'}->{annotation_src}),
                                'found annotation_src in another');

    my $func_role = 'completely fake function';
    my $exp_src1 = 'annotate_madeup_source';
    my $ann_src1 = $annoutil->_find_function_source(\%ffunc_lookup, $func_role);
    is ($ann_src1, $exp_src1, "Found function $func_role with annotation source of: $ann_src1"); 

    $func_role = 'L-carnitine dehydratase/bile acid-inducible protein';
    my $exp_src2 = 'annotate_proteins_kmer_v1';
    my $ann_src2 = $annoutil->_find_function_source(\%ffunc_lookup, $func_role);
    is ($ann_src2, $exp_src2, "Found function $func_role with annotation source of: $ann_src2"); 

    $func_role = 'non-existent function';
    my $exp_src3 = 'N/A';
    my $ann_src3 = $annoutil->_find_function_source(\%ffunc_lookup, $func_role);
    is ($ann_src3, $exp_src3, "Found function $func_role with annotation source of: $ann_src3"); 
};


#
## Global variables for the annotation process steps to share ##
#
my ($rast_ref, %rast_details01, %rast_details02, %rast_details1, %rast_details2);
my ($inputgenome01, $inputgenome02, $inputgenome1, $inputgenome2);
my ($parameters01, $parameters02, $parameters1, $parameters2);

# Test _set_parameters_by_input with genome/assembly object refs in prod
subtest '_set_parameters_by_input' => sub {
    # a genome object in workspace #65386
    $parameters01 = {
         output_genome_name => 'test_out_gn_name01',
         scientific_name => 'Clostridium botulinum',
         output_workspace => $ws,
         object_ref => $obj_65386_1
    };
    # 01. creating default genome object
    $inputgenome01 = {
        id => $parameters01->{output_genome_name},
        genetic_code => $parameters01->{genetic_code},
        scientific_name => $parameters01->{scientific_name},
        domain => $parameters01->{domain},
        contigs => [],
        features => []
    };

    if ($parameters01->{ncbi_taxon_id}) {
        $inputgenome01->{taxon_assignments} = {
            'ncbi' => '' . $parameters01->{ncbi_taxon_id}};
    }
    my $expected_params01 = {
          scientific_name => $parameters01->{scientific_name},
          output_genome_name => $parameters01->{output_genome_name},
          output_workspace => $ws,
          object_ref => $parameters01->{object_ref},
          genetic_code => 11,
          domain => 'Bacteria'
    };

    lives_ok {
        ($rast_ref, $inputgenome01) = $annoutil->_set_parameters_by_input(
                                            $parameters01, $inputgenome01);
    } "_set_parameters_by_input runs successfully on genome $obj_65386_1";
    %rast_details01 = %{ $rast_ref }; # dereference
    $parameters01 = $rast_details01{parameters};
    ok (@{$inputgenome01->{features}} > 0,
        "inputgenome has ".scalar @{$inputgenome01->{features}}." feature(s).");
    ok (@{$rast_details01{contigobj}{contigs}} > 0,
        "inputgenome has ".scalar @{$rast_details01{contigobj}{contigs}} ." contig(s).");
    is ($inputgenome01->{assembly_ref}, '19217/360049/1', 'inputgenome assembly_ref is correct.');

    # before merging with default gene call settings
    cmp_deeply($parameters01, $expected_params01, "parameters has expected input param values.");

    #  merge with the default gene call settings
    my $default_params = $annoutil->_set_default_parameters();
    $parameters01 = { %$default_params, %$parameters01 };
    $rast_details01{parameters} = $parameters01;

    my $expected_params011 = { %$default_params, %$expected_params01 };
    # after merging with default gene call settings
    cmp_deeply($parameters01, $expected_params011, "parameters has default workflows.");
    my $expected_contigID_hash = {
          'contigID_1' => 'NZ_CP028859.1',
          'contigID_2' => 'NZ_CP028860.1'
    };
    cmp_deeply($rast_details01{contigID_hash},
               $expected_contigID_hash, "contigID hash was generated correctly.");

    # another genome object in workspace #65386 that did not have a problem to run
    $parameters02 = {
         output_genome_name => 'test_out_gn_name02',
         scientific_name => 'Methanosarcina acetivorans C2A',
         output_workspace => $ws,
         object_ref => $obj_65386_2
    };

    # 02. creating default genome object
    $inputgenome02 = {
        id => $parameters02->{output_genome_name},
        genetic_code => $parameters02->{genetic_code},
        scientific_name => $parameters02->{scientific_name},
        domain => $parameters02->{domain},
        contigs => [],
        features => []
    };

    if ($parameters02->{ncbi_taxon_id}) {
        $inputgenome02->{taxon_assignments} = {
            'ncbi' => '' . $parameters02->{ncbi_taxon_id}};
    }
    my $expected_params02 = {
          scientific_name => $parameters02->{scientific_name},
          output_genome_name => $parameters02->{output_genome_name},
          output_workspace => $ws,
          object_ref => $parameters02->{object_ref},
          genetic_code => 11,
          domain => 'Archaea'
    };

    lives_ok {
        ($rast_ref, $inputgenome02) = $annoutil->_set_parameters_by_input(
                                            $parameters02, $inputgenome02);
    } "_set_parameters_by_input runs successfully on genome $obj_65386_2";
    %rast_details02 = %{ $rast_ref }; # dereference
    $parameters02 = $rast_details02{parameters};
    ok (@{$inputgenome02->{features}} > 0,
        "inputgenome has ".scalar @{$inputgenome02->{features}}." feature(s).");
    ok (@{$rast_details02{contigobj}{contigs}} > 0,
        "inputgenome has ".scalar @{$rast_details02{contigobj}{contigs}} ." contig(s).");
    is ($inputgenome02->{assembly_ref}, '19217/194865/2', 'inputgenome assembly_ref is correct.');

    # before merging with default gene call settings
    cmp_deeply($parameters02, $expected_params02, "parameters has expected input param values.");

    #  merge with the default gene call settings
    $parameters02 = { %$default_params, %$parameters02 };
    $rast_details02{parameters} = $parameters02;
    my $expected_params022 = { %$default_params, %$expected_params02};

    # after merging with default gene call settings
    cmp_deeply($parameters02, $expected_params022, "parameters has default workflows.");

    # a genome object
    $parameters1 = {
         output_genome_name => 'test_out_gn_name1',
         output_workspace => $ws,
         object_ref => $obj_Ecoli
    };
    # 1. creating default genome object
    $inputgenome1 = {
        id => $parameters1->{output_genome_name},
        genetic_code => $parameters1->{genetic_code},
        scientific_name => $parameters1->{scientific_name},
        domain => $parameters1->{domain},
        contigs => [],
        features => []
    };

    if ($parameters1->{ncbi_taxon_id}) {
        $inputgenome1->{taxon_assignments} = {
            'ncbi' => '' . $parameters1->{ncbi_taxon_id}};
    }
    my $expected_params1 = {
          'scientific_name' => 'Escherichia coli str. K-12 substr. MG1655',
          'output_genome_name' => $parameters1->{output_genome_name},
          'output_workspace' => $ws,
          'object_ref' => $parameters1->{object_ref},
          'genetic_code' => 11,
          'domain' => 'Bacteria'
    };

    lives_ok {
        ($rast_ref, $inputgenome1) = $annoutil->_set_parameters_by_input(
                                            $parameters1, $inputgenome1);
    } "_set_parameters_by_input runs successfully on genome $obj_Ecoli";
    %rast_details1 = %{ $rast_ref }; # dereference
    $parameters1 = $rast_details1{parameters};

    ok (@{$inputgenome1->{features}} > 0,
        "inputgenome has ".scalar @{$inputgenome1->{features}}." feature(s).");
    ok (@{$rast_details1{contigobj}{contigs}} > 0,
        "inputgenome has ".scalar @{$rast_details1{contigobj}{contigs}} ." contig(s).");
    is ($inputgenome1->{assembly_ref}, '2901/78/1', 'inputgenome assembly_ref is correct.');

    # before merging with default gene call settings
    cmp_deeply($parameters1, $expected_params1, 'parameters has expected input param values.');

    #  merge with the default gene call settings
    $default_params = $annoutil->_set_default_parameters();
    $parameters1 = { %$default_params, %$parameters1 };
    $rast_details1{parameters} = $parameters1;

    my $expected_params2 = { %$default_params, %$expected_params1 };

    # after merging with default gene call settings
    cmp_deeply($expected_params2, $parameters1, 'parameters has default workflows.');

    # an assembly object
    $parameters2 = {
         output_genome_name => 'test_out_gn_name2',
         output_workspace => $ws,
         object_ref => $obj_asmb
    };
    # 2 creating default genome object
    $inputgenome2 = {
        id => $parameters2->{output_genome_name},
        genetic_code => $parameters2->{genetic_code},
        scientific_name => $parameters2->{scientific_name},
        domain => $parameters2->{domain},
        contigs => [],
        features => []
    };

    if ($parameters2->{ncbi_taxon_id}) {
        $inputgenome2->{taxon_assignments} = {
            'ncbi' => '' . $parameters2->{ncbi_taxon_id}};
    }

    my $expected_params3 = {
          'object_ref' => $obj_asmb,
          'genetic_code' => undef,
          'output_genome_name' => $parameters2->{output_genome_name},
          'domain' => undef,
          'scientific_name' => undef,
          'output_workspace' => $ws
    };
    lives_ok {
        ($rast_ref, $inputgenome2) = $annoutil->_set_parameters_by_input(
                                            $parameters2, $inputgenome2);
    } "_set_parameters_by_input runs successfully on assembly $obj_asmb";
    %rast_details2 = %{$rast_ref}; # dereference
    $parameters2 = $rast_details2{parameters};

    # before merging with default gene call settings
    cmp_deeply($parameters2, $expected_params3, 'parameters has input values');

    # merge with the default gene call settings
    $parameters2 = { %$default_params, %$parameters2 };
    my $expected_params4 = { %$default_params, %$expected_params3 };
    $rast_details2{parameters} = $parameters2;

    # after merging with default gene call settings
    cmp_deeply($expected_params4, $parameters2, 'parameters has default workflows');

    ok (@{$inputgenome2->{features}} == 0, 'inputgenome (assembly) has NO features.');
    ok (@{$rast_details2{contigobj}{contigs}} > 0,
        "inputgenome has ".scalar @{$rast_details2{contigobj}{contigs}} ." contig(s).");
    is ($inputgenome2->{assembly_ref}, $obj_asmb, 'inputgenome assembly_ref is correct.');
};


# Test _set_messageNcontigs with genome/assembly object refs in prod
subtest '_set_messageNcontigs' => sub {
    # a genome object in workspace #65386
    lives_ok {
        ($rast_ref, $inputgenome01) = $annoutil->_set_messageNcontigs(
                                            \%rast_details01, $inputgenome01);
    } "_set_messageNcontigs runs successfully on genome $obj_65386_1";
    %rast_details01 = %{ $rast_ref }; # dereference
    $parameters01 = $rast_details01{parameters};
    my $msg01 = $rast_details01{message};
    my $tax01 = $rast_details01{tax_domain};

    ok (@{$inputgenome01->{features}} > 0,
        "inputgenome has ".scalar @{$inputgenome01->{features}}. " feature(s).");
    ok (@{$inputgenome01->{contigs}} > 0,
        "inputgenome has ".scalar @{$inputgenome01->{contigs}} . " contig(s).");
    ok (length($msg01) > 0, "Message for genome input has contents:\n$msg01");
    ok ($tax01 eq 'Bacteria', "tax_domain for genome input has value:$tax01");

    # another genome object in workspace #65386
    lives_ok {
        ($rast_ref, $inputgenome02) = $annoutil->_set_messageNcontigs(
                                            \%rast_details02, $inputgenome02);
    } "_set_messageNcontigs runs successfully on genome $obj_65386_2";

    %rast_details02 = %{ $rast_ref }; # dereference
    $parameters02 = $rast_details02{parameters};
    my $msg02 = $rast_details02{message};
    my $tax02 = $rast_details02{tax_domain};

    ok (@{$inputgenome02->{features}} > 0,
        "inputgenome has ".scalar @{$inputgenome02->{features}}. " feature(s).");
    ok (@{$inputgenome02->{contigs}} > 0,
        "inputgenome has ".scalar @{$inputgenome02->{contigs}} . " contig(s).");
    ok (length($msg02) > 0, "Message for genome input has contents:\n$msg02");
    ok ($tax02 eq 'Archaea', "tax_domain for genome input has value:$tax02");

    # a genome object
    lives_ok {
        ($rast_ref, $inputgenome1) = $annoutil->_set_messageNcontigs(
                                            \%rast_details1, $inputgenome1);
    } "_set_messageNcontigs runs successfully on genome $obj_Ecoli";
    %rast_details1 = %{ $rast_ref }; # dereference
    $parameters1 = $rast_details1{parameters};
    my $msg1 = $rast_details1{message};
    my $tax1 = $rast_details1{tax_domain};

    ok (@{$inputgenome1->{features}} > 0,
        "inputgenome has ".scalar @{$inputgenome1->{features}}. " feature(s).");
    ok (@{$inputgenome1->{contigs}} > 0,
        "inputgenome has ".scalar @{$inputgenome1->{contigs}} . " contig(s).");
    ok (length($msg1) > 0, "Message for genome input has contents:\n$msg1");
    ok ($tax1 eq 'Bacteria', "tax_domain for genome input has value:$tax1");

    # an assembly object
    lives_ok {
        ($rast_ref, $inputgenome2) = $annoutil->_set_messageNcontigs(
                                            \%rast_details2, $inputgenome2);
    } "_set_messageNcontigs runs successfully on assembly $obj_asmb";
    %rast_details2 = %{$rast_ref}; # dereference
    $parameters2 = $rast_details2{parameters};
    my $msg2 = $rast_details2{message};
    my $tax2 = $rast_details2{tax_domain};

    ok (@{$inputgenome2->{features}} == 0, 'inputgenome (assembly) has no features.');
    ok (@{$inputgenome2->{contigs}} > 0,
        "inputgenome has ".scalar @{$inputgenome2->{contigs}} . " contig(s).");
    ok (length($msg2) > 0, "Message for assembly input has contents.");
    ok ($tax2 eq 'U', "tax_domain for assembly input has value:$tax2");
};


# Test _set_genecall_workflow with genome/assembly object refs in prod
subtest '_set_genecall_workflow' => sub {
    my $exp_gc_workflow = {
        'stages' => [
            { 'name' => 'call_features_CDS_prodigal' },
            { 'name' => 'call_features_CDS_glimmer3',
              'glimmer3_parameters' => {
                                         'min_training_len' => '2000'
                                       }
            },
            {
              'name' => 'call_features_rRNA_SEED'
            },
            {
              'name' => 'call_features_tRNA_trnascan'
            },
            {
              'name' => 'call_selenoproteins'
            },
            {
              'name' => 'call_pyrrolysoproteins'
            },
            { 'name' => 'call_features_repeat_region_SEED',
              'repeat_region_SEED_parameters' => {
                                                   'min_length' => '100',
                                                   'min_identity' => '95'
                                                 }
            },
            { 'name' => 'call_features_crispr' }
        ]
    };

    # a genome object in workspace #65386
    lives_ok {
        $rast_ref = $annoutil->_set_genecall_workflow(
                                            \%rast_details01, $inputgenome01);
    } "_set_genecall_workflow runs successfully on genome $obj_65386_1";
    %rast_details01 = %{ $rast_ref }; # dereference
    my $genecall_workflow01 = $rast_details01{genecall_workflow};
    cmp_deeply($genecall_workflow01, $exp_gc_workflow, 'gc_workflow built correctly');

    # another genome object in workspace #65386
    lives_ok {
        $rast_ref = $annoutil->_set_genecall_workflow(
                                            \%rast_details02, $inputgenome02);
    } "_set_genecall_workflow runs successfully on genome $obj_65386_2";
    %rast_details02 = %{ $rast_ref }; # dereference
    my $genecall_workflow02 = $rast_details02{genecall_workflow};
    cmp_deeply($genecall_workflow02, $exp_gc_workflow, 'gc_workflow built correctly');

    # a genome object
    lives_ok {
        $rast_ref = $annoutil->_set_genecall_workflow(
                                            \%rast_details1, $inputgenome1);
    } "_set_genecall_workflow runs successfully on genome $obj_Ecoli";
    %rast_details1 = %{ $rast_ref }; # dereference
    my $genecall_workflow1 = $rast_details1{genecall_workflow};
    cmp_deeply($genecall_workflow1, $exp_gc_workflow, 'gc_workflow built correctly');

    # an assembly object
    my $exp_gc_workflow2 = {
        'stages' => [
            { 'name' => 'call_features_CDS_prodigal'},
            { 'name' => 'call_features_CDS_glimmer3',
              'glimmer3_parameters' => {
                                           'min_training_len' => '2000'
                                       }
            },
            { 'name' => 'call_features_rRNA_SEED' },
            { 'name' => 'call_features_tRNA_trnascan' },
            { 'name' => 'call_features_repeat_region_SEED',
              'repeat_region_SEED_parameters' => {
                                                    'min_length' => '100',
                                                    'min_identity' => '95'
                                                 }
            }
        ]
    };

    lives_ok {
        $rast_ref = $annoutil->_set_genecall_workflow(
                                            \%rast_details2, $inputgenome2);
    } "_set_genecall_workflow runs successfully on assembly $obj_asmb";
    %rast_details2 = %{$rast_ref}; # dereference
    my $msg2 = $rast_details2{message};
    my $genecall_workflow2 = $rast_details2{genecall_workflow};
    cmp_deeply($genecall_workflow2, $exp_gc_workflow2, 'gc_workflow built correctly');
};


# Test _set_annotation_workflow with genome/assembly object refs in prod
subtest '_set_annotation_workflow' => sub {
    my $exp_ann_workflow = {
        'stages' => [
            { 'name' => 'annotate_proteins_kmer_v2',
                'kmer_v2_parameters' => {
                                          'min_hits' => '5',
                                          'annotate_hypothetical_only' => 1
                                        }
            },
            { 'kmer_v1_parameters' => {
                                        'annotate_hypothetical_only' => 1,
                                        'dataset_name' => 'Release70'
                                      },
              'name' => 'annotate_proteins_kmer_v1'
            },
            { 'name' => 'annotate_proteins_similarity',
              'similarity_parameters' => {
                                             'annotate_hypothetical_only' => 1
                                           }
            },
            { 'resolve_overlapping_features_parameters' => {},
              'name' => 'resolve_overlapping_features'
            },
            {
              'name' => 'call_features_prophage_phispy'
            }
        ]
    };

    # a genome object in workspace #65386
    lives_ok {
        $rast_ref = $annoutil->_set_annotation_workflow(\%rast_details01);
    } "_set_genecall_workflow runs successfully on genome $obj_65386_1";
    %rast_details01 = %{ $rast_ref }; # dereference
    my $annomessage01 = $rast_details01{annomessage};
    my $annotate_workflow01 = $rast_details01{annotate_workflow};
    print "genome annotation-msg:\n$annomessage01";
    cmp_deeply($annotate_workflow01, $exp_ann_workflow, 'ann_workflow built correctly');

    # another genome object in workspace #65386
    lives_ok {
        $rast_ref = $annoutil->_set_annotation_workflow(\%rast_details02);
    } "_set_genecall_workflow runs successfully on genome $obj_65386_2";
    %rast_details02 = %{ $rast_ref }; # dereference
    my $annomessage02 = $rast_details02{annomessage};
    my $annotate_workflow02 = $rast_details02{annotate_workflow};
    print "genome annotation-msg:\n$annomessage02";
    cmp_deeply($annotate_workflow02, $exp_ann_workflow, 'ann_workflow built correctly');

    # a genome object
    lives_ok {
        $rast_ref = $annoutil->_set_annotation_workflow(\%rast_details1);
    } "_set_genecall_workflow runs successfully on genome $obj_Ecoli";
    %rast_details1 = %{ $rast_ref }; # dereference
    my $annomessage1 = $rast_details1{annomessage};
    my $annotate_workflow1 = $rast_details1{annotate_workflow};
    print "genome annotation-msg:\n$annomessage1";
    cmp_deeply($annotate_workflow1, $exp_ann_workflow, 'ann_workflow built correctly');

    # an assembly object
    my $exp_ann_workflow2 = {
        'stages' => [
            { 'name' => 'annotate_proteins_kmer_v2',
                'kmer_v2_parameters' => {
                                          'min_hits' => '5',
                                          'annotate_hypothetical_only' => 1
                                        }
            },
            { 'kmer_v1_parameters' => {
                                        'annotate_hypothetical_only' => 1,
                                        'dataset_name' => 'Release70'
                                      },
              'name' => 'annotate_proteins_kmer_v1'
            },
            { 'name' => 'annotate_proteins_similarity',
              'similarity_parameters' => {
                                             'annotate_hypothetical_only' => 1
                                           }
            },
            { 'resolve_overlapping_features_parameters' => {},
              'name' => 'resolve_overlapping_features'
            }
        ]
    };

    lives_ok {
        $rast_ref = $annoutil->_set_annotation_workflow(\%rast_details2);
    } "_set_genecall_workflow runs successfully on assembly $obj_asmb";
    %rast_details2 = %{$rast_ref}; # dereference
    my $annomessage2 = $rast_details2{annomessage};
    my $annotate_workflow2 = $rast_details2{annotate_workflow};
    print "assembly annotation-msg:\n$annomessage2";
    cmp_deeply($exp_ann_workflow2, $annotate_workflow2, 'ann_workflow built correctly');
};


# Test _renumber_features with genome/assembly object refs in prod
subtest '_renumber_features' => sub {
    # a genome object in workspace #65386
    lives_ok {
        ($rast_ref, $inputgenome01) = $annoutil->_renumber_features(
                        \%rast_details01, $inputgenome01);
    } "_renumber_features runs successfully on genome $obj_65386_1";
    %rast_details01 = %{ $rast_ref }; # dereference
    my $msg01 = $rast_details01{message};
    print "genome merged-msg:\n$msg01";
    ok (@{$inputgenome01->{features}} == 0,
        "_renumber_features: inputgenome has NO features.");

    # another genome object in workspace #65386
    lives_ok {
        ($rast_ref, $inputgenome02) = $annoutil->_renumber_features(
                        \%rast_details02, $inputgenome02);
    } "_renumber_features runs successfully on genome $obj_65386_2";
    %rast_details02 = %{ $rast_ref }; # dereference
    my $msg02 = $rast_details02{message};
    print "genome merged-msg:\n$msg02";
    ok (@{$inputgenome02->{features}} == 0,
        "_renumber_features: inputgenome has NO features.");

    # a genome object
    lives_ok {
        ($rast_ref, $inputgenome1) = $annoutil->_renumber_features(
                        \%rast_details1, $inputgenome1);
    } "_renumber_features runs successfully on genome $obj_Ecoli";
    %rast_details1 = %{ $rast_ref }; # dereference
    my $msg1 = $rast_details1{message};
    print "genome merged-msg:\n$msg1";
    ok (@{$inputgenome1->{features}} > 0,
        "_renumber_features: inputgenome has ".scalar @{$inputgenome1->{features}}." features.");

    # an assembly object
    lives_ok {
        ($rast_ref, $inputgenome2) = $annoutil->_renumber_features(
                        \%rast_details2, $inputgenome2);
    } "_renumber_features runs successfully on assembly $obj_asmb";
    %rast_details2 = %{$rast_ref}; # dereference
    ok (exists($rast_details2{extra_workflow}), "extra_workflow is created");
    ok (@{$inputgenome2->{features}} == 0,
        "_renumber_features: Assembly inputgenome has NO features.");
};


# Test _pre_rast_call with genome/assembly object refs in prod
subtest '_pre_rast_call' => sub {
    # a genome object in workspace #65386
    lives_ok {
        ($rast_ref, $inputgenome01) = $annoutil->_pre_rast_call(
                        \%rast_details01, $inputgenome01);
    } "_pre_rast_call runs successfully on genome $obj_65386_1";
    %rast_details01 = %{ $rast_ref }; # dereference
    my $genehash01 = $rast_details01{genehash};
    ok (keys %{$genehash01}, "Gene hash created from genome with elements.");
    if (defined($inputgenome01->{ontology_events})){
        ok (@{$inputgenome01->{ontology_events}} > 0,
            "There are ".scalar @{$inputgenome01->{ontology_events}}." ontology events for genome.");
    }
    ok (@{$inputgenome01->{features}} == 0,
        "_pre_rast_call: inputgenome has NO feature(s).");
    ok (@{$rast_details01{contigobj}{contigs}} == 2,
        "_pre_rast_call: inputgenome has 2 contig(s).");

    # another genome object in workspace #65386
    lives_ok {
        ($rast_ref, $inputgenome02) = $annoutil->_pre_rast_call(
                        \%rast_details02, $inputgenome02);
    } "_pre_rast_call runs successfully on genome $obj_65386_2";
    %rast_details02 = %{ $rast_ref }; # dereference
    my $genehash02 = $rast_details02{genehash};

    ok (keys %{$genehash02}, "Gene hash created from genome with elements.");
    if (defined($inputgenome02->{ontology_events})){
        ok (@{$inputgenome02->{ontology_events}} > 0,
            "There are ".scalar @{$inputgenome02->{ontology_events}}." ontology events for genome.");
    }
    ok (@{$inputgenome02->{features}} == 0,
        "_pre_rast_call: inputgenome has NO feature(s).");
    ok (@{$rast_details02{contigobj}{contigs}} > 0,
        "_pre_rast_call: inputgenome has 2 contig(s).");

    # a genome object
    lives_ok {
        ($rast_ref, $inputgenome1) = $annoutil->_pre_rast_call(
                        \%rast_details1, $inputgenome1);
    } "_pre_rast_call runs successfully on genome $obj_Ecoli";
    %rast_details1 = %{ $rast_ref }; # dereference
    my $genehash1 = $rast_details1{genehash};

    ok (keys %{$genehash1}, "Gene hash created from genome with elements.");
    if (defined($inputgenome1->{ontology_events})){
        ok (@{$inputgenome1->{ontology_events}} > 0,
            "There are ".scalar @{$inputgenome1->{ontology_events}}." ontology events for genome.");
    }
    ok (@{$inputgenome1->{features}} > 0,
        "_pre_rast_call: inputgenome has ".scalar @{$inputgenome1->{features}}." feature(s).");
    ok (@{$rast_details1{contigobj}{contigs}} > 0,
        "_pre_rast_call: inputgenome has ".scalar @{$rast_details1{contigobj}{contigs}} ." contig(s).");

    # an assembly object
    lives_ok {
        ($rast_ref, $inputgenome2) = $annoutil->_pre_rast_call(
                        \%rast_details2, $inputgenome2);
    } "_pre_rast_call runs successfully on assembly $obj_asmb";
    %rast_details2 = %{$rast_ref}; # dereference
    my $genehash2 = $rast_details2{genehash};

    ok (keys %{$genehash2} == 0, "Gene hash created from assembly with no elements.");
    if (defined($inputgenome2->{ontology_events})){
        ok (@{$inputgenome2->{ontology_events}} > 0,
            "There are ".scalar @{$inputgenome2->{ontology_events}}." ontology events for assembly.");
    }
    ok (@{$inputgenome2->{features}} == 0,
        "_pre_rast_call: inputgenome has NO feature(s).");
    ok (@{$rast_details2{contigobj}{contigs}} > 0,
        "_pre_rast_call: inputgenome has ".scalar @{$rast_details2{contigobj}{contigs}} ." contig(s).");
};


# Test _run_rast_workflow on annotation workflow with object refs in prod
my ($ann_genome01, $ann_genome02, $ann_genome1, $ann_genome2,
    $final_genome01, $final_genome02,$final_genome1, $final_genome2);
subtest '_run_rast_workflow_ann' => sub {
    # a genome object in workspace #65386
    lives_ok {
        $ann_genome01 = $annoutil->_run_rast_workflow(
              $inputgenome01, $rast_details01{annotate_workflow});
    } "_run_rast_workflow on annotation runs successfully on genome $obj_65386_1";

    # another genome object in workspace #65386
    lives_ok {
        $ann_genome02 = $annoutil->_run_rast_workflow(
              $inputgenome02, $rast_details02{annotate_workflow});
    } "_run_rast_workflow on annotation runs successfully on genome $obj_65386_2";

    # a genome object
    lives_ok {
        $ann_genome1 = $annoutil->_run_rast_workflow(
              $inputgenome1, $rast_details1{annotate_workflow});
    } "_run_rast_workflow on annotation runs successfully on genome $obj_Ecoli";

    # an assembly object
    lives_ok {
        $ann_genome2 = $annoutil->_run_rast_workflow(
              $inputgenome2, $rast_details2{annotate_workflow});
    } "_run_rast_workflow returns normally on assembly $obj_asmb";
};


# Test _post_rast_ann_call with genome/assembly object refs in prod
subtest '_post_rast_ann_call' => sub {
    # a genome object in workspace #65386
    my $cnt = 0;
    my $ncoding_features01 = $ann_genome01->{non_coding_features};
    for my $ncoding_ftr (@{$ncoding_features01}) {
        if(exists($ncoding_ftr->{type})) {
            $cnt++;
        }
    }
    print "**for $obj_65386_1: Total count of non_coding_features BEFORE _post_rast_ann_call:".scalar @{$ncoding_features01}."\n";
    print "**for $obj_65386_1:Count of non_coding_features WITH 'type' BEFORE _post_rast_ann_call:$cnt\n";
    lives_ok {
        $final_genome01 = $annoutil->_post_rast_ann_call(
                          $ann_genome01, $inputgenome01,
                          $rast_details01{parameters},
                          $rast_details01{contigobj});
    } "_post_rast_ann_call runs successfully on genome $obj_65386_1";
    # Focus on the 'non_coding_features' type field check
    $cnt = 0;
    $ncoding_features01 = $final_genome01->{non_coding_features};
    for my $ncoding_ftr (@{$ncoding_features01}) {
        if(exists($ncoding_ftr->{type})) {
            $cnt++;
            # print "type value: $ncoding_ftr->{type}\n";
        }
    }
    print "**for $obj_65386_1:Count of non_coding_features WITH 'type' AFTER _post_rast_ann_call:$cnt\n";

    # another genome object in workspace #65386
    $cnt = 0;
    my $ncoding_features02 = $ann_genome02->{non_coding_features};
    for my $ncoding_ftr (@{$ncoding_features02}) {
        if(exists($ncoding_ftr->{type})) {
            $cnt++;
        }
    }
    print "**for $obj_65386_2:Total count of non_coding_features BEFORE _post_rast_ann_call:".scalar @{$ncoding_features02}."\n";
    print "**for $obj_65386_2:Count of non_coding_features WITH 'type' BEFORE _post_rast_ann_call:$cnt\n";
    lives_ok {
        $final_genome02 = $annoutil->_post_rast_ann_call(
                          $ann_genome02, $inputgenome02,
                          $rast_details02{parameters},
                          $rast_details02{contigobj});
    } "_post_rast_ann_call runs successfully on genome $obj_65386_2";
    # Focus on the 'non_coding_features' type field check
    $cnt = 0;
    $ncoding_features02 = $final_genome02->{non_coding_features};
    for my $ncoding_ftr (@{$ncoding_features02}) {
        if(exists($ncoding_ftr->{type})) {
            $cnt++;
        }
    }
    print "**for $obj_65386_2:Count of non_coding_features WITH 'type' AFTER _post_rast_ann_call:$cnt\n";

    # a genome object
    lives_ok {
        $final_genome1 = $annoutil->_post_rast_ann_call(
                          $ann_genome1, $inputgenome1,
                          $rast_details1{parameters},
                          $rast_details1{contigobj});
    } "_post_rast_ann_call runs successfully on genome $obj_Ecoli";

    # an assembly object
    lives_ok {
        $final_genome2 = $annoutil->_post_rast_ann_call(
                          $ann_genome2, $inputgenome2,
                          $rast_details2{parameters},
                          $rast_details2{contigobj});
    } "_post_rast_ann_call runs successfully on assembly $obj_asmb";
};


# Test _build_seed_ontology with genome/assembly object refs in prod
subtest '_build_seed_ontology' => sub {
    # a genome object in workspace #65386
    my $nc_ftr_count = @{$final_genome01->{non_coding_features}};
    print "\n********For case $obj_65386_1*********\nBEFORE _build_seed_ontology there are $nc_ftr_count non_coding features.\n";

    lives_ok {
        ($final_genome01, $rast_ref) = $annoutil->_build_seed_ontology(
              \%rast_details01, $final_genome01, $inputgenome01);
    } "_build_seed_ontology returns normally on genome $obj_65386_1";
    %rast_details01 = %{$rast_ref}; # dereference
    $nc_ftr_count = @{$final_genome01->{non_coding_features}};
    print "\n********For case $obj_65386_1*********\nAFTER _build_seed_ontology there are $nc_ftr_count non_coding features.\n";

    # another genome object in workspace #65386
    $nc_ftr_count = @{$final_genome02->{non_coding_features}};
    print "\n********For case $obj_65386_2*********\nBEFORE _build_seed_ontology there are $nc_ftr_count non_coding features.\n";

    lives_ok {
        ($final_genome02, $rast_ref) = $annoutil->_build_seed_ontology(
              \%rast_details02, $final_genome02, $inputgenome02);
    } "_build_seed_ontology returns normally on genome $obj_65386_2";
    %rast_details02 = %{$rast_ref}; # dereference
    $nc_ftr_count = @{$final_genome02->{non_coding_features}};
    print "\n********For case $obj_65386_2*********\nAFTER _build_seed_ontology there are $nc_ftr_count non_coding features.\n";

    # a genome object
    lives_ok {
        ($final_genome1, $rast_ref) = $annoutil->_build_seed_ontology(
              \%rast_details1, $final_genome1, $inputgenome1);
    } "_build_seed_ontology returns normally on genome $obj_Ecoli";
    %rast_details1 = %{$rast_ref}; # dereference

    # an assembly object
    lives_ok {
        ($final_genome2, $rast_ref) = $annoutil->_build_seed_ontology(
              \%rast_details2, $final_genome2, $inputgenome2);
    } "_build_seed_ontology returns on assembly $obj_asmb";
    %rast_details2 = %{$rast_ref}; # dereference
};


# Test _summarize_annotation with genome/assembly object refs in prod
subtest '_summarize_annotation' => sub {
    # a genome object in workspace #65386
    my $nc_ftr_count = @{$final_genome01->{non_coding_features}};
    print "\n********For case $obj_65386_1*********\nBEFORE _summarize_annotation there are $nc_ftr_count non_coding features.\n";

    lives_ok {
        ($final_genome01, $rast_ref) = $annoutil->_summarize_annotation(
              \%rast_details01, $final_genome01, $inputgenome01);
    } "_summarize_annotation runs successfully on genome $obj_65386_1";
    %rast_details01 = %{$rast_ref}; # dereference

    $nc_ftr_count = @{$final_genome01->{non_coding_features}};
    print "\n********For case $obj_65386_1*********\nAFTER _summarize_annotation there are $nc_ftr_count non_coding features.\n";

    # another genome object in workspace #65386
    $nc_ftr_count = @{$final_genome02->{non_coding_features}};
    print "\n********For case $obj_65386_2*********\nBEFORE _summarize_annotation there are $nc_ftr_count non_coding features.\n";

    lives_ok {
        ($final_genome02, $rast_ref) = $annoutil->_summarize_annotation(
              \%rast_details02, $final_genome02, $inputgenome02);
    } "_summarize_annotation runs successfully on genome $obj_65386_2";
    %rast_details02 = %{$rast_ref}; # dereference
    $nc_ftr_count = @{$final_genome02->{non_coding_features}};
    print "\n********For case $obj_65386_2*********\nAFTER _summarize_annotation there are $nc_ftr_count non_coding features.\n";

    # a genome object
    lives_ok {
        ($final_genome1, $rast_ref) = $annoutil->_summarize_annotation(
              \%rast_details1, $final_genome1, $inputgenome1);
    } "_summarize_annotation runs successfully on genome $obj_Ecoli";
    %rast_details1 = %{$rast_ref}; # dereference

    # an assembly object
    lives_ok {
        ($final_genome2, $rast_ref) = $annoutil->_summarize_annotation(
              \%rast_details2, $final_genome2, $inputgenome2);
    } "_summarize_annotation runs successfully on assembly $obj_asmb";
    %rast_details2 = %{$rast_ref}; # dereference
};


# Test _save_annotation_results with genome/assembly object refs in prod
subtest '_save_annotation_results' => sub {
    my ($aa_out, $out_msg);

    # a genome object in workspace #65386
    my $nc_ftr_count = @{$final_genome01->{non_coding_features}};
    print "\n********For case $obj_65386_1*********\nBEFORE _save_annotation_results there are $nc_ftr_count non_coding features.\n";

    my $ncoding_features01 = $final_genome01->{non_coding_features};
    my $cnt = 0;
    for my $ncoding_ftr (@{$ncoding_features01}) {
        if(exists($ncoding_ftr->{type})) {
            $cnt++;
            if ($cnt < 20) {
                print "type value: $ncoding_ftr->{type}\n";
            }
        }
    }
    ok ($nc_ftr_count==$cnt, "All $cnt non-coding features have defined type.\n");

    throws_ok {
        ($aa_out, $out_msg) = $annoutil->_save_annotation_results(
              $final_genome01, $rast_details01{parameters});
    } qr/Can't call method/,
      "_save_annotation_results throws an error on genome $obj_65386_1";

    $nc_ftr_count = @{$final_genome01->{non_coding_features}};
    print "\n********For case $obj_65386_1*********\nAFTER _save_annotation_results there are $nc_ftr_count non_coding features.\n";

    # another genome object in workspace #65386
    throws_ok {
        ($aa_out, $out_msg) = $annoutil->_save_annotation_results(
              $final_genome02, $rast_details02{parameters});
    } qr/Can't call method/,
      "_save_annotation_results throws an error on genome $obj_65386_2";

    # a genome object
    throws_ok {
        ($aa_out, $out_msg) = $annoutil->_save_annotation_results(
              $final_genome1, $rast_details1{parameters});
    } qr/Can't call method/,
      "_save_annotation_results throws an error on genome $obj_Ecoli";

    # an assembly object
    throws_ok {
        ($aa_out, $out_msg) = $annoutil->_save_annotation_results(
              \%rast_details2, $final_genome2, $inputgenome2);
    } qr/Can't call method/,
      "_save_annotation_results throws an error on assembly $obj_asmb";
};


#
## variables for testing _build_workflows, _run_rast_workflow
## and _pre_rast_call
my ($gc_wf_ret0, $gc_inputgenome0, $gc_wf0, %gc_details0,
    $gc_wf_ret1, $gc_inputgenome1, $gc_wf1, %gc_details1,
    $gc_wf_ret2, $gc_inputgenome2, $gc_wf2, %gc_details2);

# Test _build_workflows with genome/assembly object refs in prod
subtest '_build_workflows' => sub {
    # a genome object from https://narrative.kbase.us/narrative/ws.65386.obj.104
    my $params0 = {
         output_genome_name => 'build_gcwf_name0',
         output_workspace => $ws,
         object_ref => $obj_65386_1
    };
    lives_ok {
        ($gc_wf_ret0, $gc_inputgenome0) = $annoutil->_build_workflows($params0);
    } "_build_workflows returns normally on $obj_65386_1\n";
    %gc_details0 = %{ $gc_wf_ret0 };
    $gc_wf0 = $gc_details0{genecall_workflow};
    ok (exists($gc_details0{genecall_workflow}), "generate genecall workflow:\n".Dumper($gc_wf0));
    ok (@{$gc_inputgenome0->{features}} == 0, 'inputgenome has NO features.');
    ok (@{$gc_inputgenome0->{contigs}} > 0,
        "inputgenome has ".scalar @{$gc_inputgenome0->{contigs}}." contig(s).");

    my $params1 = {
         output_genome_name => 'build_gcwf_gn_name',
         output_workspace => $ws,
         object_ref => $obj_Ecoli
    };

    lives_ok {
        ($gc_wf_ret1, $gc_inputgenome1) = $annoutil->_build_workflows($params1);
    } "_build_workflows returns normally on $obj_Ecoli\n";
    %gc_details1 = %{ $gc_wf_ret1 };
    $gc_wf1 = $gc_details1{genecall_workflow};
    ok (exists($gc_details1{genecall_workflow}), "generate genecall workflow:\n".Dumper($gc_wf1));
    ok (@{$gc_inputgenome1->{features}} > 0,
        "inputgenome has ".scalar @{$gc_inputgenome1->{features}}." features.");
    ok (@{$gc_inputgenome1->{contigs}} > 0,
        "inputgenome has ".scalar @{$gc_inputgenome1->{contigs}}." contig(s).");

    my $params2 = {
         output_genome_name => 'build_gcwf_asmb_name',
         output_workspace => $ws,
         object_ref => $obj_asmb
    };
    lives_ok {
        ($gc_wf_ret2, $gc_inputgenome2) = $annoutil->_build_workflows($params2);
    } "_build_workflows returns normally on $obj_asmb\n";
    %gc_details2 = %{ $gc_wf_ret2 };
    $gc_wf2 = $gc_details2{genecall_workflow};
    ok (exists($gc_details2{genecall_workflow}), "generate genecall workflow:\n".Dumper($gc_wf2));
    ok (@{$gc_inputgenome2->{features}} == 0, 'inputgenome (assembly) has no features.');
    ok (@{$gc_inputgenome2->{contigs}} > 0,
        "inputgenome (assembly) has ".scalar @{$gc_inputgenome2->{contigs}}." contig(s).");
};



# Test _run_rast_workflow with genome/assembly object refs in prod
subtest '_run_rast_workflow' => sub {
    # a genome object from workspace #65386
    my $rast_gn0;
    lives_ok {
        $rast_gn0 = $annoutil->_run_rast_workflow($gc_inputgenome0, $gc_wf0);
    } '_run_rast_workflow call returns ERROR due to kmer data absence or other causes.';
    ok (@{$rast_gn0->{features}} == 0, "Returned genome has NO features.\n");
    cmp_deeply($rast_gn0, $gc_inputgenome0, 'rast workflow will not run locally');

    # a genome object
    my $rast_gn1;
    lives_ok {
        $rast_gn1 = $annoutil->_run_rast_workflow($gc_inputgenome1, $gc_wf1);
    } '_run_rast_workflow returns ERROR due to kmer data absence or other causes.';
    ok (@{$rast_gn1->{features}} > 0,
        "Returned genome has ".scalar @{$rast_gn1->{features}}." features.");
    cmp_deeply($rast_gn1, $gc_inputgenome1, 'rast workflow will not run locally');

    # an assembly object
    my $rast_gn2;
    lives_ok {
        $rast_gn2 = $annoutil->_run_rast_workflow($gc_inputgenome2, $gc_wf2);
    } '_run_rast_workflow call returns ERROR due to kmer data absence or other causes.';
    ok (@{$rast_gn2->{features}} == 0, "Returned genome has NO features.");
    cmp_deeply($rast_gn2, $gc_inputgenome2, 'rast workflow will not run locally');
};


# Test _run_rast_genecalls with genome/assembly object refs in prod
subtest '_run_rast_genecalls' => sub {
    # a genome object from https://narrative.kbase.us/narrative/ws.65386.obj.104
    my $params0 = {
         output_genome_name => 'build_gcwf_name0',
         output_workspace => $ws,
         object_ref => $obj_65386_1
    };
    my ($rast_ref0, $gc_gn0);
    lives_ok {
        ($gc_gn0, $rast_ref0) = $annoutil->_run_rast_genecalls($params0);
    } '_run_rast_genecalls returns ERROR due to kmer data absence or other causes.';
    ok (@{$gc_gn0->{features}} == 0, "Returned genome has NO features.");
    cmp_deeply($gc_gn0, $gc_inputgenome0, 'rast workflow will not run locally');

    # a genome object
    my $params1 = {
         output_genome_name => 'build_gcwf_gn_name',
         output_workspace => $ws,
         object_ref => $obj_Ecoli
    };
    my ($rast_ref1, $gc_gn1);
    lives_ok {
        ($gc_gn1, $rast_ref1) = $annoutil->_run_rast_genecalls($params1);
    } '_run_rast_genecalls returns ERROR due to kmer data absence or other causes.';
    ok (@{$gc_gn1->{features}} > 0, "Returned genome has features.");
    cmp_deeply($gc_gn1, $gc_inputgenome1, 'rast workflow will not run locally');

    # an assembly object
    my $params2 = {
         output_genome_name => 'build_gcwf_asmb_name',
         output_workspace => $ws,
         object_ref => $obj_asmb
    };
    my ($rast_ref2, $gc_gn2);
    lives_ok {
        ($gc_gn2, $rast_ref2) = $annoutil->_run_rast_genecalls($params2);
    } '_run_rast_genecalls returns ERROR due to kmer data absence or other causes.';
    ok (@{$gc_gn2->{features}} == 0, "Returned genome has NO features.");
    cmp_deeply($gc_gn2, $gc_inputgenome2, 'rast workflow will not run locally');
};


subtest 'Impl_annotate_genome' => sub {
    my $obj_asmb1 = '1234/56/7';
    my $assembly_obj_name = "Acidilobus_sp._CIS.fna";
    my $assembly_ref = prepare_assembly($assembly_obj_name);
    my $genome_obj_name = 'Acidilobus_sp_CIS';

    my $parms={
        "input_contigset" => $assembly_obj_name,
        "workspace" => $ws,
        "output_genome" => 'Acidilobus_sp_7',
        "scientific_name" => 'Acidilobus sp 7',
        "domain" => 'A',
        "genetic_code" => '4',
        "call_features_CDS_prodigal" => '1',
    };
    my $rast_ann;
    throws_ok {
        $rast_ann = $rast_impl->annotate_genome($parms);
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
    } qr/Error invoking method call_/,
      "test Impl annotate_genome on an assembly died.";
};

subtest '_validate_KB_objref' => sub {
	my $object_ref = 'qzhang:narrative_1581052755332/short_one_metagenome';
	my $passed_test = $annoutil->_validate_KB_objref($object_ref);
	ok ($passed_test, "$object_ref is a valid workspace object.\n");

	$object_ref = 'qzhang:narrative_1581052755332/short_one_metagenome/1';
	$passed_test = $annoutil->_validate_KB_objref($object_ref);
	ok ($passed_test, "$object_ref is a valid workspace object.\n");

	$object_ref = '52755332/short_one/1';
	$passed_test = $annoutil->_validate_KB_objref($object_ref);
	ok ($passed_test, "$object_ref is a valid workspace object.\n");

	$object_ref = '5332/345/3';
	$passed_test = $annoutil->_validate_KB_objref($object_ref);
	ok ($passed_test, "$object_ref is a valid workspace object.\n");

	$object_ref = '5332/3wda9/123';
	$passed_test = $annoutil->_validate_KB_objref($object_ref);
	ok ($passed_test, "$object_ref is a valid workspace object.\n");

	$object_ref = '5332/39';
	$passed_test = $annoutil->_validate_KB_objref($object_ref);
	ok ($passed_test, "$object_ref is a valid workspace object.\n");

	$object_ref = '5332/3wda9';
	$passed_test = $annoutil->_validate_KB_objref($object_ref);
	ok ($passed_test, "$object_ref is a valid workspace object.\n");

	$object_ref = '5332/3wda9/a';
	$passed_test = $annoutil->_validate_KB_objref($object_ref);
	ok ($passed_test == 0, "$object_ref is an invalid workspace object.\n");
};


subtest '_check_annotation_params' => sub {
    my $obj = '1234/56/7';

    my $missing_params = "Missing required parameters for annotating genome.\n";
    my $req1 = "'output_workspace' is required for running rast_genome.\n";
    my $req2 = "'object_ref' is required for running rast_genome.\n";

    throws_ok {
        $annoutil->_check_annotation_params()
    } qr/$missing_params/,
        '_check_annotation_params dies without params';

    throws_ok {
        $annoutil->_check_annotation_params( {} )
    } qr/$missing_params/,
        '_check_annotation_params dies with an empty hashref';

    throws_ok {
        my $p = {output_workspace => $ws,
                 output_genome_name => $outgn_name};
        print "input parameter=\n". Dumper($p);
        $annoutil->_check_annotation_params($p)
    } qr/$req2/,
        '_check_annotation_params dies with no object_ref';

    throws_ok {
        my $p = {output_workspace => $ws,
                 output_genome_name => $outgn_name,
                 object_ref => ''};
        print "input parameter=\n". Dumper($p);
        $annoutil->_check_annotation_params($p)
    } qr/$req2/,
        '_check_annotation_params dies with blank object_ref';

    throws_ok {
        $annoutil->_check_annotation_params(
            {object_ref => $obj,
             output_genome_name => $outgn_name})
    } qr/$req1/,
        '_check_annotation_params_metag dies with no outpout_workspace';

    throws_ok {
        $annoutil->_check_annotation_params(
            {workspace => $ws,
             output_genome_name => $outgn_name,
             obect_ref => $obj})
    } qr/$req1/,
        '_check_annotation_params dies with wrong workspace key';

    throws_ok {
        $annoutil->_check_annotation_params(
            {output_workspace => '',
             output_genome_name => $outgn_name,
             obect_ref => $obj})
    } qr/$req1/,
        '_check_annotation_params dies with blank workspace name';

    lives_ok {
        $annoutil->_check_annotation_params(
            {output_workspace => $ws,
             output_genome_name => $outgn_name,
             object_ref => 'abc/1/2'});
    } '_check_annotation_params object_ref check ok';

    lives_ok {
        $annoutil->_check_annotation_params(
            {output_workspace => 'ab:c',
             output_genome_name => $outgn_name,
             object_ref => '456/1/2'});
    } '_check_annotation_params workspace name check ok';

    # _check_annotation_params passed
    my $expected = {
             output_workspace => $ws,
             output_genome_name => $outgn_name,
             object_ref => '456/1/2'};
    my $set_default_ok = '_check_annotation_params sets the default value for output_genome_name.';

    my $ret = $annoutil->_check_annotation_params(
            {output_workspace => $ws,
             output_genome_name => undef,
             object_ref => '456/1/2'});
    ok ($ret->{output_genome_name} eq $expected->{output_genome_name},
        'When undefined, '.$set_default_ok);

    $ret = $annoutil->_check_annotation_params(
            {output_workspace => $ws,
             output_genome_name => '',
             object_ref => '456/1/2'});
    ok ($ret->{output_genome_name} eq $expected->{output_genome_name},
        'When blank, '.$set_default_ok);
};


subtest '_check_bulk_annotation_params' => sub {
    my $error_message = qr/ERROR:Missing required inputs/;
    my $error_mand = qr/Mandatory arguments missing/;

    my $params = {
        "output_GenomeSet_name" => "out_genomeSet"
    };
    throws_ok {
        my $ret_parms1 = $annoutil->_check_bulk_annotation_params($params);
    } qr/'output_workspace' is required/,
      'Missing required parameter output_workspace die correctly'
      or diag explain $params;

    $params = {
        "output_workspace" => get_ws_name()
    };
    my $expected = {
       'input_assemblies' => [],
       'input_genomes' => [],
       'input_text' => '',
       'output_GenomeSet_name' => 'rasted_GenomeSet_name'
    };

    my $set_default_ok = '_check_annotation_params sets the default value for output_GenomeSet_name.';
    $params->{input_text} = '';
    my $ret_parms2 = $annoutil->_check_bulk_annotation_params($params);
    ok ($ret_parms2->{output_GenomeSet_name} eq $expected->{output_GenomeSet_name},
        "When undefined, ".$set_default_ok);
    ok ($ret_parms2->{output_workspace} eq $params->{output_workspace},
        "output_workspace is defined");
    ok (!@{$ret_parms2->{input_genomes}} && !@{$ret_parms2->{input_assemblies}},
        "No input genome or assembly was specified as input.");
    print Dumper($ret_parms2);

    $params->{input_genomes} = [];
    $params->{input_assemblies} = [];
    $params->{input_text} = '';
    my $ret_parms3 = $annoutil->_check_bulk_annotation_params($params);
    ok (!@{$ret_parms3->{input_genomes}} && !@{$ret_parms3->{input_assemblies}},
        "No input genome or assembly was specified as input.");

    $params->{input_genomes} = ["48109/9/1"]; # array of a prod object
    $params->{input_text} = '';
    my $ret_parms4 = $annoutil->_check_bulk_annotation_params($params);
    ok ($ret_parms4->{input_genomes} eq $params->{input_genomes},
        "Input genome array is not empty.");

    $params->{input_genomes} = []; # array of a prod object
    $params->{input_text} = '48109/9/1;123/4/5';
    my $ret_parms5 = $annoutil->_check_bulk_annotation_params($params);
    ok ($ret_parms5->{input_text} eq $params->{input_text},
        "Input text is not empty.");
    ok (!@{$ret_parms5->{input_genomes}} && !@{$ret_parms5->{input_assemblies}},
        "No input genome or assembly was specified as input.");

    $params->{input_genomes} = "48109/9/1"; # non-array
    $params->{input_text} = '';
    my $ret_parms6 = $annoutil->_check_bulk_annotation_params($params);
    ok ($ret_parms6->{input_genomes}->[0] eq "48109/9/1",
        "Non array input genome converted into array");

    $params->{input_genomes} = [$obj_Ecoli]; # array of prod objects
    $params->{input_assemblies} = [$obj_asmb]; # array of prod objects
    $params->{input_text} = '';
    my $ret_parms7 = $annoutil->_check_bulk_annotation_params($params);
    ok (@{$ret_parms7->{input_genomes}} && @{$ret_parms7->{input_assemblies}},
        "Both input_genomes and input_assemblies arrays are not empty.");

    $params->{input_genomes} = [$obj_Echinacea, $obj_Ecoli]; # array of prod objects
    $params->{input_assemblies} = [];
    $params->{input_text} = '';
    my $ret_parms8 = $annoutil->_check_bulk_annotation_params($params);
    ok (@{$ret_parms8->{input_genomes}}==2 && @{$ret_parms8->{input_assemblies}}==0,
        "The input_genomes array has 2 elements while input_assemblies is empty.");

    $params->{input_assemblies} = [$obj_asmb_refseq, $obj_asmb]; # array of prod objects
    $params->{input_genomes} = [];
    $params->{input_text} = '';
    my $ret_parms9 = $annoutil->_check_bulk_annotation_params($params);
    ok (@{$ret_parms9->{input_genomes}}==0 && @{$ret_parms9->{input_assemblies}}==2,
        "The input_assemblies array has 2 elements while input_genomes is empty.");

    $params->{input_genomes} = [$obj_Echinacea, $obj_Ecoli]; # array of prod objects
    my $ret_parms10 = $annoutil->_check_bulk_annotation_params($params);
    ok (@{$ret_parms10->{input_genomes}}==2 && @{$ret_parms10->{input_assemblies}}==2,
        "Both input_genomes and input_assemblies arrays have 2 elements.");
};


subtest '_parseNwrite_gff' => sub {
    # testing get the gff from a genome using obj ids from prod ONLY
    my $gff_fpath;
    my $attr_delimiter = '=';
    lives_ok {
        $gff_fpath = $annoutil->_write_gff_from_genome($obj_Echinacea);
    } 'Writing gff from a genome runs ok';

    ok((-e $gff_fpath), "GFF file created for $obj_Echinacea.\n");
    ok ((-s $gff_fpath), "GFF file written for $obj_Echinacea to $gff_fpath.\n");

    print "First 20 lines of the GFF file written from genome $obj_Echinacea:\n";
    $annoutil->_print_fasta_gff(0, 20, $gff_fpath);

    my $temp_gff_file = catfile($tmp_write_dir, 'test_gff_file.gff');
    copy $gff_fpath, $temp_gff_file;
    print "First 20 lines of the GFF file copied from genome $gff_fpath:\n";
    $annoutil->_print_fasta_gff(0, 20, $temp_gff_file);

    my $gff_contents1;
    lives_ok {
        ($gff_contents1, $attr_delimiter) = $annoutil->_parse_gff($gff_fpath, '=');
    } "Testing _parse_gff on $gff_fpath succeeded.";
    ok( @{$gff_contents1} >0, "Parsing GFF on $gff_fpath returns result.\n");

    my $gff_contents2;
    lives_ok {
        ($gff_contents2, $attr_delimiter) = $annoutil->_parse_gff($temp_gff_file, '=');
    } "Testing _parse_gff on $temp_gff_file succeeded.";
    ok( @{$gff_contents1} >0, "Parsing GFF on $temp_gff_file returns result.\n");

    is_deeply($gff_contents1, $gff_contents2, 'GFF data structures should be the same!');

    my $Echinacea_gff_file = "data/Echinacea_purpurea_1762.gff";
    my $gff_contents3;
    lives_ok {
        ($gff_contents3, $attr_delimiter) = $annoutil->_parse_gff($Echinacea_gff_file, '=');
    } "Testing _parse_gff on file $Echinacea_gff_file succeeded.";
    ok( @{$gff_contents3} >0, "Parsing GFF on $Echinacea_gff_file returns result.\n");
    print "Parsed ". scalar @{$gff_contents3}." GFF contents.\n";

    is_deeply($gff_contents1, $gff_contents3, 'GFF data structures should be the same 2!');

    my $test_gff_file_written;
    lives_ok {
        $test_gff_file_written = catfile($rast_genome_dir, 'test_written.gff');
        $annoutil->_write_gff($gff_contents1, $test_gff_file_written , '=');
    } "Writing the gff contents back to a gff file is ok";
    ok ( (-s $test_gff_file_written), "GFF file written to $test_gff_file_written.\n");
};


subtest '_write_gff_from_genome' => sub {
    # testing get the gff from a genome using obj ids from prod ONLY
    my $gff_fpath;
    lives_ok {
        $gff_fpath = $annoutil->_write_gff_from_genome($obj_Ecoli);
    } 'Writing gff from a genome runs ok';

    ok((-e $gff_fpath), "GFF file created for $obj_Ecoli.\n");
    ok ((-s $gff_fpath), "GFF file written for $obj_Ecoli.\n");

    print "First 10 lines of the GFF file:\n";
    $annoutil->_print_fasta_gff(0, 10, $gff_fpath);

    lives_ok {
        $gff_fpath = $annoutil->_write_gff_from_genome($obj_Echinacea);
    } 'Writing gff from a genome runs ok';

    ok((-e $gff_fpath), "GFF file created for $obj_Echinacea.\n");
    ok ((-s $gff_fpath), "GFF file written for $obj_Echinacea.\n");

    print "ALL lines of the GFF file:\n";
    $annoutil->_print_fasta_gff(0, 2000, $gff_fpath);
};


subtest 'anno_utils_rast_genome' => sub {
    # testing anno_utils rast_genome using obj ids from prod ONLY
    my @regexes = ( qr/ERROR calling rast run_pipeline/,
                    qr/Can\'t call method/);
    my $allre= "(".join("|",@regexes).")";
    my ($parms, $rast_ret, $genome_obj);
    $parms = {
        "object_ref" => $obj_Ecoli,
        "output_genome_name" => "rasted_ecoli_prod",
        "output_workspace" => $ws
    };

    throws_ok {
        $rast_ret = $annoutil->rast_genome($parms);
    } qr/$allre/,
        '$annoutil->rast_genome returns ERROR due to kmer data absence or other causes.';
    $parms = {
        "object_ref" => $obj_Echinacea,
        "output_genome_name" => "rasted_Echinace_prod",
        "output_workspace" => $ws
    };

    throws_ok {
        $rast_ret = $annoutil->rast_genome($parms);
    } qr/$allre/,
        '$annoutil->rast_genome returns ERROR due to kmer data absence or other causes.';
    $parms = {
        "object_ref" => $obj_asmb,
        "output_genome_name" => "rasted_assembly",
        "output_workspace" => $ws
    };

    lives_ok {
        $rast_ret = $annoutil->rast_genome($parms);
    } '$annoutil->rast_genome returns the original input because of empty features.';
    if(defined($rast_ret) && defined($rast_ret->{output_genome_ref})) {
        ok (($rast_ret->{output_genome_ref} =~ m/[^\\w\\|._-]/), "rast_genome returns a VALID ref: $rast_ret->{output_genome_ref}");
    }
};


## testing Impl_rast_genome_assembly using obj ids from prod ONLY
subtest 'Impl_rast_genome_assembly' => sub {
    my $parms = {
        "object_ref" => $obj_Ecoli,
        "output_genome_name" => "rasted_genome",
        "output_workspace" => $ws,
        "create_report" => 1
    };
    my $rast_ret;
    lives_ok {
        $rast_ret = $rast_impl->rast_genome_assembly($parms);
    } 'Impl rast_genome call returns normally on genome.';
    ok ($rast_ret->{output_genome_ref} =~ m/[^\\w\\|._-]/,
        "rast_genome_assembly returns a VALID ref: $rast_ret->{output_genome_ref}");

    $parms = {
        "object_ref" => $obj_asmb,
        "output_genome_name" => "rasted_assembly",
        "output_workspace" => $ws
    };
    lives_ok {
        $rast_ret = $rast_impl->rast_genome_assembly($parms);
    } 'Impl rast_genome call returns without annotation due to local assembly did not run.';
    ok (!defined($rast_ret ->{output_genome_ref}),
        "due to local annotation on assembly with empty features, no rast was run.");
};


## testing _get_bulk_rast_parameters using obj ids from prod ONLY
subtest '_get_bulk_rast_parameters' => sub {
    my $parms = {
          'input_genomes' => [
                               '55141/242/1',
                               '55141/212/1'
                             ],
          'create_report' => 0,
          'output_workspace' => 'test_RAST_SDK_1592592098000',
          'output_GenomeSet_name' => 'out_genomeSet',
          'input_assemblies' => [
                                  '55141/266/3',
                                  '55141/243/1'
                                ],
          'input_text' => '55141/266/3;55141/212/1;63171/394/1'
    };
    my @refs = split(';', $parms->{input_text});
    push(@refs, @{$parms->{input_assemblies}});
    push(@refs, @{$parms->{input_genomes}});

    my $params;
    lives_ok {
        $params = $annoutil->_get_bulk_rast_parameters($parms);
    } "anno_utils _get_bulk_rast_parameters call returns normally.";

    for my $ref (@refs) {
        ok ($annoutil->_value_in_array($ref, $params), "$ref is included in the parameters\n");
    }

    # testing invalid object name/id inputs
    $parms->{input_text} = '55141/266/3;55141/212/1;63171/394/1;abc/3/2;333/22/a';
    @refs = split(';', $parms->{input_text});

    lives_ok {
        $params = $annoutil->_get_bulk_rast_parameters($parms);
    } "anno_utils _get_bulk_rast_parameters call returns normally.";
    for my $ref (@refs) {
        if ($ref ne 'abc/3/2' && $ref ne '333/22/a') {
            ok ($annoutil->_value_in_array($ref, $params), "$ref is included in the parameters\n");
		}
        else {
            ok ($annoutil->_value_in_array($ref, $params) == 0, "invalid reference $ref is excluded from the parameters\n");
		}
    }
};


## testing Impl_rast_genomes_assemblies using obj ids from prod ONLY
subtest 'rast_genomes_assemblies' => sub {
    my $params = {
        "output_GenomeSet_name" => "out_genomeSet"
    };

    lives_ok {
        $params->{output_workspace} = get_ws_name();
        $params->{input_genomes} = [$obj_Ecoli]; # array of prod objects
        $params->{input_assemblies} = [$obj_asmb]; # array of prod objects
        $params->{input_text} = '';
        my $ret_ann6 = $rast_impl->rast_genomes_assemblies($params);
    } "rast_impl rast_genomes_assemblies call on 1N1 returns normally.";

    lives_ok {
        $params->{output_workspace} = get_ws_name();
        #$params->{input_genomes} = ["31020/5/1"]; # array of an appdev object
        $params->{input_genomes} = [$obj_Echinacea, $obj_Ecoli]; # array of prod objects
        $params->{input_assemblies} = [];
        $params->{input_text} = '';
        my $ret_ann7 = $rast_impl->rast_genomes_assemblies($params);
    } "rast_impl rast_genomes_assemblies call on array of 2 genomes returns normally.";

    lives_ok {
        $params->{output_workspace} = get_ws_name();
        $params->{input_assemblies} = [$obj_asmb_refseq, $obj_asmb]; # array of prod objects
        $params->{input_genomes} = [];
        $params->{input_text} = '';
        my $ret_ann8 = $rast_impl->rast_genomes_assemblies($params);
    } "rast_impl rast_genomes_assemblies call on array of 2 assemblies returns normally.";

    lives_ok {
        $params->{output_workspace} = get_ws_name();
        $params->{input_genomes} = [$obj_Echinacea, $obj_Ecoli]; # array of prod objects
        $params->{input_assemblies} = [$obj_asmb_refseq, $obj_asmb]; # array of prod objects
        #$params->{input_text} = '55141/266/3;55141/212/1;63171/394/1';
        $params->{input_text} = '36230/12/5;36230/13/5; 36230/14/5';
        my $ret_ann9 = $rast_impl->rast_genomes_assemblies($params);
    } "rast_impl rast_genomes_assemblies call on two arrays returns normally.";
};


# Test checking annotate_genomes input params for empty input_genomes and blank/undef genome_text
subtest 'annotation_genomes_throw_messages' => sub {
    my $error_message = qr/ERROR:Missing required inputs/;

    my $params = {
        "output_genome" => "out_genome_name",
        "workspace" => get_ws_name()
    };
    throws_ok {
        $params->{genome_text} = '';
        my $ret_ann1 = $rast_impl->annotate_genomes($params);
    } $error_message,
      'Blank genome_text plus undef input_genoms die correctly'
      or diag explain $params;

    $params = {
        "output_genome" => "out_genome_name",
        "workspace" => get_ws_name()
    };
    throws_ok {
        $params->{input_genomes} = [];
        my $ret_ann2 = $rast_impl->annotate_genomes($params);
    } $error_message,
      'Empty input_genomes plus undef genome_text die correctly'
      or diag explain $params;

    $params = {
        "output_genome" => "out_genome_name",
        "workspace" => get_ws_name()
    };
    throws_ok {
        $params->{input_genomes} = [];
        $params->{genome_text} = '';
        my $ret_ann3 = $rast_impl->annotate_genomes($params);
    } $error_message,
      'Blank genome_text AND empty input_genoms die correctly'
      or diag explain $params;
};

#----- For checking the stats of a given obj id in prod ONLY-----#
my $stats_ok = 'stats generation runs ok.\n';

subtest '_generate_stats_from_aa & from_gffContents' => sub {
    my ($gff_contents, $attr_delimiter) = ([], '=');
    my $gff_path = '';

    my $aa_stats_ok = 'stats generation from AA runs ok.\n';
    my $gff_stats_ok = 'stats generation from GFF runs ok.\n';

    my %ret_stats = ();

    # $obj_Ecoli
    lives_ok {
        %ret_stats = $annoutil->_generate_stats_from_aa($obj_Ecoli);
    } $aa_stats_ok;
    print "Stats from annoated genome on $obj_Ecoli:\n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generated from annotated genome $obj_Ecoli.");

    # $obj_Ecoli_ann
    %ret_stats = ();
    lives_ok {
        %ret_stats = $annoutil->_generate_stats_from_aa($obj_Ecoli_ann);
    } $aa_stats_ok;
    print "Stats from annoated genome on $obj_Ecoli_ann:\n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generated from annotated genome $obj_Ecoli_ann.");

    %ret_stats = ();
    lives_ok {
        $gff_path = $annoutil->_write_gff_from_genome($obj_Ecoli_ann);
        ($gff_contents, $attr_delimiter) = $annoutil->_parse_gff($gff_path, $attr_delimiter);
        %ret_stats = $annoutil->_generate_stats_from_gffContents($gff_contents);
    } $gff_stats_ok;
    #print "Stats from GFF on $obj_Ecoli_ann: \n".Dumper(\%ret_stats);
    is(keys %ret_stats, 2, "Statistics generated from gffContents on $obj_Ecoli_ann.");

    # $obj_asmb_ann
    %ret_stats = ();
    lives_ok {
        %ret_stats = $annoutil->_generate_stats_from_aa($obj_asmb_ann);
    } $aa_stats_ok;
    print "Stats from annoated assembly on $obj_asmb_ann:\n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generated from annotated assembly $obj_asmb_ann.");

    %ret_stats = ();
    lives_ok {
        $gff_path = $annoutil->_write_gff_from_genome($obj_asmb_ann);
        ($gff_contents, $attr_delimiter) = $annoutil->_parse_gff($gff_path, $attr_delimiter);
        %ret_stats = $annoutil->_generate_stats_from_gffContents($gff_contents);
    } $gff_stats_ok;
    #print "Stats from GFF on $obj_asmb_ann: \n".Dumper(\%ret_stats);
    is(keys %ret_stats, 2, "Statistics generated from gffContents on $obj_asmb_ann.");
};

# testing generate_genome_report using obj ids from prod ONLY
subtest 'generate_genome_report1' => sub {
    my $stats_ok = 'stats generation runs ok.\n';

    my $gff_path1 = $annoutil->_write_gff_from_genome($obj_Ecoli);
    my $gff_path2 = $annoutil->_write_gff_from_genome($obj_Ecoli_ann);

    my ($gff_contents1, $attr_delimiter) = ([], '=');
    my %ret_stats1;
    lives_ok {
        ($gff_contents1, $attr_delimiter) = $annoutil->_parse_gff($gff_path1, $attr_delimiter);
        %ret_stats1 = $annoutil->_generate_stats_from_gffContents($gff_contents1);
    } $stats_ok;
    is(keys %ret_stats1, 2, "_generate_stats_from_gffContents on $obj_Ecoli should return non-empty.\n");

    my $gff_contents2 = [];
    my %ret_stats2;
    lives_ok {
        ($gff_contents2, $attr_delimiter) = $annoutil->_parse_gff($gff_path2, $attr_delimiter);
        %ret_stats2 = $annoutil->_generate_stats_from_gffContents($gff_contents2);
    } $stats_ok;
    is(keys %ret_stats2, 2, "_generate_stats_from_gffContents on $obj_Ecoli_ann should return non-empty.\n");
    ok(exists($ret_stats2{gene_role_map}), '_generate_stats_from_gffContents stats contains gene_roles.');
    ok(exists($ret_stats2{function_roles}), '_generate_stats_from_gffContents stats contains function roles.');

    my %ftr_tab = $annoutil->_get_feature_function_lookup($test_ftrs);
    #print "\nFeature lookup:\n".Dumper(\%ftr_tab);
    my $test_msg = "Test message for reporting of annotation details";
    my $ftr_cnt = scalar @{$test_ftrs};
    my $ret_rpt = $annoutil->_generate_genome_report(
                $obj_Ecoli_ann, $gff_contents2, \%ftr_tab, $ftr_cnt, $test_msg);
    ok( exists($ret_rpt->{report_ref}), 'Report generation returns report_ref.');
    ok( exists($ret_rpt->{report_name}), 'Report generation returns report_name.');
    ok( exists($ret_rpt->{output_genome_ref}), 'Report generation returns output_gemome_ref.');
};


## testing generate_genome_report using obj ids from prod ONLY
subtest 'generate_genome_report2' => sub {
    my $stats_ok = 'stats generation runs ok.\n';

    my $gff_path = $annoutil->_write_gff_from_genome($obj_asmb_ann);

    my ($gff_contents, $attr_delimiter) = ([], '=');

    my %ret_stats;
    lives_ok {
        ($gff_contents, $attr_delimiter) = $annoutil->_parse_gff($gff_path, $attr_delimiter);
        %ret_stats = $annoutil->_generate_stats_from_gffContents($gff_contents);
    } $stats_ok;
    is(keys %ret_stats, 2, "_generate_stats_from_gffContents on $obj_asmb_ann should return non-empty.\n");
    ok(exists($ret_stats{gene_role_map}), '_generate_stats_from_gffContents stats contains gene_roles.');
    ok(exists($ret_stats{function_roles}), '_generate_stats_from_gffContents stats contains function roles.');

    my %ftr_tab = $annoutil->_get_feature_function_lookup($test_ftrs);
    my $test_msg = "Test message for reporting of annotation details";
    my $ftr_cnt = scalar @{$test_ftrs};
    my $ret_rpt = $annoutil->_generate_genome_report(
                $obj_asmb_ann, $gff_contents, \%ftr_tab, $ftr_cnt, $test_msg);
    ok( exists($ret_rpt->{report_ref}), 'Report generation returns report_ref.');
    ok( exists($ret_rpt->{report_name}), 'Report generation returns report_name.');
    ok( exists($ret_rpt->{output_genome_ref}), 'Report generation returns output_gemome_ref.');
};

# testing _generate_stats_from_aa using obj ids from prod ONLY
subtest '_generate_stats_from_aa' => sub {
    my %ret_stats = $annoutil->_generate_stats_from_aa($obj_Ecoli_ann);
    ok(keys %ret_stats, "Statistics generation from $obj_Ecoli_ann returns result.");

    %ret_stats = $annoutil->_generate_stats_from_aa($obj_Echinacea_ann);
    ok(keys %ret_stats, "Statistics generation from $obj_Echinacea_ann returns result.");
};

# Testing two small uniq functions
subtest 'annoutil_uniq_functions' => sub {
    my @words_array = qw(foo bar baz foo zorg baz);
    my @expected_array = qw(foo bar baz zorg);

    # uniq arrays
    my @data = @words_array;
    my @uniq_arr = $annoutil->_uniq(@data);
    cmp_deeply(sort @expected_array, sort @uniq_arr, 'unique array is correct');

    # uniq refs to arrays
    my $data = \@words_array;
    my $uniq_ref = $annoutil->_uniq_ref($data);
    my @sorted_ret = sort @{$uniq_ref};
    cmp_deeply(sort @expected_array, @sorted_ret, 'unique array ref is correct');
};

done_testing();


my $err = undef;
if ($@) {
    $err = $@;
}
eval {
    if (defined($ws)) {
        $ws_client->delete_workspace({workspace => $ws});
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
