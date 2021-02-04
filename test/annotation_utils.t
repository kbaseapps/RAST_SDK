use Test::Most;
use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catfile splitpath);
use File::Copy;
use Carp qw(croak);
use File::Compare;
use Config::Simple;
use Bio::KBase::AuthToken;
use Clone 'clone';

use RASTTestUtils;

use installed_clients::WorkspaceClient;
use installed_clients::GenomeFileUtilClient;
use RAST_SDK::RAST_SDKImpl;

use_ok "RAST_SDK::AnnotationUtils";

## global variables
my $token       = $ENV{ 'KB_AUTH_TOKEN' };
my $config_file = $ENV{ 'KB_DEPLOYMENT_CONFIG' };
my $config      = Config::Simple->new( $config_file )->get_block( 'RAST_SDK' );
my $auth_token  = Bio::KBase::AuthToken->new(
    token         => $token,
    ignore_authrc => 1,
    auth_svc      => $config->{ 'auth-service-url' }
);

my $ws_client     = RASTTestUtils::get_ws_client();
my $ws_name       = RASTTestUtils::get_ws_name();

my $call_back_url = $ENV{ SDK_CALLBACK_URL };
my $ctx           = LocalCallContext->new( $token, $auth_token->user_id );
$RAST_SDK::RAST_SDKServer::CallContext = $ctx;

my $gbff_file = catfile('/kb/module/test', 'data/Clostridium_botulinum.gbff');
my $fa_LOng = catfile('/kb/module/test', 'data/LOng_contig_names.fa');
my $Echinacea_gff_file = catfile('/kb/module/test', 'data/Echinacea_purpurea_1762.gff');

my $out_name = 'annotated_metag';
my $outgn_name = 'rast_annotated_genome';

my $rast_impl = RAST_SDK::RAST_SDKImpl->new();
my $annoutil  = RAST_SDK::AnnotationUtils->new( $config, $ctx );

my $scratch = $config->{ 'scratch' };    #'/kb/module/work/tmp';
my $rast_genome_dir = $annoutil->_create_rast_subdir( $scratch, "genome_annotation_dir_" );

sub get_feature_locations {
    my ($ctgID_hash, $ftr_arr) = @_;
    my $ret_locs  = [];
    for my $ftr (@{$ftr_arr}) {
        for my $loc (@{$ftr->{location}} ) {
            if( $ctgID_hash->{$loc->[0]} ) {
                $loc->[0] = $ctgID_hash->{$loc->[0]};
            }
        }
        push @{$ret_locs}, $ftr->{location};
    }
    return $ret_locs;
}

##-----------------Test Blocks--------------------##

my $obj_Echinacea = "55141/242/1";  # prod genome
my $obj_Echinacea_ann = "55141/247/1";  # prod genome
my $obj_Ecoli = "55141/212/1";  # prod genome
my $obj_Ecoli_ann = "55141/252/1";  # prod genome
my $obj_asmb = "55141/243/1";  # prod assembly
my $obj_asmb_ann = "55141/244/1";  # prod assembly
my $obj_asmb_refseq = "55141/266/3";  # prod assembly
my $obj1 = "37798/14/1";  # appdev
my $obj2 = "37798/15/1";  # appdev
my $obj3 = "55141/77/1";  # prod KBaseGenomeAnnotations.Assembly

my $asmb_fasta;

lives_ok {
    $asmb_fasta = $annoutil->_get_fasta_from_assembly($obj_asmb);
} 'Got FASTA from saved assembly file';

unless ( $asmb_fasta ) {
    fail 'Could not set up test data; skipping all tests';
    RASTTestUtils::clean_up();
    done_testing;
    exit 1;
}

#####RAST_SDK Module test objects #####
my $obj2_1 = "63171/315/1";
my $obj_65386_1 = '65386/2/1';  # same as 63171/436/1, i.e., GCF_003058445.1/Clostridium_botulinum
my $obj_65386_2 = '65386/12/1'; # Methanosarcina_acetivorans_C2A
my $obj_refseq_GCF = '63171/483/1';  # refseq_GCF_GCF_900128725.1
my $GEBA_1003_asmb = '63171/564/1';
my $GEBA_1003_asmb_ann = '63171/569/1';

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

=begin
# test the _map_location_contigIDs function
subtest '_map_location_contigIDs' => sub {
    my $arr_with_locations = [
        {location => [[
            'contigID_1',
            123,
            '+',
            456
        ]]},
        {location => [[
            'contigID_2',
            123,
            '+',
            456
        ]]},
        {location => [[
            'NZ_CP028858.1',
            777,
            '+',
            898
        ]]},
        {location => [[
            'contigID_4',
            1000,
            '+',
            1234
        ]]},
        {location => [
          [
            'contigID_3',
            259992,
            '-',
            36
          ],
          [
            'contigID_3',
            259915,
            '-',
            36
          ]
        ]}
    ];
    my $contigID_hash = {
          'contigID_1' => 'NZ_CP028859.1',
          'contigID_2' => 'NZ_CP028860.1',
          'contigID_3' => 'NC_003552.1',
          'contigID_4' => 'CP0035411Candidatus_Carsonella_ruddii_CE_isolate_Thao2000_complete_genome'
    };
    my $arr_mapped = $annoutil->_map_location_contigIDs($contigID_hash, $arr_with_locations);
    my $expected = [
        {location => [[
            'NZ_CP028859.1',
            123,
            '+',
            456
        ]]},
        {location => [[
            'NZ_CP028860.1',
            123,
            '+',
            456
        ]]},
        {location => [[
            'NZ_CP028858.1',
            777,
            '+',
            898
        ]]},
        {location => [[
            'CP0035411Candidatus_Carsonella_ruddii_CE_isolate_Thao2000_complete_genome',
            1000,
            '+',
            1234
        ]]},
        {location => [
          [
            'NC_003552.1',
            259992,
            '-',
            36
          ],
          [
            'NC_003552.1',
            259915,
            '-',
            36
          ]
        ]}
    ];

    cmp_deeply $arr_mapped, $expected, 'remapping location contig_ids correctly';
};


## Re-mapping the contigIDs back to their original (long) names
subtest '_remap_contigIDs' => sub {
    my $contigID_hash = {};
    my $gn0 = {
            contig_ids => ['fake_contig_id_1', 'fake_contig_id_2'],
    };
    my $gn1 = $annoutil->_remap_contigIDs($contigID_hash, $gn0);
    cmp_deeply($gn1, $gn0, "No contigID mapping happened--contigID_hash was empty");

    $contigID_hash = {
          'contigID_1' => 'NZ_CP028859.1',
          'contigID_2' => 'NZ_CP028860.1',
          'contigID_3' => 'NZ_CP028859.1',
          'contigID_4' => 'CP0035411Candidatus_Carsonella_ruddii_CE_isolate_Thao2000_complete_genome'
    };
    $gn0 = { contig_ids => [] };
    $gn1 = $annoutil->_remap_contigIDs($contigID_hash, $gn0);
    cmp_deeply $gn1->{ contig_ids }, [], 'No contigID remapping: empty contig ID array';

    my $gn = {
           contig_ids => ['contigID_1', 'contigID_2', 'contigID_3', 'contigID_4']
    };
    $gn = $annoutil->_remap_contigIDs($contigID_hash, $gn);
    cmp_deeply
    $gn->{ contig_ids },
    [  $contigID_hash->{ contigID_1 }, $contigID_hash->{ contigID_2 }, $contigID_hash->{ contigID_3 }, $contigID_hash->{ contigID_4 } ],
    'contig IDs remapped correctly';
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
         output_workspace => $ws_name,
         object_ref => $obj_Ecoli
    };

    my $obj_ref = $parameters->{object_ref};
    my $ret_gn;
    lives_ok {
        $ret_gn = $annoutil->_get_genome($obj_ref);
    } '_get_genome runs successfully';
    ok (@{$ret_gn->{features}} > 0, 'Genome has features!');
    is ($ret_gn->{assembly_ref}, '2901/78/1', 'found genome assembly ref');
};

subtest '_get_contigs' => sub {
    my $parameters = {
         output_genome_name => 'test_out_gn_name',
         output_workspace => $ws_name,
         object_ref => $obj_Ecoli
    };

    my $obj_ref = $parameters->{object_ref};
    my $obj = $annoutil->_fetch_object_data($obj_ref);

    my ($ret_contig_obj, $cid_hash);
    lives_ok {
        ($ret_contig_obj, $cid_hash) = $annoutil->_get_contigs($obj_ref);
    } "_get_contigs returns empty contigs on genome object $obj_ref";
    ok (@{$ret_contig_obj->{contigs}} eq 0, 'returned object has no contigs!');
    ok (keys %$cid_hash == 0, 'no contigID_hash!');

    lives_ok {
        ($ret_contig_obj, $cid_hash) = $annoutil->_get_contigs($obj->{assembly_ref});
    } "_get_contigs runs successfully on assembly_ref $obj->{assembly_ref}";
    ok (@{$ret_contig_obj->{contigs}} > 0, 'returned object has contigs!');
    ok (keys %$cid_hash, 'There is a contigID_hash:'.Dumper($cid_hash));

    lives_ok {
        ($ret_contig_obj, $cid_hash) = $annoutil->_get_contigs($obj_asmb);
    } "_get_contigs runs successfully on assembly $obj_asmb";
    ok (@{$ret_contig_obj->{contigs}} > 0, 'returned object has contigs!');
    ok (keys %$cid_hash, 'There is a contigID_hash:'.Dumper($cid_hash));

    lives_ok {
        ($ret_contig_obj, $cid_hash) = $annoutil->_get_contigs($obj_refseq_GCF);
    } "_get_contigs runs successfully on assembly $obj_refseq_GCF";
    ok (@{$ret_contig_obj->{contigs}} eq 0, 'returned object has no contigs!');
    ok (keys %$cid_hash == 0, 'no contigID_hash!');
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
=cut

#
## Global variables for the annotation process steps to share ##
#
my ($rast_ref00, $rast_ref01, $rast_ref02, $rast_ref1, $rast_ref2);
my ($inputgenome00, $inputgenome01, $inputgenome02, $inputgenome1, $inputgenome2);
my ($parameters00, $parameters01, $parameters02, $parameters1, $parameters2);

# Test _set_parameters_by_input with genome/assembly object refs in prod
subtest '_set_parameters_by_input' => sub {
    # genome objects in workspace #65386
    $parameters00 = {
         output_genome_name => 'test_out_gn_name00',
         scientific_name => 'Buchnera aphidicola',
         output_workspace => $ws_name,
         object_ref => $obj_refseq_GCF
    };
    # 00. creating default genome object
    $inputgenome00 = {
        id => $parameters00->{output_genome_name},
        genetic_code => $parameters00->{genetic_code},
        scientific_name => $parameters00->{scientific_name},
        domain => $parameters00->{domain},
        contigs => [],
        features => []
    };

    if ($parameters00->{ncbi_taxon_id}) {
        $inputgenome00->{taxon_assignments} = {
            'ncbi' => '' . $parameters00->{ncbi_taxon_id}};
    }
    my $expected_params00 = {
          scientific_name => $parameters00->{scientific_name},
          output_genome_name => $parameters00->{output_genome_name},
          output_workspace => $ws_name,
          object_ref => $parameters00->{object_ref},
          genetic_code => 11,
          domain => 'Bacteria'
    };

    lives_ok {
        ($rast_ref00, $inputgenome00) = $annoutil->_set_parameters_by_input(
                                            $parameters00, $inputgenome00);
    } "_set_parameters_by_input runs successfully on genome $obj_refseq_GCF";
    $parameters00 = $rast_ref00->{parameters};
    ok (@{$inputgenome00->{features}} == 410, "inputgenome has 410 feature(s).");
    ok (@{$rast_ref00->{contigobj}{contigs}} == 3, "inputgenome has 3 contig(s).");
    is ($inputgenome00->{assembly_ref}, '19217/1/7', 'inputgenome assembly_ref is correct.');

    # before merging with default gene call settings
    cmp_deeply($parameters00, $expected_params00, "parameters has expected input param values.");

    #  merge with the default gene call settings
    my $default_params = $annoutil->_set_default_parameters();
    $parameters00 = { %$default_params, %$parameters00 };
    $rast_ref00->{parameters} = $parameters00;

    my $expected_params001 = { %$default_params, %$expected_params00 };
    # after merging with default gene call settings
    cmp_deeply($parameters00, $expected_params001, "parameters has default workflows.");

    my $expected_contigID_hash = {
          'contigID_2' => 'NZ_LT667501.1',
          'contigID_1' => 'NZ_LT667500.1',
          'contigID_3' => 'NZ_LT667502.1'
        };

    cmp_deeply($rast_ref00->{contigID_hash}, $expected_contigID_hash,
               "contigID hash on $obj_refseq_GCF was generated correctly.");

    $parameters01 = {
         output_genome_name => 'test_out_gn_name01',
         scientific_name => 'Clostridium botulinum',
         output_workspace => $ws_name,
         object_ref => $obj_65386_1  # $obj_refseq_GCF
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
          output_workspace => $ws_name,
          object_ref => $parameters01->{object_ref},
          genetic_code => 11,
          domain => 'Bacteria'
    };

    lives_ok {
        ($rast_ref01, $inputgenome01) = $annoutil->_set_parameters_by_input(
                                            $parameters01, $inputgenome01);
    } "_set_parameters_by_input runs successfully on genome $obj_65386_1";
    $parameters01 = $rast_ref01->{parameters};
    ok (@{$inputgenome01->{features}} > 0,
        "inputgenome has ".scalar @{$inputgenome01->{features}}." feature(s).");
    ok (@{$rast_ref01->{contigobj}{contigs}} > 0,
        "inputgenome has ".scalar @{$rast_ref01->{contigobj}{contigs}} ." contig(s).");
    is ($inputgenome01->{assembly_ref}, '19217/360049/1', 'inputgenome assembly_ref is correct.');

    # before merging with default gene call settings
    cmp_deeply($parameters01, $expected_params01, "parameters has expected input param values.");

    #  merge with the default gene call settings
    $default_params = $annoutil->_set_default_parameters();
    $parameters01 = { %$default_params, %$parameters01 };
    $rast_ref01->{parameters} = $parameters01;

    my $expected_params011 = { %$default_params, %$expected_params01 };
    # after merging with default gene call settings
    cmp_deeply($parameters01, $expected_params011, "parameters has default workflows.");
    $expected_contigID_hash = {
          'contigID_1' => 'NZ_CP028859.1',
          'contigID_2' => 'NZ_CP028860.1'
    };
    cmp_deeply($rast_ref01->{contigID_hash}, $expected_contigID_hash,
               "contigID hash on $obj_65386_1 was generated correctly.");

    # another genome object in workspace #65386 that did not have a problem to run
    $parameters02 = {
         output_genome_name => 'test_out_gn_name02',
         scientific_name => 'Methanosarcina acetivorans C2A',
         output_workspace => $ws_name,
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
          output_workspace => $ws_name,
          object_ref => $parameters02->{object_ref},
          genetic_code => 11,
          domain => 'Archaea'
    };

    lives_ok {
        ($rast_ref02, $inputgenome02) = $annoutil->_set_parameters_by_input(
                                            $parameters02, $inputgenome02);
    } "_set_parameters_by_input runs successfully on genome $obj_65386_2";
    $parameters02 = $rast_ref02->{parameters};
    ok (@{$inputgenome02->{features}} > 0,
        "inputgenome has ".scalar @{$inputgenome02->{features}}." feature(s).");
    ok (@{$rast_ref02->{contigobj}{contigs}} > 0,
        "inputgenome has ".scalar @{$rast_ref02->{contigobj}{contigs}} ." contig(s).");
    is ($inputgenome02->{assembly_ref}, '19217/194865/2', 'inputgenome assembly_ref is correct.');

    # before merging with default gene call settings
    cmp_deeply($parameters02, $expected_params02, "parameters has expected input param values.");

    #  merge with the default gene call settings
    $parameters02 = { %$default_params, %$parameters02 };
    $rast_ref02->{parameters} = $parameters02;
    my $expected_params022 = { %$default_params, %$expected_params02};

    # after merging with default gene call settings
    cmp_deeply($parameters02, $expected_params022, "parameters has default workflows.");

    $expected_contigID_hash = {
          'contigID_1' => 'NC_003552.1'
    };
    cmp_deeply($rast_ref02->{contigID_hash}, $expected_contigID_hash,
               "contigID hash on $obj_65386_2 was generated correctly.");

    # a genome object
    $parameters1 = {
         output_genome_name => 'test_out_gn_name1',
         output_workspace => $ws_name,
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
          'output_workspace' => $ws_name,
          'object_ref' => $parameters1->{object_ref},
          'genetic_code' => 11,
          'domain' => 'Bacteria'
    };

    lives_ok {
        ($rast_ref1, $inputgenome1) = $annoutil->_set_parameters_by_input(
                                            $parameters1, $inputgenome1);
    } "_set_parameters_by_input runs successfully on genome $obj_Ecoli";
    $parameters1 = $rast_ref1->{parameters};

    ok (@{$inputgenome1->{features}} > 0,
        "inputgenome has ".scalar @{$inputgenome1->{features}}." feature(s).");
    ok (@{$rast_ref1->{contigobj}{contigs}} > 0,
        "inputgenome has ".scalar @{$rast_ref1->{contigobj}{contigs}} ." contig(s).");
    is ($inputgenome1->{assembly_ref}, '2901/78/1', 'inputgenome assembly_ref is correct.');

    # before merging with default gene call settings
    cmp_deeply($parameters1, $expected_params1, 'parameters has expected input param values.');

    # merge with the default gene call settings
    $default_params = $annoutil->_set_default_parameters();
    $parameters1 = { %$default_params, %$parameters1 };
    $rast_ref1->{parameters} = $parameters1;

    my $expected_params2 = { %$default_params, %$expected_params1 };

    # after merging with default gene call settings
    cmp_deeply($expected_params2, $parameters1, 'parameters has default workflows.');

    # an assembly object
    $parameters2 = {
         output_genome_name => 'test_out_gn_name2',
         output_workspace => $ws_name,
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
        object_ref => $obj_asmb,
        genetic_code => undef,
        output_genome_name => $parameters2->{output_genome_name},
        domain => undef,
        scientific_name => undef,
        output_workspace => $ws_name
    };
    lives_ok {
        ($rast_ref2, $inputgenome2) = $annoutil->_set_parameters_by_input(
                                            $parameters2, $inputgenome2);
    } "_set_parameters_by_input runs successfully on assembly $obj_asmb";
    $parameters2 = $rast_ref2->{parameters};

    # before merging with default gene call settings
    cmp_deeply($parameters2, $expected_params3, 'parameters has input values');

    # merge with the default gene call settings
    $parameters2 = { %$default_params, %$parameters2 };
    my $expected_params4 = { %$default_params, %$expected_params3 };
    $rast_ref2->{parameters} = $parameters2;

    # after merging with default gene call settings
    cmp_deeply($expected_params4, $parameters2, 'parameters has default workflows');

    ok (@{$inputgenome2->{features}} == 0, 'inputgenome (assembly) has NO features.');
    ok (@{$rast_ref2->{contigobj}{contigs}} > 0,
        "inputgenome has ".scalar @{$rast_ref2->{contigobj}{contigs}} ." contig(s).");
    is ($inputgenome2->{assembly_ref}, $obj_asmb, 'inputgenome assembly_ref is correct.');
};

# Test _set_messageNcontigs with genome/assembly object refs in prod
subtest '_set_messageNcontigs' => sub {
    # a refseq genome object in workspace #63171
    lives_ok {
        ($rast_ref00, $inputgenome00) = $annoutil->_set_messageNcontigs(
                                            $rast_ref00, $inputgenome00);
    } "_set_messageNcontigs runs successfully on genome $obj_refseq_GCF";
    $parameters00 = $rast_ref00->{parameters};
    my $msg00 = $rast_ref00->{message};
    my $tax00 = $rast_ref00->{tax_domain};

    ok (@{$inputgenome00->{features}} > 0,
        "inputgenome has ".scalar @{$inputgenome00->{features}}. " feature(s).");
    ok (@{$inputgenome00->{contigs}} > 0,
        "inputgenome has ".scalar @{$inputgenome00->{contigs}} . " contig(s).");
    ok (length($msg00) > 0, "Message for genome input has contents:\n$msg00");
    ok ($tax00 eq 'Bacteria', "tax_domain for genome input has value:$tax00");

    # a genome object in workspace #65386
    lives_ok {
        ($rast_ref01, $inputgenome01) = $annoutil->_set_messageNcontigs(
                                            $rast_ref01, $inputgenome01);
    } "_set_messageNcontigs runs successfully on genome $obj_65386_1";
    $parameters01 = $rast_ref01->{parameters};
    my $msg01 = $rast_ref01->{message};
    my $tax01 = $rast_ref01->{tax_domain};

    ok (@{$inputgenome01->{features}} > 0,
        "inputgenome has ".scalar @{$inputgenome01->{features}}. " feature(s).");
    ok (@{$inputgenome01->{contigs}} > 0,
        "inputgenome has ".scalar @{$inputgenome01->{contigs}} . " contig(s).");
    ok (length($msg01) > 0, "Message for genome input has contents:\n$msg01");
    ok ($tax01 eq 'Bacteria', "tax_domain for genome input has value:$tax01");

    # another genome object in workspace #65386
    lives_ok {
        ($rast_ref02, $inputgenome02) = $annoutil->_set_messageNcontigs(
                                            $rast_ref02, $inputgenome02);
    } "_set_messageNcontigs runs successfully on genome $obj_65386_2";

    $parameters02 = $rast_ref02->{parameters};
    my $msg02 = $rast_ref02->{message};
    my $tax02 = $rast_ref02->{tax_domain};

    ok (@{$inputgenome02->{features}} > 0,
        "inputgenome has ".scalar @{$inputgenome02->{features}}. " feature(s).");
    ok (@{$inputgenome02->{contigs}} > 0,
        "inputgenome has ".scalar @{$inputgenome02->{contigs}} . " contig(s).");
    ok (length($msg02) > 0, "Message for genome input has contents:\n$msg02");
    ok ($tax02 eq 'Archaea', "tax_domain for genome input has value:$tax02");

    # a genome object
    lives_ok {
        ($rast_ref1, $inputgenome1) = $annoutil->_set_messageNcontigs(
                                            $rast_ref1, $inputgenome1);
    } "_set_messageNcontigs runs successfully on genome $obj_Ecoli";
    $parameters1 = $rast_ref1->{parameters};
    my $msg1 = $rast_ref1->{message};
    my $tax1 = $rast_ref1->{tax_domain};

    ok (@{$inputgenome1->{features}} > 0,
        "inputgenome has ".scalar @{$inputgenome1->{features}}. " feature(s).");
    ok (@{$inputgenome1->{contigs}} > 0,
        "inputgenome has ".scalar @{$inputgenome1->{contigs}} . " contig(s).");
    ok (length($msg1) > 0, "Message for genome input has contents:\n$msg1");
    ok ($tax1 eq 'Bacteria', "tax_domain for genome input has value:$tax1");

    # an assembly object
    lives_ok {
        ($rast_ref2, $inputgenome2) = $annoutil->_set_messageNcontigs(
                                            $rast_ref2, $inputgenome2);
    } "_set_messageNcontigs runs successfully on assembly $obj_asmb";
    $parameters2 = $rast_ref2->{parameters};
    my $msg2 = $rast_ref2->{message};
    my $tax2 = $rast_ref2->{tax_domain};

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
            { 'name' => 'call_features_crispr' },
            { 'name' => 'call_features_CDS_glimmer3',
              'glimmer3_parameters' => {
                                         'min_training_len' => '2000'
                                       }
            },
            { 'name' => 'call_features_CDS_prodigal' }
        ]
    };

    # a refseq genome object in workspace #63171
    lives_ok {
        $rast_ref00 = $annoutil->_set_genecall_workflow(
                                   $rast_ref00, $inputgenome00);
    } "_set_genecall_workflow runs successfully on genome $obj_refseq_GCF";
    my $genecall_workflow00 = $rast_ref00->{genecall_workflow};
    cmp_deeply($genecall_workflow00, $exp_gc_workflow,
                'gc_workflow built correctly');

    # a genome object in workspace #65386
    lives_ok {
        $rast_ref01 = $annoutil->_set_genecall_workflow(
                                   $rast_ref01, $inputgenome01);
    } "_set_genecall_workflow runs successfully on genome $obj_65386_1";
    my $genecall_workflow01 = $rast_ref01->{genecall_workflow};
    cmp_deeply($genecall_workflow01, $exp_gc_workflow,
                'gc_workflow built correctly');

    # another genome object in workspace #65386
    lives_ok {
        $rast_ref02 = $annoutil->_set_genecall_workflow(
                                   $rast_ref02, $inputgenome02);
    } "_set_genecall_workflow runs successfully on genome $obj_65386_2";
    my $genecall_workflow02 = $rast_ref02->{genecall_workflow};
    cmp_deeply($genecall_workflow02, $exp_gc_workflow, 'gc_workflow built correctly');

    # a genome object
    lives_ok {
        $rast_ref1 = $annoutil->_set_genecall_workflow(
                                  $rast_ref1, $inputgenome1);
    } "_set_genecall_workflow runs successfully on genome $obj_Ecoli";
    my $genecall_workflow1 = $rast_ref1->{genecall_workflow};
    cmp_deeply($genecall_workflow1, $exp_gc_workflow, 'gc_workflow built correctly');

    # an assembly object
    my $exp_gc_workflow2 = {
        'stages' => [
            { 'name' => 'call_features_rRNA_SEED' },
            { 'name' => 'call_features_tRNA_trnascan' },
            { 'name' => 'call_features_repeat_region_SEED',
              'repeat_region_SEED_parameters' => {
                                                    'min_length' => '100',
                                                    'min_identity' => '95'
                                                 }
            },
            { 'name' => 'call_features_CDS_glimmer3',
              'glimmer3_parameters' => {
                                           'min_training_len' => '2000'
                                       }
            },
            { 'name' => 'call_features_CDS_prodigal'}
        ]
    };

    lives_ok {
        $rast_ref2 = $annoutil->_set_genecall_workflow(
                                  $rast_ref2, $inputgenome2);
    } "_set_genecall_workflow runs successfully on assembly $obj_asmb";
    my $msg2 = $rast_ref2->{message};
    my $genecall_workflow2 = $rast_ref2->{genecall_workflow};
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

    # a refseq genome object in workspace #63171
    lives_ok {
        $rast_ref00 = $annoutil->_set_annotation_workflow($rast_ref00);
    } "_set_genecall_workflow runs successfully on genome $obj_refseq_GCF";
    my $annomessage00 = $rast_ref00->{annomessage};
    my $annotate_workflow00 = $rast_ref00->{annotate_workflow};
    print "genome annotation-msg:\n$annomessage00";
    cmp_deeply($annotate_workflow00, $exp_ann_workflow, 'ann_workflow built correctly');

    # a genome object in workspace #65386
    lives_ok {
        $rast_ref01 = $annoutil->_set_annotation_workflow($rast_ref01);
    } "_set_genecall_workflow runs successfully on genome $obj_65386_1";
    my $annomessage01 = $rast_ref01->{annomessage};
    my $annotate_workflow01 = $rast_ref01->{annotate_workflow};
    print "genome annotation-msg:\n$annomessage01";
    cmp_deeply($annotate_workflow01, $exp_ann_workflow, 'ann_workflow built correctly');

    # another genome object in workspace #65386
    lives_ok {
        $rast_ref02 = $annoutil->_set_annotation_workflow($rast_ref02);
    } "_set_genecall_workflow runs successfully on genome $obj_65386_2";
    my $annomessage02 = $rast_ref02->{annomessage};
    my $annotate_workflow02 = $rast_ref02->{annotate_workflow};
    print "genome annotation-msg:\n$annomessage02";
    cmp_deeply($annotate_workflow02, $exp_ann_workflow, 'ann_workflow built correctly');

    # a genome object
    lives_ok {
        $rast_ref1 = $annoutil->_set_annotation_workflow($rast_ref1);
    } "_set_genecall_workflow runs successfully on genome $obj_Ecoli";
    my $annomessage1 = $rast_ref1->{annomessage};
    my $annotate_workflow1 = $rast_ref1->{annotate_workflow};
    print "genome annotation-msg:\n$annomessage1";
    cmp_deeply($annotate_workflow1, $exp_ann_workflow, 'ann_workflow built correctly');

    # an assembly object
    my $exp_ann_workflow_asmb = {
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
        $rast_ref2 = $annoutil->_set_annotation_workflow($rast_ref2);
    } "_set_genecall_workflow runs successfully on assembly $obj_asmb";
    my $annomessage2 = $rast_ref2->{annomessage};
    my $annotate_workflow2 = $rast_ref2->{annotate_workflow};
    print "assembly annotation-msg:\n$annomessage2";
    cmp_deeply($exp_ann_workflow_asmb, $annotate_workflow2,
                'ann_workflow built correctly for an assembly');
};

# Test _renumber_features with genome/assembly object refs in prod
subtest '_renumber_features' => sub {
    # a refseq genome object in workspace #63171
    lives_ok {
        ($rast_ref00, $inputgenome00) = $annoutil->_renumber_features(
                        $rast_ref00, $inputgenome00);
    } "_renumber_features runs successfully on genome $obj_refseq_GCF";
    my $msg00 = $rast_ref00->{message};
    print "genome merged-msg:\n$msg00";
    ok (@{$inputgenome00->{features}} == 0,
        "_renumber_features: inputgenome has NO features.");

    # a genome object in workspace #65386
    lives_ok {
        ($rast_ref01, $inputgenome01) = $annoutil->_renumber_features(
                        $rast_ref01, $inputgenome01);
    } "_renumber_features runs successfully on genome $obj_65386_1";
    my $msg01 = $rast_ref01->{message};
    print "genome merged-msg:\n$msg01";
    ok (@{$inputgenome01->{features}} == 0,
        "_renumber_features: inputgenome has NO features.");

    # another genome object in workspace #65386
    lives_ok {
        ($rast_ref02, $inputgenome02) = $annoutil->_renumber_features(
                        $rast_ref02, $inputgenome02);
    } "_renumber_features runs successfully on genome $obj_65386_2";
    my $msg02 = $rast_ref02->{message};
    print "genome merged-msg:\n$msg02";
    ok (@{$inputgenome02->{features}} == 0,
        "_renumber_features: inputgenome has NO features.");

    # a genome object
    lives_ok {
        ($rast_ref1, $inputgenome1) = $annoutil->_renumber_features(
                        $rast_ref1, $inputgenome1);
    } "_renumber_features runs successfully on genome $obj_Ecoli";
    my $msg1 = $rast_ref1->{message};
    print "genome merged-msg:\n$msg1";
    ok (@{$inputgenome1->{features}} > 0,
        "_renumber_features: inputgenome has ".scalar @{$inputgenome1->{features}}." features.");

    # an assembly object
    lives_ok {
        ($rast_ref2, $inputgenome2) = $annoutil->_renumber_features(
                        $rast_ref2, $inputgenome2);
    } "_renumber_features runs successfully on assembly $obj_asmb";
    ok (exists($rast_ref2->{renumber_workflow}), "renumber_workflow is created");
    ok (@{$inputgenome2->{features}} == 0,
        "_renumber_features: Assembly inputgenome has NO features.");
};

# Test _pre_rast_call with genome/assembly object refs in prod
subtest '_pre_rast_call' => sub {
    # a refseq genome object in workspace #63171
    lives_ok {
        ($rast_ref00, $inputgenome00) = $annoutil->_pre_rast_call(
                        $rast_ref00, $inputgenome00);
    } "_pre_rast_call runs successfully on genome $obj_refseq_GCF";

    ok (@{$inputgenome00->{features}} == 0,
        "_pre_rast_call: inputgenome has NO feature(s).");
    ok (@{$inputgenome00->{non_coding_features}} == 106,
        "_pre_rast_call: inputgenome has 106 non_coding_feature(s).");
    ok (keys %{ $rast_ref00->{genehash} }, "genehash created.");
    ok (@{$inputgenome00->{ontology_events}} == 1,
        "There is 1 ontology event for input genome.");
    ok (@{$rast_ref00->{contigobj}{contigs}} == 3,
        "_pre_rast_call: inputgenome has 3 contig(s).");

    # a genome object in workspace #65386
    lives_ok {
        ($rast_ref01, $inputgenome01) = $annoutil->_pre_rast_call(
                        $rast_ref01, $inputgenome01);
    } "_pre_rast_call runs successfully on genome $obj_65386_1";

    ok (@{$inputgenome01->{features}} == 0,
        "_pre_rast_call: inputgenome has NO feature(s).");
    ok (@{$inputgenome01->{non_coding_features}} == 255,
        "_pre_rast_call: inputgenome has 255 non_coding_feature(s).");
    ok (keys %{ $rast_ref01->{genehash} },
        "Gene hash created from genome with elements.");
    ok (@{$rast_ref01->{contigobj}{contigs}} == 2,
        "_pre_rast_call: inputgenome has 2 contig(s).");

    # another genome object in workspace #65386
    lives_ok {
        ($rast_ref02, $inputgenome02) = $annoutil->_pre_rast_call(
                        $rast_ref02, $inputgenome02);
    } "_pre_rast_call runs successfully on genome $obj_65386_2";

    ok (keys %{ $rast_ref02->{genehash} },
        "Gene hash created from genome with elements.");
    ok (@{$inputgenome02->{features}} == 0,
        "_pre_rast_call: inputgenome has NO feature(s).");
    ok (@{$inputgenome02->{non_coding_features}} == 138,
        "_pre_rast_call: inputgenome has 138 non_coding_feature(s).");
    ok (@{$rast_ref02->{contigobj}{contigs}} > 0,
        "_pre_rast_call: inputgenome has 2 contig(s).");

    # a genome object
    lives_ok {
        ($rast_ref1, $inputgenome1) = $annoutil->_pre_rast_call(
                        $rast_ref1, $inputgenome1);
    } "_pre_rast_call runs successfully on genome $obj_Ecoli";

    ok (keys %{ $rast_ref1->{genehash} }, "Gene hash created from genome $obj_Ecoli.");
    ok (@{$inputgenome1->{features}} == 358,
        "_pre_rast_call: inputgenome has 358 feature(s).");
    ok (@{$inputgenome1->{non_coding_features}} == 0,
        "_pre_rast_call: inputgenome has 0 non_coding_feature(s).");
    ok (@{$rast_ref1->{contigobj}{contigs}} == 1,
        "_pre_rast_call: inputgenome has 1 contig(s).");

    # an assembly object
    lives_ok {
        ($rast_ref2, $inputgenome2) = $annoutil->_pre_rast_call(
                        $rast_ref2, $inputgenome2);
    } "_pre_rast_call runs successfully on assembly $obj_asmb";

    ok ( scalar keys %{ $rast_ref2->{genehash} } == 0,
        "No genehash created from assembly $obj_asmb.");
    ok (@{$inputgenome2->{features}} == 0,
        "_pre_rast_call: inputgenome has NO feature(s).");
    ok (!defined($inputgenome2->{non_coding_features}),
        "_pre_rast_call: inputgenome has NO non_coding_feature(s).");
    ok (@{$rast_ref2->{contigobj}{contigs}} == 1,
        "_pre_rast_call: $obj_asmb has 1 contig.");
};


# Test _run_rast_workflow on annotation workflow with object refs in prod
my ($ann_genome00, $ann_genome01, $ann_genome02, $ann_genome1, $ann_genome2,
    $final_genome00, $final_genome01, $final_genome02, $final_genome1, $final_genome2);
subtest '_run_rast_workflow_ann' => sub {
    # a refseq genome object in workspace #63171
    lives_ok {
        $ann_genome00 = $annoutil->_run_rast_workflow(
              $inputgenome00, $rast_ref00->{annotate_workflow});
    } "_run_rast_workflow on annotation runs successfully on genome $obj_refseq_GCF";

    # another genome object in workspace #65386
    lives_ok {
        $ann_genome02 = $annoutil->_run_rast_workflow(
              $inputgenome02, $rast_ref02->{annotate_workflow});
    } "_run_rast_workflow on annotation runs successfully on genome $obj_65386_2";

    # a genome object in workspace #65386
    lives_ok {
        $ann_genome01 = $annoutil->_run_rast_workflow(
              $inputgenome01, $rast_ref01->{annotate_workflow});
    } "_run_rast_workflow on annotation runs successfully on genome $obj_65386_1";

    # another genome object in workspace #65386
    lives_ok {
        $ann_genome02 = $annoutil->_run_rast_workflow(
              $inputgenome02, $rast_ref02->{annotate_workflow});
    } "_run_rast_workflow on annotation runs successfully on genome $obj_65386_2";

    # a genome object
    lives_ok {
        $ann_genome1 = $annoutil->_run_rast_workflow(
              $inputgenome1, $rast_ref1->{annotate_workflow});
    } "_run_rast_workflow on annotation runs successfully on genome $obj_Ecoli";

    # an assembly object
    lives_ok {
        $ann_genome2 = $annoutil->_run_rast_workflow(
              $inputgenome2, $rast_ref2->{annotate_workflow});
    } "_run_rast_workflow returns normally on assembly $obj_asmb";
};


## test remapping of contig_ids##
subtest '_check_contigID_mapping' => sub {
    my $gn_pairs = [
        {gn => $ann_genome00, rd_ref => $rast_ref00},
        {gn => $ann_genome01, rd_ref => $rast_ref01},
        {gn => $ann_genome02, rd_ref => $rast_ref02},
        {gn => $ann_genome1, rd_ref => $rast_ref1},
        {gn => $ann_genome2, rd_ref => $rast_ref2}
    ];

    for my $p (@$gn_pairs) {
        my %rd = %{ $p->{rd_ref} };
        my $genome_before_remapping = $p->{gn};

        my $ctgID_hash = $rd{contigID_hash};
        ok( keys %$ctgID_hash, "There is a contigID_hash" );

        ## Run the remapping
        my $genome_after_remapping = $annoutil->_remap_contigIDs( $ctgID_hash, $genome_before_remapping );
        print "\n=========after remapping the rasted genome on $genome_before_remapping->{id}:\n";

        my $exp_ctg_ids = [];
        if( $genome_after_remapping->{contig_ids} ) {
            for my $cid (@{$genome_before_remapping->{contig_ids}}) {
                if( $ctgID_hash->{$cid} ) {
                    push @{$exp_ctg_ids}, $ctgID_hash->{$cid};
                } else {
                    push @{$exp_ctg_ids}, $cid;
                }
            }
            cmp_deeply $genome_after_remapping->{contig_ids}, $exp_ctg_ids,
                        'contig ids remapped correctly';
        }

        $exp_ctg_ids = [];
        my $result_ctg_ids = [];
        if( $genome_after_remapping->{contigs} ) {
            for my $ctg_before (@{$genome_before_remapping->{contigs}}) {
                if( $ctgID_hash->{$ctg_before->{id}} ) {
                    push @{$exp_ctg_ids}, $ctgID_hash->{$ctg_before->{id}};
                } else {
                    push @{$exp_ctg_ids}, $ctg_before->{id};
                }
            }
            for my $ctg_after (@{$genome_after_remapping->{contigs}}) {
                push @{$result_ctg_ids}, $ctg_after->{id};
            }
            cmp_deeply $result_ctg_ids, $exp_ctg_ids,
                        'contig ids of contigs array remapped correctly';
        }

        my ($exp_locs, $result_locs);
        for my $feature_type ( qw( features non_coding_features cdss mrnas ) ) {
            next unless $genome_after_remapping->{$feature_type};
            $result_locs = [];
            for my $ftr_after ( @{ $genome_after_remapping->{$feature_type} } ) {
                push @{$result_locs}, $ftr_after->{location};
            }
            $exp_locs = get_feature_locations($ctgID_hash, $genome_before_remapping->{$feature_type});

            cmp_deeply $result_locs, $exp_locs,
                        "contig ids of $feature_type locations remapped correctly";
        }
    }
};

# Test _post_rast_ann_call with genome/assembly object refs in prod
subtest '_post_rast_ann_call' => sub {
    # a refseq genome object in workspace #63171
    my $cnt = 0;
    my $ncoding_features00 = $ann_genome00->{non_coding_features};
    for my $ncoding_ftr (@{$ncoding_features00}) {
        if(exists($ncoding_ftr->{type})) {
            $cnt++;
        }
    }
    print "**for $obj_refseq_GCF: Total count of non_coding_features BEFORE _post_rast_ann_call:".scalar @{$ncoding_features00}."\n";
    print "**for $obj_refseq_GCF:Count of non_coding_features WITH 'type' BEFORE _post_rast_ann_call:$cnt\n";

    lives_ok {
        $final_genome00 = $annoutil->_post_rast_ann_call(
                          $ann_genome00, $inputgenome00, $rast_ref00);
    } "_post_rast_ann_call runs successfully on genome $obj_refseq_GCF";
    # Focus on the 'non_coding_features' type field check
    $cnt = 0;
    $ncoding_features00 = $final_genome00->{non_coding_features};
    for my $ncoding_ftr (@{$ncoding_features00}) {
        if(exists($ncoding_ftr->{type})) {
            $cnt++;
            # print "type value: $ncoding_ftr->{type}\n";
        }
    }
    print "**for $obj_refseq_GCF:Count of non_coding_features WITH 'type' AFTER _post_rast_ann_call:$cnt\n";

    # a genome object in workspace #65386
    $cnt = 0;
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
                          $ann_genome01, $inputgenome01, $rast_ref01);
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
                          $ann_genome02, $inputgenome02, $rast_ref02);
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
                          $ann_genome1, $inputgenome1, $rast_ref1);
    } "_post_rast_ann_call runs successfully on genome $obj_Ecoli";

    # an assembly object
    lives_ok {
        $final_genome2 = $annoutil->_post_rast_ann_call(
                          $ann_genome2, $inputgenome2, $rast_ref2);
    } "_post_rast_ann_call runs successfully on assembly $obj_asmb";
};


# Test _build_seed_ontology with genome/assembly object refs in prod
subtest '_build_seed_ontology' => sub {
    # a refseq genome object in workspace #63171
    my $nc_ftr_count = @{$final_genome00->{non_coding_features}};
    print "\n********For case $obj_refseq_GCF*********\nBEFORE _build_seed_ontology there are $nc_ftr_count non_coding features.\n";

    lives_ok {
        ($final_genome00, $rast_ref00) = $annoutil->_build_seed_ontology(
              $rast_ref00, $final_genome00, $inputgenome00);
    } "_build_seed_ontology returns normally on genome $obj_refseq_GCF";

    $nc_ftr_count = @{$final_genome00->{non_coding_features}};
    print "\n********For case $obj_refseq_GCF*********\nAFTER _build_seed_ontology there are $nc_ftr_count non_coding features.\n";

    # a genome object in workspace #65386
    lives_ok {
        ($final_genome01, $rast_ref01) = $annoutil->_build_seed_ontology(
              $rast_ref01, $final_genome01, $inputgenome01);
    } "_build_seed_ontology returns normally on genome $obj_65386_1";

    # another genome object in workspace #65386
    lives_ok {
        ($final_genome02, $rast_ref02) = $annoutil->_build_seed_ontology(
              $rast_ref02, $final_genome02, $inputgenome02);
    } "_build_seed_ontology returns normally on genome $obj_65386_2";

    # a genome object
    lives_ok {
        ($final_genome1, $rast_ref1) = $annoutil->_build_seed_ontology(
              $rast_ref1, $final_genome1, $inputgenome1);
    } "_build_seed_ontology returns normally on genome $obj_Ecoli";

    # an assembly object
    $nc_ftr_count = @{$final_genome2->{non_coding_features}};
    print "\n********For case $obj_asmb*********\nBEFORE _build_seed_ontology there are $nc_ftr_count non_coding features.\n";

    lives_ok {
        ($final_genome2, $rast_ref2) = $annoutil->_build_seed_ontology(
              $rast_ref2, $final_genome2, $inputgenome2);
    } "_build_seed_ontology returns on assembly $obj_asmb";

    $nc_ftr_count = @{$final_genome2->{non_coding_features}};
    print "\n********For case $obj_asmb*********\nAFTER _build_seed_ontology there are $nc_ftr_count non_coding features.\n";
};


# Test _summarize_annotation with genome/assembly object refs in prod
subtest '_summarize_annotation' => sub {
    # a refseq genome object in workspace #63171
    my $nc_ftr_count = @{$final_genome00->{non_coding_features}};
    print "\n********For case $obj_refseq_GCF*********\nBEFORE _summarize_annotation there are $nc_ftr_count non_coding features.\n";

    lives_ok {
        ($final_genome00, $rast_ref00) = $annoutil->_summarize_annotation(
              $rast_ref00, $final_genome00, $inputgenome00);
    } "_summarize_annotation runs successfully on genome $obj_refseq_GCF";

    $nc_ftr_count = @{$final_genome00->{non_coding_features}};
    print "\n********For case $obj_refseq_GCF*********\nAFTER _summarize_annotation there are $nc_ftr_count non_coding features.\n";

    # a genome object in workspace #65386
    lives_ok {
        ($final_genome01, $rast_ref01) = $annoutil->_summarize_annotation(
              $rast_ref01, $final_genome01, $inputgenome01);
    } "_summarize_annotation runs successfully on genome $obj_65386_1";

    # another genome object in workspace #65386
    lives_ok {
        ($final_genome02, $rast_ref02) = $annoutil->_summarize_annotation(
              $rast_ref02, $final_genome02, $inputgenome02);
    } "_summarize_annotation runs successfully on genome $obj_65386_2";

    # a genome object
    lives_ok {
        ($final_genome1, $rast_ref1) = $annoutil->_summarize_annotation(
              $rast_ref1, $final_genome1, $inputgenome1);
    } "_summarize_annotation runs successfully on genome $obj_Ecoli";

    # an assembly object
    lives_ok {
        ($final_genome2, $rast_ref2) = $annoutil->_summarize_annotation(
              $rast_ref2, $final_genome2, $inputgenome2);
    } "_summarize_annotation runs successfully on assembly $obj_asmb";
};

=begin
# Test _reformat_feature_aliases for RAST annotated objects in prod
subtest '_reformat_feature_aliases' => sub {
    my $ret_cds;

    # a refseq genome object in workspace #63171
    foreach my $cds (@{$final_genome00->{cdss}}[0..2]) {
        ok (!defined($cds->{aliases}), "Aliases data is undef, so nothing will be changed.");
        $ret_cds = $annoutil->_reformat_feature_aliases($cds);
        cmp_deeply $ret_cds, $cds, "No change.";
    }

    # a genome object in workspace #65386
    foreach my $cds (@{$final_genome01->{cdss}}[0..2]) {
        ok (ref($cds->{aliases}->[0]) =~ /ARRAY/, "Aliases data is array of arrays, so nothing will be changed.");
        $ret_cds = $annoutil->_reformat_feature_aliases($cds);
        cmp_deeply $ret_cds, $cds, "no change.";
    }

    # a genome object in workspace #65386
    foreach my $cds (@{$final_genome02->{cdss}}[0..2]) {
        ok (ref($cds->{aliases}->[0]) =~ /ARRAY/, "Aliases data is array of arrays, so nothing will be changed.");
        $ret_cds = $annoutil->_reformat_feature_aliases($cds);
        cmp_deeply $ret_cds, $cds, "no change.";
    }

    # a genome object
    foreach my $cds (@{$final_genome1->{cdss}}[0..2]) {
        print "Before: the aliases data:\n".Dumper($cds->{aliases});
        ok (ref($cds->{aliases}->[0]) !~ /ARRAY/, "Aliases data is array of strings, so it will be changed.");
		my $expected = [];
		for my $als (@{$cds->{aliases}})  {
			my @ary = ('alias', $als);
			push(@{$expected}, \@ary);
		}
        $ret_cds = $annoutil->_reformat_feature_aliases($cds);
        cmp_deeply $ret_cds->{aliases}, $expected, "Aliases changed correctly.";
    }

    # rasted from an assembly object
    foreach my $cds (@{$final_genome2->{cdss}}[0..2]) {
        ok (!defined($cds->{aliases}), "Because of no aliases data, nothing will be changed.");
        $ret_cds = $annoutil->_reformat_feature_aliases($cds);
        cmp_deeply $ret_cds, $cds, "no change.";
    }
};


# Test _fillRequiredFields for RAST annotated objects in prod
subtest '_fillRequiredFields' => sub {
    # a refseq genome object in workspace #63171
    my $genome_clone = clone($final_genome00);
    my $ret_gn = $annoutil->_fillRequiredFields($genome_clone);
    my @this_not_that = ();
    foreach (keys %$ret_gn) {
        push(@this_not_that, $_) unless exists $final_genome00->{$_};
    }
    ok (!@this_not_that, "No new fields added.");

    # a genome object in workspace #65386
    @this_not_that = ();
    $genome_clone = clone($final_genome01);
    $ret_gn = $annoutil->_fillRequiredFields($genome_clone);
    foreach (keys %$ret_gn) {
        push(@this_not_that, $_) unless exists $final_genome01->{$_};
    }
    ok (!@this_not_that, "No new fields added.");

    # a genome object in workspace #65386
    @this_not_that = ();
    $genome_clone = clone($final_genome02);
    $ret_gn = $annoutil->_fillRequiredFields($genome_clone);
    foreach (keys %$ret_gn) {
        push(@this_not_that, $_) unless exists $final_genome02->{$_};
    }
    ok (!@this_not_that, "No new fields added.");
    print "Diff_arr=".Dumper(\@this_not_that);

    # a genome object
    @this_not_that = ();
    $genome_clone = clone($final_genome1);
    $ret_gn = $annoutil->_fillRequiredFields($genome_clone);
    foreach (keys %$ret_gn) {
        push(@this_not_that, $_) unless exists $final_genome1->{$_};
    }
    ok (@this_not_that, "New fields added:".Dumper(\@this_not_that));

    # rasted from an assembly object
    @this_not_that = ();
    $genome_clone = clone($final_genome2);
    $ret_gn = $annoutil->_fillRequiredFields($genome_clone);
    foreach (keys %$ret_gn) {
        push(@this_not_that, $_) unless exists $final_genome2->{$_};
    }
    ok (@this_not_that, "New fields added:".Dumper(\@this_not_that));
};


# Test _create_onto_terms for RAST annotated objects in prod
subtest '_create_onto_terms' => sub {
    # a refseq genome object in workspace #63171
    my $genome_clone = clone($final_genome00);
    my $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'features');
    my $termSize = keys %$ret_terms;
    ok (!$termSize, "No features to create ontology terms.");

    $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'cdss');
    $termSize = keys %$ret_terms;
    ok ($termSize, "$termSize cdss ontology terms created.");

    $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'non_coding_features');
    $termSize = keys %$ret_terms;
    ok ($termSize, "$termSize non_coding_features ontology terms created.");

    # a genome object in workspace #65386
    $genome_clone = clone($final_genome01);
    $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'features');
    $termSize = keys %$ret_terms;
    ok (!$termSize, "No features to create ontology terms.");

    $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'cdss');
    $termSize = keys %$ret_terms;
    ok ($termSize, "$termSize cdss ontology terms created.");

    $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'non_coding_features');
    $termSize = keys %$ret_terms;
    ok ($termSize, "$termSize non_coding_features ontology terms created.");

    # a genome object in workspace #65386
    $genome_clone = clone($final_genome02);
    $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'features');
    $termSize = keys %$ret_terms;
    ok (!$termSize, "No features to create ontology terms.");

    $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'cdss');
    $termSize = keys %$ret_terms;
    ok ($termSize, "$termSize cdss ontology terms created.");

    $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'non_coding_features');
    $termSize = keys %$ret_terms;
    ok ($termSize, "$termSize non_coding_features ontology terms created.");

    # a genome object
    $genome_clone = clone($final_genome1);
    my $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'features');
    $termSize = keys %$ret_terms;
    ok (!$termSize, "No features to create ontology terms.");

    $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'cdss');
    $termSize = keys %$ret_terms;
    ok ($termSize, "$termSize cdss ontology terms created.");

    $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'non_coding_features');
    $termSize = keys %$ret_terms;
    ok (!$termSize, "$termSize non_coding_features ontology terms created.");

    # rasted from an assembly object
    $genome_clone = clone($final_genome2);
    $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'features');
    $termSize = keys %$ret_terms;
    ok (!$termSize, "No features to create ontology terms.");

    $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'cdss');
    $termSize = keys %$ret_terms;
    ok (!$termSize, "$termSize cdss ontology terms created.");

    $ret_terms = $annoutil->_create_onto_terms($genome_clone, 'non_coding_features');
    $termSize = keys %$ret_terms;
    ok (!$termSize, "$termSize non_coding_features ontology terms created.");
};


# Test _build_ontology_events for RAST annotated objects in prod
subtest '_build_ontology_events' => sub {
    my ($ret_gn, $pm, $evts, $genome_clone);

    # a refseq genome object in workspace #63171
    $genome_clone = clone($final_genome00);
    $pm = $rast_ref00->{parameters};
    lives_ok {
        $ret_gn = $annoutil->_build_ontology_events($genome_clone, $pm);
    } "_build_ontology_events on rasted genome from $obj_refseq_GCF";
    $evts = $ret_gn->{events};
    ok( $evts, "_build_ontology_events returns events.");

    # a genome object in workspace #65386
    $genome_clone = clone($final_genome01);
    $pm = $rast_ref01->{parameters};
    lives_ok {
        $ret_gn = $annoutil->_build_ontology_events($genome_clone, $pm);
    } "_build_ontology_events on rasted genome from $obj_65386_1";
    $evts = $ret_gn->{events};
    ok( $evts, "_build_ontology_events returns events.");

    # another genome object in workspace #65386
    $genome_clone = clone($final_genome02);
    $pm = $rast_ref02->{parameters};
    lives_ok {
        $ret_gn = $annoutil->_build_ontology_events($genome_clone, $pm);
    } "_build_ontology_events on rasted genome from $obj_65386_2";
    $evts = $ret_gn->{events};
    ok( $evts, "_build_ontology_events returns events.");

    # a genome object
    $genome_clone = clone($final_genome1);
    $pm = $rast_ref1->{parameters};
    lives_ok {
        $ret_gn = $annoutil->_build_ontology_events($genome_clone, $pm);
    } "_build_ontology_events on rasted genome from $obj_Ecoli";
    $evts = $ret_gn->{events};
    ok( $evts, "_build_ontology_events returns events.");

    # from an assembly object
    $genome_clone = clone($final_genome2);
    $pm = $rast_ref2->{parameters};
    lives_ok {
        $ret_gn = $annoutil->_build_ontology_events($genome_clone, $pm);
    } "_build_ontology_events on rasted genome annotated from $obj_asmb";
    $evts = $ret_gn->{events};
    ok( $evts, "_build_ontology_events returns events.");
};
=cut

# Test _save_annotation_results with genome/assembly object refs in prod
subtest '_save_annotation_results' => sub {
    my ($save_ret, $out_msg);

    # a refseq genome object in workspace #63171
    my $nc_ftr_count = @{$final_genome00->{non_coding_features}};
    print "\n********For case $obj_refseq_GCF*********\n".
          "BEFORE _save_annotation_results there are $nc_ftr_count non_coding features.\n";

    my $ncoding_features00 = $final_genome00->{non_coding_features};
    my $cnt = 0;
    for my $ncoding_ftr (@{$ncoding_features00}) {
        if(exists($ncoding_ftr->{type})) {
            $cnt++;
            if ($cnt < 20) {
                print "type value: $ncoding_ftr->{type}\n";
            }
        }
    }
    ok ($nc_ftr_count==$cnt, "All $cnt non-coding features have defined field of type.\n");
    #print "***BEFORE OntSer saving, RAST annotation object data****\n".Dumper($final_genome00);

    lives_ok {
        ($save_ret, $out_msg) = $annoutil->_save_annotation_results(
                                            $final_genome00, $rast_ref00);
    } "_save_annotation_results finished on genome $obj_refseq_GCF and returned results.";
    ok (exists($save_ret->{ref}), "_save_annotation_results succeeded.");
    #my $saved_anno_data = $annoutil->_fetch_object_data($save_ret->{ref});
    #print "***One OntSer saved RAST annotation object data****\n".Dumper($saved_anno_data);

    # a genome object in workspace #65386
    $nc_ftr_count = @{$final_genome01->{non_coding_features}};
    lives_ok {
        ($save_ret, $out_msg) = $annoutil->_save_annotation_results(
                                            $final_genome01, $rast_ref01);
    } "_save_annotation_results on genome $obj_65386_1 returned as expected.";
    ok (exists($save_ret->{ref}), "_save_annotation_results succeeded.");

    # another genome object in workspace #65386
    lives_ok {
        ($save_ret, $out_msg) = $annoutil->_save_annotation_results(
                                            $final_genome02, $rast_ref02);
    } "_save_annotation_results finished on genome $obj_65386_2 returned as expected.";
    ok (exists($save_ret->{ref}), "_save_annotation_results succeeded.");

    # a genome object
    lives_ok {
        ($save_ret, $out_msg) = $annoutil->_save_annotation_results(
                                            $final_genome1, $rast_ref1);
    } "_save_annotation_results finished on genome $obj_Ecoli returned as expected";
    ok (exists($save_ret->{ref}), "_save_annotation_results succeeded.");

    # RAST annotation of an assembly object
    #print "******Before saving, the RASTed genome data:\n".Dumper($final_genome2);
    lives_ok {
        ($save_ret, $out_msg) = $annoutil->_save_annotation_results(
                                            $final_genome2, $rast_ref2);
    } "_save_annotation_results on assembly $obj_asmb returned expected result.";
    ok (exists($save_ret->{ref}), "_save_annotation_results succeeded.");
    my $saved_anno_data = $annoutil->_fetch_object_data($save_ret->{ref});
    #print "***One OntSer saved RAST annotation object data****\n".Dumper($saved_anno_data);
};

=begin
#
## variables for testing _build_workflows, _run_rast_workflow
## _run_rast_genecalls and _pre_rast_call
my ($gc_wf_ret0, $gc_inputgenome0, $gc_wf0,
    $gc_wf_ret1, $gc_inputgenome1, $gc_wf1,
    $gc_wf_ret2, $gc_inputgenome2, $gc_wf2);

# Test _build_workflows with genome/assembly object refs in prod
subtest '_build_workflows' => sub {
    # a genome object from https://narrative.kbase.us/narrative/ws.65386.obj.104
    my $params0 = {
         output_genome_name => 'build_gcwf_name0',
         output_workspace => $ws_name,
         object_ref => $obj_65386_1
    };
    lives_ok {
        ($gc_wf_ret0, $gc_inputgenome0) = $annoutil->_build_workflows($params0);
    } "_build_workflows returns normally on $obj_65386_1\n";
    $gc_wf0 = $gc_wf_ret0->{genecall_workflow};
    ok (exists($gc_wf_ret0->{genecall_workflow}), "generate genecall workflow:\n".Dumper($gc_wf0));
    ok (@{$gc_inputgenome0->{features}} == 0, 'inputgenome has NO features.');
    ok (@{$gc_inputgenome0->{contigs}} > 0,
        "inputgenome has ".scalar @{$gc_inputgenome0->{contigs}}." contig(s).");

    my $params1 = {
         output_genome_name => 'build_gcwf_gn_name',
         output_workspace => $ws_name,
         object_ref => $obj_Ecoli
    };

    lives_ok {
        ($gc_wf_ret1, $gc_inputgenome1) = $annoutil->_build_workflows($params1);
    } "_build_workflows returns normally on $obj_Ecoli\n";
    $gc_wf1 = $gc_wf_ret1->{genecall_workflow};
    ok (exists($gc_wf_ret1->{genecall_workflow}), "generate genecall workflow:\n".Dumper($gc_wf1));
    ok (@{$gc_inputgenome1->{features}} > 0,
        "inputgenome has ".scalar @{$gc_inputgenome1->{features}}." features.");
    ok (@{$gc_inputgenome1->{contigs}} > 0,
        "inputgenome has ".scalar @{$gc_inputgenome1->{contigs}}." contig(s).");

    my $params2 = {
         output_genome_name => 'build_gcwf_asmb_name',
         output_workspace => $ws_name,
         object_ref => $obj_asmb
    };
    lives_ok {
        ($gc_wf_ret2, $gc_inputgenome2) = $annoutil->_build_workflows($params2);
    } "_build_workflows returns normally on $obj_asmb\n";
    $gc_wf2 = $gc_wf_ret2->{genecall_workflow};
    ok (exists($gc_wf_ret2->{genecall_workflow}), "generate genecall workflow:\n".Dumper($gc_wf2));
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
    my $gc_gn0;
    lives_ok {
        $gc_gn0 = $annoutil->_run_rast_genecalls($gc_inputgenome0, $gc_wf0);
    } '_run_rast_genecalls returns ERROR due to kmer data absence or other causes.';
    ok (@{$gc_gn0->{features}} == 0, "Returned genome has NO features.");
    cmp_deeply($gc_gn0, $gc_inputgenome0, 'rast workflow will not run locally');

    # a genome object
    my $gc_gn1;
    lives_ok {
        $gc_gn1 = $annoutil->_run_rast_genecalls($gc_inputgenome1, $gc_wf1);
    } '_run_rast_genecalls returns ERROR due to kmer data absence or other causes.';
    ok (@{$gc_gn1->{features}} > 0, "Returned genome has features.");
    cmp_deeply($gc_gn1, $gc_inputgenome1, 'rast workflow will not run locally');

    # an assembly object
    my $gc_gn2;
    lives_ok {
        $gc_gn2 = $annoutil->_run_rast_genecalls($gc_inputgenome2, $gc_wf2);
    } '_run_rast_genecalls returns ERROR due to kmer data absence or other causes.';
    ok (@{$gc_gn2->{features}} == 0, "Returned genome has NO features.");
    cmp_deeply($gc_gn2, $gc_inputgenome2, 'rast workflow will not run locally');
};


subtest '_validate_KB_objref_name' => sub {
	my $obj = 'qzhang:narrative_1581052755332/short_one_metagenome';
	my $pass_test = $annoutil->_validate_KB_objref_name($obj);
	ok ($pass_test->{check_passed}, "$obj is a valid workspace object.\n");
	ok ($pass_test->{is_ref}, "$obj is a valid workspace object reference.\n");

	$obj = 'qzhang:narrative_1581052755332/short_one_metagenome/1';
	$pass_test = $annoutil->_validate_KB_objref_name($obj);
	ok ($pass_test->{check_passed}, "$obj is a valid workspace object.\n");
	ok ($pass_test->{is_ref}, "$obj is a valid workspace object reference.\n");

	$obj = '52755332/short_one/1';
	$pass_test = $annoutil->_validate_KB_objref_name($obj);
	ok ($pass_test->{check_passed}, "$obj is a valid workspace object.\n");
	ok ($pass_test->{is_ref}, "$obj is a valid workspace object reference.\n");

	$obj = '5332/345/3';
	$pass_test = $annoutil->_validate_KB_objref_name($obj);
	ok ($pass_test->{check_passed}, "$obj is a valid workspace object.\n");
	ok ($pass_test->{is_ref}, "$obj is a valid workspace object reference.\n");

	$obj = '5332/3wda9/123';
	$pass_test = $annoutil->_validate_KB_objref_name($obj);
	ok ($pass_test->{check_passed}, "$obj is a valid workspace object.\n");
	ok ($pass_test->{is_ref}, "$obj is a valid workspace object reference.\n");

	$obj = '5332/39';
	$pass_test = $annoutil->_validate_KB_objref_name($obj);
	ok ($pass_test->{check_passed}, "$obj is a valid workspace object.\n");
	ok ($pass_test->{is_ref}, "$obj is a valid workspace object reference.\n");

	$obj = '5332/3wda9';
	$pass_test = $annoutil->_validate_KB_objref_name($obj);
	ok ($pass_test->{check_passed}, "$obj is a valid workspace object.\n");
	ok ($pass_test->{is_ref}, "$obj is a valid workspace object reference.\n");

	$obj = '5332/3wda9/a';
	$pass_test = $annoutil->_validate_KB_objref_name($obj);
	ok ($pass_test->{check_passed} == 0 && $pass_test->{is_ref} == 0 &&
            $pass_test->{is_name} == 0,
            "$obj failed workspace object id test.\n");

	$obj = 'abc/def12/f';
	$pass_test = $annoutil->_validate_KB_objref_name($obj);
	ok ($pass_test->{check_passed} == 0 && $pass_test->{is_ref} == 0 &&
            $pass_test->{is_name} == 0,
            "$obj failed workspace object id test.\n");

	$obj = 'Clostridium_old.RAST';
	$pass_test = $annoutil->_validate_KB_objref_name($obj);
	ok ($pass_test->{check_passed}, "$obj is a valid workspace object.\n");
	ok ($pass_test->{is_name}, "$obj is a valid workspace object name.\n");

	$obj = '123/abc/|fg/8';
	$pass_test = $annoutil->_validate_KB_objref_name($obj);
	ok ($pass_test->{check_passed} == 0, "$obj is an invalid workspace object.\n");
	ok ($pass_test->{is_name} == 0, "$obj is an invalid workspace object name.\n");
};


## Combine several workflows into one
subtest '_combine_workflows' => sub {
    my $wf1 = {stages=>[{name => "call_features_rRNA_SEED"},
                        {name => "call_features_tRNA_trnascan"}]};
    my $wf2 = {stages=>[{
                         name => "call_features_repeat_region_SEED",
                         "repeat_region_SEED_parameters" => {
                                 "min_identity" => "95",
                                 "min_length" => "100"}}]};
    my $wf3 = {stages=>[{name => "call_features_CDS_glimmer3",
                         "glimmer3_parameters" => {
                                 "min_training_len" => "2000"}}]};
    my $wf4 = {stages=>[{
                         name => "annotate_proteins_kmer_v2",
                         "kmer_v2_parameters" => {
                                 "min_hits" => "5",
                                 "annotate_hypothetical_only" => 1}}]};
    my $wf5 = {stages=>[{name => "renumber_features"}]};

    my $wf = {};
    $wf->{genecall_workflow} = $wf1;
    $wf->{annotate_workflow} = $wf2;
    $wf->{renumber_workflow} = $wf5;
    my $exp_wf1 = {stages=>[]};
    push @{$exp_wf1->{stages}}, @{$wf1->{stages}};
    push @{$exp_wf1->{stages}}, @{$wf2->{stages}};
    push @{$exp_wf1->{stages}}, @{$wf5->{stages}};
    my $ret_wf = $annoutil->_combine_workflows($wf);
    cmp_deeply($exp_wf1, $ret_wf, "workflow1 combined as expected.");

    $wf->{genecall_workflow} = $wf3;
    $wf->{annotate_workflow} = $wf4;
    $wf->{renumber_workflow} = undef;
    my $exp_wf2 = {stages=>[]};
    push @{$exp_wf2->{stages}}, @{$wf3->{stages}};
    push @{$exp_wf2->{stages}}, @{$wf4->{stages}};
    $ret_wf = $annoutil->_combine_workflows($wf);
    cmp_deeply($exp_wf2, $ret_wf, "workflow2 combined as expected.");

    $wf->{genecall_workflow} = $wf1;
    $wf->{annotate_workflow} = $wf4;
    $wf->{renumber_workflow} = {stages=>[]};
    my $exp_wf3 = {stages=>[]};
    push @{$exp_wf3->{stages}}, @{$wf1->{stages}};
    push @{$exp_wf3->{stages}}, @{$wf4->{stages}};
    $ret_wf = $annoutil->_combine_workflows($wf);
    cmp_deeply($exp_wf3, $ret_wf, "workflow3 combined as expected.");
};

## _fetch_object_info with additonal argument added
subtest '_fetch_object_info' => sub {
    my $obj = '63171/441/1';
    my $pass_test = $annoutil->_validate_KB_objref_name($obj);
    ok ($pass_test->{check_passed}, "$obj passed id format check.\n");
    ok ($pass_test->{is_ref}, "$obj is a valid workspace object reference.\n");
    my $ret = $annoutil->_fetch_object_info($obj, $pass_test);
    is ( $ret->[0], 441, "correct object id found");
    is ( $ret->[4], 1, "correct object version found");
    is ( $ret->[6], 63171, "correct workspace id found");
    is ( $ret->[7], 'qzhang:narrative_1590705141087', "correct workspace found");

    my $ws = $ret->[7];
    my $obj1 = 'Carsonella_assembly';
    $pass_test = $annoutil->_validate_KB_objref_name($obj1);
    ok ($pass_test->{check_passed}, "$obj1 passed id format check.\n");
    ok ($pass_test->{is_name}, "$obj1 is a valid workspace object name.\n");
    $ret = $annoutil->_fetch_object_info($obj1, $pass_test, $ws);
    is ( $ret->[0], 243, "correct object id found");
    is ( $ret->[4], 1, "correct object version found");
    is ( $ret->[6], 63171, "correct workspace id found");
    is ( $ret->[7], $ws, "correct workspace found");

    my $obj2 = 'Clostridium_old.RAST';
    $pass_test = $annoutil->_validate_KB_objref_name($obj1);
    ok ($pass_test->{check_passed}, "$obj2 passed id format check.\n");
    ok ($pass_test->{is_name}, "$obj2 is a valid workspace object name.\n");
    $ret = $annoutil->_fetch_object_info($obj2, $pass_test, $ws);
    is ( $ret, undef, "object not found in workspace $ws, so return undef");
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
        my $p = {output_workspace => $ws_name,
                 output_genome_name => $outgn_name};
        print "input parameter=\n". Dumper($p);
        $annoutil->_check_annotation_params($p)
    } qr/$req2/,
        '_check_annotation_params dies with no object_ref';

    throws_ok {
        my $p = {output_workspace => $ws_name,
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
            {workspace => $ws_name,
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
            {output_workspace => $ws_name,
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
             output_workspace => $ws_name,
             output_genome_name => $outgn_name,
             object_ref => '456/1/2'};
    my $set_default_ok = '_check_annotation_params sets the default value for output_genome_name.';

    my $ret = $annoutil->_check_annotation_params(
            {output_workspace => $ws_name,
             output_genome_name => undef,
             object_ref => '456/1/2'});
    ok ($ret->{output_genome_name} eq $expected->{output_genome_name},
        'When undefined, '.$set_default_ok);

    $ret = $annoutil->_check_annotation_params(
            {output_workspace => $ws_name,
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
        "output_workspace" => $ws_name
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

    $params->{input_genomes} = []; # array of a prod object
    $params->{input_text} = '48109/9/1; 123/4/5;\t6789/12/8; |\n897/65/2';
    my $exp_input_text = '48109/9/1;123/4/5;6789/12/8;|;897/65/2';
    my $ret_parms5a = $annoutil->_check_bulk_annotation_params($params);
    print $ret_parms5a->{input_text};
    ok ($ret_parms5a->{input_text} eq $exp_input_text,
        "Got rid of spaces in the string of input_text.");

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

    my $temp_gff_file = catfile($scratch, 'temp_gff_file.gff');
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

    my $gff_contents3;
    lives_ok {
        ($gff_contents3, $attr_delimiter) = $annoutil->_parse_gff($Echinacea_gff_file, '=');
    } "Testing _parse_gff on file $Echinacea_gff_file succeeded.";
    ok( @{$gff_contents3} >0, "Parsing GFF on $Echinacea_gff_file returns result.\n");
    print "Parsed ". scalar @{$gff_contents3}." GFF contents.\n";

    is_deeply($gff_contents1, $gff_contents3, 'GFF data structures should be the same 2!');
};


# testing get the gff from a genome using obj ids from prod ONLY
subtest '_write_gff_from_genome' => sub {
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

    print "First 20 lines of the GFF file:\n";
    $annoutil->_print_fasta_gff(0, 20, $gff_fpath);
};

subtest 'anno_utils_rast_genome' => sub {
    # testing anno_utils rast_genome using obj ids from prod ONLY
    my ($parms, $rast_ret, $genome_obj, $rast_gn_data);
    $parms = {
        "object_ref" => $obj_Ecoli,
        "output_genome_name" => "rasted_ecoli_prod",
        "output_workspace" => $ws_name
    };
    lives_ok {
        $rast_ret = $annoutil->rast_genome($parms);
    } "annoutil->rast_genome call returns normally on genome $obj_Ecoli.";

    $parms = {
        "object_ref" => $obj_Echinacea,
        "output_genome_name" => "rasted_Echinace_prod",
        "output_workspace" => $ws_name
    };
    lives_ok {
        $rast_ret = $annoutil->rast_genome($parms);
    } "annoutil->rast_genome call returns normally on genome $obj_Echinacea.";

    $parms = {
        "object_ref" => $obj_asmb,
        "output_genome_name" => "rasted_assembly",
        "output_workspace" => $ws_name
    };

    lives_ok {
        $rast_ret = $annoutil->rast_genome($parms);
    } "annoutil->rast_genome call returns normally on assembly $obj_asmb.";
    if(defined($rast_ret) && defined($rast_ret->{output_genome_ref})) {
        ok (($rast_ret->{output_genome_ref} =~ m/[^\\w\\|._-]/), "rast_genome returns a VALID ref: $rast_ret->{output_genome_ref}");
    }
};
=cut

=begin
## testing Impl_rast_genome_assembly using $GEBA_1003_asmb from prod
 # to investigate the extra features
subtest 'Impl_rast_genome_assembly1' => sub {
    my $rast_ret = {};
    my $parms = {
        "object_ref" => $GEBA_1003_asmb,
        "output_genome_name" => "rasted_GEBA_1003_asmb",
        "output_workspace" => $ws_name,
        "create_report" => 1
    };

    lives_ok {
        $rast_ret = $rast_impl->rast_genome_assembly($parms);
    } "Impl rast_genome call on $GEBA_1003_asmb returns.";
    ok ( !defined($rast_ret->{output_genome_ref}),
         "rast_genome_assembly returned an undef object ref." );

    $parms = {
        "object_ref" => "63171/266/3",  # Ecoli_refseq_assembly
        "output_genome_name" => "rasted_Ecoli_refseq_assembly",
        "output_workspace" => $ws_name,
        "create_report" => 1
    };
    lives_ok {
        $rast_ret = $rast_impl->rast_genome_assembly($parms);
    } "Impl rast_genome call on $parms->{object_ref} returns.";
    ok ( !defined($rast_ret->{output_genome_ref}),
         "rast_genome_assembly returned an undef object ref." );

    $parms = {
        "object_ref" => $obj_Ecoli,  # Ecoli refseq genome
        "output_genome_name" => "rasted_Ecoli_refseq_genome",
        "output_workspace" => $ws_name,
        "create_report" => 1
    };
    lives_ok {
        $rast_ret = $rast_impl->rast_genome_assembly($parms);
        print "$parms->{object_ref} rast_genome_assembly returns:\n".Dumper($rast_ret);
    } "Impl rast_genome_assembly call on $parms->{object_ref} returns.";
    ok ($rast_ret->{output_genome_ref}, "success on $parms->{object_ref}");

    $parms = {
        "object_ref" => $obj_Echinacea,  # Echinacea genome
        "output_genome_name" => "rasted_Echinacea_genome",
        "output_workspace" => $ws_name,
        "create_report" => 1
    };
    lives_ok {
        $rast_ret = $rast_impl->rast_genome_assembly($parms);
        print "$parms->{object_ref} rast_genome_assembly returns:\n".Dumper($rast_ret);
    } "Impl rast_genome_assembly call on $parms->{object_ref} returns.";
    ok ($rast_ret->{output_genome_ref}, "success on $parms->{object_ref}");

    $parms = {
        "object_ref" => "63171/528/1",  # test_checkAll_Ecoli_Sept23
        "output_genome_name" => "rasted_test_chkAll",
        "output_workspace" => $ws_name,
        "create_report" => 1
    };
    lives_ok {
        $rast_ret = $rast_impl->rast_genome_assembly($parms);
        print "$parms->{object_ref} rast_genome_assembly returns:\n".Dumper($rast_ret);
    } "Impl rast_genome_assembly call on $parms->{object_ref} returns.";
    ok ( !defined($rast_ret->{output_genome_ref}),
         "rast_genome_assembly returned an undef object ref." );
};
=cut

=begin
## testing Impl_rast_genome_assembly using obj ids from prod ONLY
subtest 'Impl_rast_genome_assembly2' => sub {
    my $parms = {
        "object_ref" => $obj_Ecoli,
        "output_genome_name" => "rasted_genome",
        "output_workspace" => $ws_name,
        "create_report" => 1
    };
    my $rast_ret = {};
    lives_ok {
        $rast_ret = $rast_impl->rast_genome_assembly($parms);
    } "Impl rast_genome_assembly call returns normally on genome $obj_Ecoli";
    ok (keys %{ $rast_ret }, "rast_genome_assembly returns");
    ok ($rast_ret->{output_genome_ref} =~ m/[^\\w\\|._-]/,
        "rast_genome_assembly returns a valid ref");
    is ($rast_ret->{output_workspace}, $parms->{output_workspace}, "workspace is correct.");
    ok (exists($rast_ret->{report_ref}), "report created with report_ref");
    ok (exists($rast_ret->{report_name}), "report created with report_name");

    $parms = {
        "object_ref" => $obj_asmb,
        "output_genome_name" => "rasted_assembly",
        "output_workspace" => $ws_name
    };
    lives_ok {
        $rast_ret = $rast_impl->rast_genome_assembly($parms);
    } 'Impl rast_genome_assembly call returns.';
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

    # testing input text string with [\s\r\t\n|;]
    $parms->{input_genomes} = [];
    $parms->{input_assemblies} = [];
    $parms->{input_text} = '55141/242/1;55141/266/3;55141/212/1;|;63171/394/1';
    my $input_arr = ['55141/242/1', '55141/266/3', '55141/212/1', '63171/394/1'];
    lives_ok {
        $params = $annoutil->_get_bulk_rast_parameters($parms);
    } "anno_utils _get_bulk_rast_parameters call with input text containing '\\s\\r\\t\\n|;' returns normally.";
    my $param_objrefs = [];
    for my $param (@{ $params }) {
        push $param_objrefs, $param->{object_ref};
    }
    cmp_deeply( $param_objrefs, set(@$input_arr),
                "bulk input_text parameter has been correctly parsed.");
};


## testing Impl_rast_genomes_assemblies using obj ids from prod ONLY
subtest 'Impl_rast_genomes_assemblies' => sub {
    my $params = {
        "output_GenomeSet_name" => "out_genomeSet"
    };

    lives_ok {
        $params->{output_workspace} = $ws_name;
        $params->{input_genomes} = [$obj_Ecoli]; # array of prod objects
        $params->{input_assemblies} = [$obj_asmb]; # array of prod objects
        $params->{input_text} = '';
        my $ret_ann6 = $rast_impl->rast_genomes_assemblies($params);
    } "rast_impl rast_genomes_assemblies call on 1N1 returns normally.";

    lives_ok {
        $params->{output_workspace} = $ws_name;
        #$params->{input_genomes} = ["31020/5/1"]; # array of an appdev object
        $params->{input_genomes} = [$obj_Echinacea, $obj_Ecoli]; # array of prod objects
        $params->{input_assemblies} = [];
        $params->{input_text} = '';
        my $ret_ann7 = $rast_impl->rast_genomes_assemblies($params);
    } "rast_impl rast_genomes_assemblies call on array of 2 genomes returns normally.";

    lives_ok {
        $params->{output_workspace} = $ws_name;
        $params->{input_assemblies} = [$obj_asmb_refseq, $obj_asmb]; # array of prod objects
        $params->{input_genomes} = [];
        $params->{input_text} = '';
        my $ret_ann8 = $rast_impl->rast_genomes_assemblies($params);
    } "rast_impl rast_genomes_assemblies call on array of 2 assemblies returns normally.";

    lives_ok {
        $params->{output_workspace} = $ws_name;
        $params->{input_genomes} = [$obj_Echinacea, $obj_Ecoli]; # array of prod objects
        $params->{input_assemblies} = [$obj_asmb_refseq, $obj_asmb]; # array of prod objects
        #$params->{input_text} = '55141/266/3;55141/212/1;63171/394/1';
        $params->{input_text} = '36230/12/5;36230/13/5; 36230/14/5';
        my $ret_ann9 = $rast_impl->rast_genomes_assemblies($params);
    } "rast_impl rast_genomes_assemblies call on two arrays returns normally.";
};
=cut

=begin
subtest 'Impl_annotate_genome' => sub {
    my $assembly_obj_name = "Acidilobus_sp._CIS.fna";
    my $assembly_ref = RASTTestUtils::prepare_assembly($assembly_obj_name);
    my $genome_obj_name = 'Acidilobus_sp_CIS';

    my $parms={
        "input_contigset" => $assembly_obj_name,
        "workspace" => $ws_name,
        "output_genome" => 'Acidilobus_sp_7',
        "scientific_name" => 'Acidilobus sp 7',
        "domain" => 'A',
        "genetic_code" => '4',
        "call_features_CDS_prodigal" => '1',
    };
    my $rast_ann;
    throws_ok {
        $rast_ann = $rast_impl->annotate_genome($parms);
        my $genome_ref = $ws_name . "/" . $genome_obj_name;
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


# Test checking annotate_genomes input params for empty input_genomes and blank/undef genome_text
subtest 'annotation_genomes_throw_messages' => sub {
    my $error_message = qr/ERROR:Missing required inputs/;

    my $params = {
        "output_genome" => "out_genome_name",
        "workspace" => $ws_name
    };
    throws_ok {
        $params->{genome_text} = '';
        my $ret_ann1 = $rast_impl->annotate_genomes($params);
    } $error_message,
      'Blank genome_text plus undef input_genoms die correctly'
      or diag explain $params;

    $params = {
        "output_genome" => "out_genome_name",
        "workspace" => $ws_name
    };
    throws_ok {
        $params->{input_genomes} = [];
        my $ret_ann2 = $rast_impl->annotate_genomes($params);
    } $error_message,
      'Empty input_genomes plus undef genome_text die correctly'
      or diag explain $params;

    $params = {
        "output_genome" => "out_genome_name",
        "workspace" => $ws_name
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

## testing generate_genome_report on $GEBA_1003_asmb_ann
#  to check on extra/unwanted features
subtest 'generate_genome_report3' => sub {
    my $stats_ok = 'stats generation runs ok.\n';

    my $gff_path = $annoutil->_write_gff_from_genome($GEBA_1003_asmb_ann);

    my ($gff_contents, $attr_delimiter) = ([], '=');

    my (%gff_stats, %obj_stats);
    lives_ok {
        ($gff_contents, $attr_delimiter) = $annoutil->_parse_gff($gff_path, $attr_delimiter);
         %gff_stats = $annoutil->_generate_stats_from_gffContents($gff_contents);
         %obj_stats = $annoutil->_generate_stats_from_aa($GEBA_1003_asmb_ann);
    } $stats_ok;

    is(keys %gff_stats, 2, "_generate_stats_from_gffContents on $GEBA_1003_asmb_ann should return non-empty.\n");
    ok(exists($gff_stats{gene_role_map}), '_generate_stats_from_gffContents stats contains gene_roles.');
    ok(exists($gff_stats{function_roles}), '_generate_stats_from_gffContents stats contains function roles.');
    ok( $obj_stats{gc_content}, '_generate_stats_from_aa stats contains GC_Content.');
    ok( $obj_stats{contig_count}, '_generate_stats_from_aa stats contains contig_count.');
    ok( $obj_stats{num_features}, '_generate_stats_from_aa stats contains num_features.');

    my %ftr_tab = $annoutil->_get_feature_function_lookup($test_ftrs);
    my $test_msg = "Test message for reporting of annotation details";
    my $ftr_cnt = scalar @{$test_ftrs};

    my $ret_rpt = $annoutil->_generate_genome_report(
                $GEBA_1003_asmb_ann, $gff_contents, \%ftr_tab, $ftr_cnt, $test_msg);
    ok( exists($ret_rpt->{report_ref}), "Report generation on $GEBA_1003_asmb_ann returns report_ref.");
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

## Checking the rmrnas and cdss features of genome $GEBA_1003_asmb
#  and compare it with the rasted $GEBA_1003_asmb_ann.
subtest '_generate_stats_from_aa_GEBA' => sub {
    my ($gff_contents, $attr_delimiter) = ([], '=');
    my $gff_path = '';

    my $aa_stats_ok = 'stats generation from AA runs ok.\n';

    my %ret_stats = ();
    my %gn_stats = ();

    # $GEBA_1003_asmb
	my $gn_data = $annoutil->_fetch_object_data($GEBA_1003_asmb);
	$gn_stats{gc_content} = $gn_data->{gc_content};
	$gn_stats{contig_count} = $gn_data->{num_contigs};
	$gn_stats{feature_counts} = $gn_data->{feature_counts};
	if( $gn_data->{num_features} ) {
		$gn_stats{num_features} = $gn_data->{num_features};
	} else {
		$gn_stats{num_features} = $gn_data->{feature_counts}->{CDS};
	}

    lives_ok {
        %ret_stats = $annoutil->_generate_stats_from_aa($GEBA_1003_asmb);
    } $aa_stats_ok;
    ok(keys %ret_stats, "Statistics generated from annotated genome $GEBA_1003_asmb.");
    ok( $gn_stats{gc_content} == $ret_stats{gc_content},
        "GC Content of assembly matched in two different methods.");
    ok( $gn_stats{contig_count} == $ret_stats{contig_count},
        "Contig content of assembly matched in two different methods.");
    ok( !defined($gn_stats{num_features}) && !defined($ret_stats{num_features}),
        "num_features of assembly undefined.");
    ok( !defined($gn_stats{feature_counts}) && !defined($ret_stats{feature_counts}),
        "feature_counts of assembly undefined.");

    # $GEBA_1003_asmb_ann
    %ret_stats = ();
	$gn_data = $annoutil->_fetch_object_data($GEBA_1003_asmb_ann);
	$gn_stats{gc_content} = $gn_data->{gc_content};
	$gn_stats{contig_count} = $gn_data->{num_contigs};
	$gn_stats{feature_counts} = $gn_data->{feature_counts};
	if( $gn_data->{num_features} ) {
		$gn_stats{num_features} = $gn_data->{num_features};
	} else {
		$gn_stats{num_features} = $gn_data->{feature_counts}->{CDS};
	}
    lives_ok {
        %ret_stats = $annoutil->_generate_stats_from_aa($GEBA_1003_asmb_ann);
    } $aa_stats_ok;
    ok(keys %ret_stats, "Statistics generated from annotated genome $GEBA_1003_asmb_ann.");
    ok( $gn_stats{gc_content} == $ret_stats{gc_content},
        "GC Content of assembly matched in two different methods.");
    ok( $gn_stats{contig_count} == $ret_stats{contig_count},
        "Contig content of assembly matched in two different methods.");
    ok( $gn_stats{num_features} == $ret_stats{num_features},
        "num_features of annotated genome correct.");
    cmp_deeply( $gn_stats{feature_counts}, $ret_stats{feature_counts},
                "feature_counts of matched.");
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
=cut

=begin
#
## testing bulk_rast_genomes using obj ids from my own workspaces
## When GFU.save_one_genome failed to save, no rasted genome object(s) is created,
## so an empty hash is returned.
#
subtest 'bulk_rast_genomes' => sub {
    my $parms = {
          'input_genomes' => [
                               '55141/242/1',
                               '55141/212/1'
                             ],
          'create_report' => 0,
          'output_workspace' => $ws_name,
          'output_GenomeSet_name' => 'out_genomeSet',
          'input_assemblies' => [
                                  '55141/266/3',
                                  '55141/243/1'
                                ],
          'input_text' => '55141/266/3;55141/212/1;63171/394/1'
    };

    my $bulk_ann_ret;
    lives_ok {
        $bulk_ann_ret = $annoutil->bulk_rast_genomes($parms);
    } "annoutil->bulk_rast_genomes call on array of 2 genomes, array of 2 assemblies and input text string returns successfully.";
    #print "bulk_rast_genomes on input_genomes array returned:\n".Dumper($bulk_ann_ret);
    ok($bulk_ann_ret->{report_ref}, "Annotation report generated!!");
    ok($bulk_ann_ret->{output_genomeSet_ref}, "Annotated genomeSet saved!");
};


#
## testing bulk_rast_genomes using obj ids from public workspace id of 19217
#
subtest 'bulk_rast_genomes' => sub {
    my $params = {
        "input_assemblies" => [],
        "output_GenomeSet_name" => "bulk_genomeSet",
        "output_workspace" => $ws_name
    };
    my $refseq_gn1 = "19217/172902/1";
    my $refseq_gn2 = "19217/330276/1";
    my $refseq_gn3 = "19217/131931/1";  # failed to save in the narrative test
    my $refseq_gn4 = "19217/143319/1";  # failed to save in the narrative test
    my ($rfsq_ann1, $rfsq_ann2);

    lives_ok {
        $params->{input_genomes} = [$refseq_gn1, $refseq_gn2, $refseq_gn3, $refseq_gn4];
        $params->{input_text} = '';
        $rfsq_ann1 = $annoutil->bulk_rast_genomes($params);
    } "annoutil->bulk_rast_genomes call on array of 2 genomes returns.";
    #print "bulk_rast_genomes on input_genomes array returned:\n".Dumper($rfsq_ann1);
    ok($rfsq_ann1->{report_ref}, "Annotation report generated!!");
    ok($rfsq_ann1->{output_genomeSet_ref}, "Annotated genomeSet saved!");

    lives_ok {
        $params->{input_genomes} = [];
        $params->{input_text} = "19217/172902/1;19217/330276/1";
        $rfsq_ann2 = $annoutil->bulk_rast_genomes($params);
    } "annoutil->bulk_rast_genomes call on string of 2 genomes returns normally.";
    #print "bulk_rast_genomes on input_text string returned:\n".Dumper($rfsq_ann2);
    ok($rfsq_ann2->{report_ref}, "Annotation report generated!!");
    ok($rfsq_ann2->{output_genomeSet_ref}, "Annotated genomeSet saved!");
};
=cut

RASTTestUtils::clean_up();

done_testing();
