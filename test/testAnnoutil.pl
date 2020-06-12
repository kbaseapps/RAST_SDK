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
my $gff2 = 'data/metag_test/59111.assembled.gff';
my $fasta_scrt = 'fasta_file.fa';
my $gff_scrt = 'gff_file.gff';
my $prodigal_cmd = '/kb/runtime/bin/prodigal';

my $rast_impl = new RAST_SDK::RAST_SDKImpl();
my $annoutil = new anno_utils($config, $ctx);

my $scratch = $config->{'scratch'}; #'/kb/module/work/tmp';
my $rast_genome_dir = $annoutil->_create_rast_subdir($scratch, "genome_annotation_dir_");


=begin
sub genome_to_fasta {
    my($gn_ref) = @_;

    my $gfu = new installed_clients::GenomeFileUtilClient($call_back_url);

    my $fasta_result = $gfu->genome_proteins_to_fasta({"genome_ref" => $gn_ref});
    print "First 10 lines of the FASTA file from gfu->genome_proteins_to_fasta:\n";
    $annoutil->_print_fasta_gff(0, 10, $fasta_result->{file_path});
    return $fasta_result->{file_path};
}

sub generate_genome {
    my($ws, $gn_name, $gbff) = @_;
    my $gbff_path = catfile($rast_genome_dir, $gbff);

    copy($gbff, $gbff_path) || croak "Copy file failed: $!\n";

    my $gfu = new installed_clients::GenomeFileUtilClient($call_back_url);

    my $gn = $gfu->genbank_to_genome({
        "file" => {'path' => $gbff_path},
        "genome_name" => $gn_name,
        "workspace_name" => $ws,
        "generate_missing_genes" => 1
    });
    return $gn;
}
=cut

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
 'annotate_maedup_source',
 '1583389302.95394',
 '6c670c83-2a11-49ff-97bf-b1c3e2121f33'
 ]
 ],
 }];



=begin
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
        print "genome object returned on $obj_ref:\n".Dumper(keys %$ret_gn);
    } '_get_genome runs successfully';
    ok (@{$ret_gn->{features}} > 0, 'Genome has features!');
    is ($ret_gn->{assembly_ref}, '2901/78/1', 'found genome assembly ref';
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
=cut

=begin
#
## Global variables for the annotation process steps to share ##
#
my ($rast_ref, %rast_details1, %rast_details2);
my ($inputgenome1, $inputgenome2);
my ($parameters1, $parameters2);

# Test _set_parameters_by_input with genome/assembly object refs in prod
subtest '_set_parameters_by_input' => sub {
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
    } '_set_parameters_by_input runs successfully on genome';
    %rast_details1 = %{ $rast_ref }; # dereference
    $parameters1 = $rast_details1{parameters};

    ok (@{$inputgenome1->{features}} > 0, 'inputgenome has features.');
    cmp_deeply($expected_params1, $parameters1, 'parameters are correct');
    ok (@{$rast_details1{contigobj}{contigs}} == 1, 'inputgenome has 1 contig.');
    is ($inputgenome1->{assembly_ref}, '2901/78/1', 'inputgenome assembly_ref is correct.');

    # before merging with default gene call settings
    cmp_deeply($expected_params1, $parameters1, 'parameters has input param values.');

    #  merge with the default gene call settings
    my $default_params = $annoutil->_set_default_parameters();
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
    } '_set_parameters_by_input runs successfully on assembly';
    %rast_details2 = %{$rast_ref}; # dereference
    print "parameters on assembly:\n".Dumper($rast_details2{parameters});
    $parameters2 = $rast_details2{parameters};

    # before merging with default gene call settings
    cmp_deeply($expected_params3, $parameters2, 'parameters has input values');

    # merge with the default gene call settings
    $parameters2 = { %$default_params, %$parameters2 };
    my $expected_params4 = { %$default_params, %$expected_params3 };
    $rast_details2{parameters} = $parameters2;

    # after merging with default gene call settings
    cmp_deeply($expected_params4, $parameters2, 'parameters has default workflows');

    ok (@{$inputgenome2->{features}} == 0, 'inputgenome (assembly) has no features.');
    ok (@{$rast_details2{contigobj}{contigs}} == 1, 'inputgenome has 1 contig.');
    is ($inputgenome2->{assembly_ref}, $obj_asmb, 'inputgenome assembly_ref is correct.');
};


# Test _set_messageNcontigs with genome/assembly object refs in prod
subtest '_set_messageNcontigs' => sub {
    # a genome object
    lives_ok {
        ($rast_ref, $inputgenome1) = $annoutil->_set_messageNcontigs(
                                            \%rast_details1, $inputgenome1);
    } '_set_messageNcontigs runs successfully on genome';
    %rast_details1 = %{ $rast_ref }; # dereference
    $parameters1 = $rast_details1{parameters};
    my $msg1 = $rast_details1{message};
    my $tax1 = $rast_details1{tax_domain};

    ok (@{$inputgenome1->{features}} > 0, 'inputgenome has feature(s).');
    ok (@{$inputgenome1->{contigs}} > 0, 'inputgenome has contig(s).');
    ok (length($msg1) > 0, "Message for genome input has contents:\n$msg1");
    ok ($tax1 eq 'Bacteria', "tax_domain for genome input has value:$tax1");

    # an assembly object
    lives_ok {
        ($rast_ref, $inputgenome2) = $annoutil->_set_messageNcontigs(
                                            \%rast_details2, $inputgenome2);
    } '_set_messageNcontigs runs successfully on assembly';
    %rast_details2 = %{$rast_ref}; # dereference
    $parameters2 = $rast_details2{parameters};
    my $msg2 = $rast_details2{message};
    my $tax2 = $rast_details2{tax_domain};

    ok (@{$inputgenome2->{features}} == 0, 'inputgenome (assembly) has no features.');
    ok (@{$inputgenome2->{contigs}} > 0, 'inputgenome (assembly) has contig(s).');
    ok (length($msg2) > 0, "Message for assembly input has contents.");
    ok ($tax2 eq 'U', "tax_domain for assembly input has value:$tax2");

};


# Test _set_genecall_workflow with genome/assembly object refs in prod
subtest '_set_genecall_workflow' => sub {
    # a genome object
    lives_ok {
        $rast_ref = $annoutil->_set_genecall_workflow(
                                            \%rast_details1, $inputgenome1);
    } '_set_genecall_workflow runs successfully on genome';
    %rast_details1 = %{ $rast_ref }; # dereference
    my $genecall_workflow1 = $rast_details1{genecall_workflow};

    my $exp_gc_workflow1 = {
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
    cmp_deeply($exp_gc_workflow1, $genecall_workflow1, 'gc_workflow built correctly');

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
    } '_set_genecall_workflow runs successfully on assembly';
    %rast_details2 = %{$rast_ref}; # dereference
    my $msg2 = $rast_details2{message};
    my $genecall_workflow2 = $rast_details2{genecall_workflow};
    cmp_deeply($exp_gc_workflow2, $genecall_workflow2, 'gc_workflow built correctly');

};


# Test _set_annotation_workflow with genome/assembly object refs in prod
subtest '_set_annotation_workflow' => sub {
    # a genome object
    lives_ok {
        $rast_ref = $annoutil->_set_annotation_workflow(\%rast_details1);
    } '_set_genecall_workflow runs successfully on genome';
    %rast_details1 = %{ $rast_ref }; # dereference
    my $annomessage1 = $rast_details1{annomessage};
    my $annotate_workflow1 = $rast_details1{annotate_workflow};

    my $exp_ann_workflow1 = {
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
    print "genome annotation-msg:\n$annomessage1";
    cmp_deeply($exp_ann_workflow1, $annotate_workflow1, 'ann_workflow built correctly');

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
    } '_set_genecall_workflow runs successfully on assembly';
    %rast_details2 = %{$rast_ref}; # dereference
    my $annomessage2 = $rast_details2{annomessage};
    my $annotate_workflow2 = $rast_details2{annotate_workflow};
    print "assembly annotation-msg:\n$annomessage2";
    cmp_deeply($exp_ann_workflow2, $annotate_workflow2, 'ann_workflow built correctly');

};


# Test _renumber_features with genome/assembly object refs in prod
subtest '_renumber_features' => sub {
    # a genome object
    lives_ok {
        ($rast_ref, $inputgenome1) = $annoutil->_renumber_features(
                        \%rast_details1, $inputgenome1);
    } '_renumber_features runs successfully on genome';
    %rast_details1 = %{ $rast_ref }; # dereference
    my $msg1 = $rast_details2{message};
    print "genome merged-msg:\n$msg1";
    ok (@{$inputgenome1->{features}} > 0,
        "Genome inputgenome has ".scalar @{$inputgenome1->{features}}." features.");

    # an assembly object
    lives_ok {
        ($rast_ref, $inputgenome2) = $annoutil->_renumber_features(
                        \%rast_details2, $inputgenome2);
    } '_renumber_features runs successfully on assembly';
    %rast_details2 = %{$rast_ref}; # dereference
    my $msg2 = $rast_details2{message};
    print "genome merged-msg:\n$msg2";
    ok (@{$inputgenome2->{features}} == 0, "Assembly inputgenome has no features.");
};


# Test _pre_rast_call with genome/assembly object refs in prod
subtest '_pre_rast_call' => sub {
    # a genome object
    lives_ok {
        ($rast_ref, $inputgenome1) = $annoutil->_pre_rast_call(
                        \%rast_details1, $inputgenome1);
    } '_pre_rast_call runs successfully on genome';
    %rast_details1 = %{ $rast_ref }; # dereference
    my $genehash1 = $rast_details1{genehash};

    ok (keys %{$genehash1}, "Gene hash created from genome with elements.");
    if (defined($inputgenome1->{ontology_events})){
        ok (@{$inputgenome1->{ontology_events}} > 0,
            "There are ".scalar @{$inputgenome1->{ontology_events}}." ontology events for genome.");
    }
    ok (@{$inputgenome1->{features}} > 0, 'inputgenome has feature(s).');
    ok (@{$inputgenome1->{contigs}} > 0, 'inputgenome has contig(s).');

    # an assembly object
    lives_ok {
        ($rast_ref, $inputgenome2) = $annoutil->_pre_rast_call(
                        \%rast_details2, $inputgenome2);
    } '_renumber_features runs successfully on assembly';
    %rast_details2 = %{$rast_ref}; # dereference
    my $genehash2 = $rast_details2{genehash};

    ok (keys %{$genehash2} == 0, "Gene hash created from assembly with no elements.");
    if (defined($inputgenome2->{ontology_events})){
        ok (@{$inputgenome2->{ontology_events}} > 0,
            "There are ".scalar @{$inputgenome2->{ontology_events}}." ontology events for assembly.");
    }
    ok (@{$inputgenome2->{features}} == 0, 'inputgenome (assembly) has no features.');
    ok (@{$inputgenome2->{contigs}} > 0, 'inputgenome (assembly) has contig(s).');
};
=cut

=begin
#
## variables for testing _build_workflows, _run_rast_workflow
## and _pre_rast_call
my ($gc_wf_ret1, $gc_inputgenome1, $gc_wf1, %gc_details1,
    $gc_wf_ret2, $gc_inputgenome2, $gc_wf2, %gc_details2);

# Test _build_workflows with genome/assembly object refs in prod
subtest '_build_workflows' => sub {
    my $params1 = {
         output_genome_name => 'build_gcwf_gn_name',
         output_workspace => $ws,
         object_ref => $obj_Ecoli
    };

    lives_ok {
        ($gc_wf_ret1, $gc_inputgenome1) = $annoutil->_build_workflows($params1);
    } '_build_workflows returns normally';
    %gc_details1 = %{ $gc_wf_ret1 };
    $gc_wf1 = $gc_details1{genecall_workflow};
    ok (exists($gc_details1{genecall_workflow}), "generate genecall workflow:\n".Dumper($gc_wf1));
    ok (@{$gc_inputgenome1->{features}} > 0, 'inputgenome has some features.');
    ok (@{$gc_inputgenome1->{contigs}} > 0, 'inputgenome has contig(s).');

    my $params2 = {
         output_genome_name => 'build_gcwf_asmb_name',
         output_workspace => $ws,
         object_ref => $obj_asmb
    };
    lives_ok {
        ($gc_wf_ret2, $gc_inputgenome2) = $annoutil->_build_workflows($params2);
    } '_build_workflows returns normally';
    %gc_details2 = %{ $gc_wf_ret2 };
    $gc_wf2 = $gc_details2{genecall_workflow};
    ok (exists($gc_details2{genecall_workflow}), "generate genecall workflow:\n".Dumper($gc_wf2));
    ok (@{$gc_inputgenome2->{features}} == 0, 'inputgenome (assembly) has no features.');
    ok (@{$gc_inputgenome2->{contigs}} > 0, 'inputgenome (assembly) has contig(s).');
};


# Test _run_rast_workflow with genome/assembly object refs in prod
subtest '_run_rast_workflow' => sub {
    # a genome object
    my $rast_gn1;
    lives_ok {
        $rast_gn1 = $annoutil->_run_rast_workflow($gc_inputgenome1, $gc_wf1);
    } 'RAST run_pipeline call returns ERROR due to kmer data absence or other causes.';
    ok (@{$rast_gn1->{features}} > 0, "Returned genome has features.");
    cmp_deeply($rast_gn1, $gc_inputgenome1, 'rast workflow will not run locally');

    # an assembly object
    my $rast_gn2;
    lives_ok {
        $rast_gn2 = $annoutil->_run_rast_workflow($gc_inputgenome2, $gc_wf2);
    } 'RAST run_pipeline call returns ERROR due to kmer data absence or other causes.';
    ok (@{$rast_gn2->{features}} == 0, "Returned genome has NO features.");
    cmp_deeply($rast_gn2, $gc_inputgenome2, 'rast workflow will not run locally');
};


# Test _run_rast_genecalls with genome/assembly object refs in prod
subtest '_run_rast_genecalls' => sub {
    # a genome object
    my $params1 = {
         output_genome_name => 'gc_gn_name',
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
         output_genome_name => 'gc_asmb_name',
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
=cut

=begin
#
## Tesing only the annotation part of RAST
#
subtest '_run_rast_annotation' => sub {
    my $inputgenome = {
        features => []
    };
    foreach my $gene (sort keys %$protein_seqs){
        push(@{$inputgenome->{features}},{
            id => $gene,
            protein_translation => $protein_seqs->{$gene}
        });
    }

    throws_ok {
        my $rast_ret = $annoutil->_run_rast_annotation($inputgenome);
    } qr/ERROR calling rast run_pipeline/,
        'RAST run_pipeline call returns ERROR due to kmer data absence or other causes.';
};

# Test _annotate_process_allInOne with genome/assembly object refs in prod
subtest '_annotate_process_allInOne' => sub {
    my @regexes = ( qr/ERROR calling rast run_pipeline/,
                    qr/Can\'t call method/);
    my $allre= "(".join("|",@regexes).")";
    my $params = {
         output_genome_name => 'anno_gn_name',
         output_workspace => $ws,
         object_ref => $obj_Ecoli
    };
    my ($ain1_out, $out_msg);
    throws_ok {
        ($ain1_out, $out_msg) = $annoutil->_annotate_process_allInOne($params);
    } qr/$allre/,
        '_annotate_process_allInOne returns ERROR due to kmer data absence or other causes.';

    my $params2 = {
         output_genome_name => 'anno_asmb_name',
         output_workspace => $ws,
         object_ref => $obj_asmb
    };

    ## This should be fixed by using prodigal gene calling first...
    throws_ok {
        ($ain1_out, $out_msg) = $annoutil->_annotate_process_allInOne($params2);
    } qr/$allre/,
        '_annotate_process_allInOne returns ERROR due to kmer data absence or other causes.';
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


subtest '_parse_prodigal_results' => sub {
    my $prd_out_path = catfile($rast_genome_dir, 'prodigal_output.gff');
    my $trans_path = catfile($rast_genome_dir, 'protein_translation');
    copy($trans_file, $trans_path) || croak "Copy file failed: $!\n";

    # Prodigal generate a GFF output file
    my $out_type = 'gff';
    copy($ecoli_gff, $prd_out_path) || croak "Copy file failed: $!\n";
    my ($prd_results, %trans_tab) = $annoutil->_parse_prodigal_results(
                          $trans_path, $prd_out_path, $out_type);
    ok( @{$prd_results} >0, "Prodigal GFF parsing returns result.");
    ok( keys %trans_tab , "Prodigal GFF parsing returns translation table.");

    # Prodigal generate an SCO output file
    $out_type = 'sco';
    copy($ecoli_sco, $prd_out_path) || croak "Copy file failed: $!\n";
    ($prd_results, %trans_tab) = $annoutil->_parse_prodigal_results(
                          $trans_path, $prd_out_path, $out_type);
    ok( @{$prd_results} >0, "Prodigal SCO parsing returns result.");
    ok( keys %trans_tab , "Prodigal GFF parsing returns translation table.");
};

subtest '_prodigal_gene_call' => sub {
    my $p_input = $fasta1;
    my $md = 'meta';
    my $out_type = 'gff';
    my $gff_filename = catfile($rast_genome_dir, 'genome.gff');
    my $trans = catfile($rast_genome_dir, 'protein_translation');
    my $nuc = catfile($rast_genome_dir, 'nucleotide_seq');
    my $out_file = catfile($rast_genome_dir, 'prodigal_output').'.'.$out_type;

    my $prd_gene_results;
    lives_ok {
        ($out_file, $prd_gene_results) = $annoutil->_prodigal_gene_call(
                               $p_input, $trans, $nuc, $out_file, $out_type, $md);
    } 'Prodigal finished run 1.';
    ok( @{$prd_gene_results}=0, "Prodigal gene call on $p_input returns 0 result.");

    ## Check if the GFF file from Prodigal is tab delimited
    print "***********First 10 lines of prodigal gff file for $p_input:\n";
    $annoutil->_print_fasta_gff(0, 10, $out_file);

    $p_input = $ecoli_fasta;
    lives_ok {
        ($out_file, $prd_gene_results) = $annoutil->_prodigal_gene_call(
                               $p_input, $trans, $nuc, $out_file, $out_type, $md);
    } 'Prodigal finished run 2.';
    ok( @{$prd_gene_results}, "Prodigal gene call on $p_input returns result.");
    print "Prodigal gene call results2:\n".Dumper(@{$prd_gene_results}[0..10]);

    ## Check if the GFF file from Prodigal is tab delimited
    print "***********First 10 lines of prodigal gff file for $p_input:\n";
    $annoutil->_print_fasta_gff(0, 10, $out_file);

    $p_input = $fasta4;
    lives_ok {
        ($out_file, $prd_gene_results) = $annoutil->_prodigal_gene_call(
                               $p_input, $trans, $nuc, $out_file, $out_type, $md);
    } 'Prodigal finished run 3.';
    ok( @{$prd_gene_results}, "Prodigal gene call on $p_input returns result.");
    print "Prodigal gene call results3:\n".Dumper(@{$prd_gene_results}[0..10]);

    $p_input = $asmb_fasta;
    lives_ok {
        ($out_file, $prd_gene_results) = $annoutil->_prodigal_gene_call(
                               $p_input, $trans, $nuc, $out_file, $out_type, $md);
    } 'Prodigal finished run 4.';
    ok( @{$prd_gene_results}, "Prodigal gene call on $p_input returns result.");
    print "Prodigal gene call results3:\n".Dumper(@{$prd_gene_results}[0..10]);

    ## Check if the GFF file from Prodigal is tab delimited
    print "***********First 10 lines of prodigal gff file for $p_input:\n";
    $annoutil->_print_fasta_gff(0, 10, $out_file);
};

# test _glimmer3_gene_call
subtest '_glimmer3_gene_call' => sub {
    my $glimmer3_ok = "Glimmer3 gene call runs ok.";
    my $glimmer3_notOk = "ERROR";

    my $glimmer3_ret;
    throws_ok {
        $glimmer3_ret = $annoutil->_glimmer3_gene_call($fasta1);
    } qr/$glimmer3_notOk/,
        '_glimmer3_gene_call errors with contigs too short';

    lives_ok {
        $glimmer3_ret = $annoutil->_glimmer3_gene_call($ecoli_fasta);
    } $glimmer3_ok;

    lives_ok {
        $glimmer3_ret = $annoutil->_glimmer3_gene_call($fasta4);
    } $glimmer3_ok;
    ok( @{$glimmer3_ret} > 0, "_glimmer3_gene_call on $fasta4 returns gene call result.\n");
    print "Glimmer3 gene call results:\n". Dumper(@{$glimmer3_ret}[0..10]);
};

subtest '_prodigal_then_glimmer3' => sub {
    my $fa_input = $fasta4; # $ecoli_fasta;
    my $md = 'meta';
    my $out_type = 'gff';
    my $gff_filename = catfile($rast_genome_dir, 'genome.gff');
    my $trans = catfile($rast_genome_dir, 'protein_translation');
    my $nuc = catfile($rast_genome_dir, 'nucleotide_seq');
    my $out_file = catfile($rast_genome_dir, 'prodigal_output').'.'.$out_type;

    my ($pNg_gene_results, $pNg_gff_file);
    lives_ok {
        ($pNg_gff_file, $pNg_gene_results) = $annoutil->_prodigal_then_glimmer3(
                               $fa_input, $trans, $nuc, $out_file, $out_type, $md);
    } "_prodigal_then_glimmer3 finished run 1.";
    ok( @{$pNg_gene_results} > 0, "_prodigal_then_glimmer3 on $fa_input returns result.");
    print "_prodigal_then_glimmer3 on $fa_input results:\n".Dumper(@{$pNg_gene_results}[0..10]);

    ## Check if the GFF file from Prodigal is tab delimited
    # print "***********First 10 lines of prodigalNglimmer3 gff file for $fa_input:\n";
    # $annoutil->_print_fasta_gff(0, 10, $pNg_gff_file);

    $fa_input = $asmb_fasta;
    lives_ok {
        ($pNg_gff_file, $pNg_gene_results) = $annoutil->_prodigal_then_glimmer3(
                               $fa_input, $trans, $nuc, $out_file, $out_type, $md);
    } "_prodigal_then_glimmer3 finished run 2.";
    ok( @{$pNg_gene_results} > 0, "_prodigal_then_glimmer3 on $asmb_fasta returns result.");
    print "_prodigal_then_glimmer3 on $fa_input results:\n".Dumper(@{$pNg_gene_results}[0..10]);
};
=cut


=begin
## a CI object
my $ci_obj_id = '47032/4/8';

subtest 'annoutil_write_fasta_from_genome' => sub {
    # testing get the fasta from a genome using obj ids from prod ONLY
    my $fasta_fpath;
    lives_ok {
        $fasta_fpath = $annoutil->_write_fasta_from_genome($ci_obj_id);
    } 'Writing fasta from a genome runs ok';

    ok ((-s $fasta_fpath), "fasta file written for $ci_obj_id.\n");
    print "First 10 lines of the FASTA file:\n";
    $annoutil->_print_fasta_gff(0, 10, $fasta_fpath);
};
=cut

=begin
## testing rast-annotating genome functions
subtest '_write_fasta_from_genome' => sub {
    # testing get the fasta from a genome using obj ids from prod ONLY
    my $fasta_fpath;
    lives_ok {
        $fasta_fpath = $annoutil->_write_fasta_from_genome($obj_Ecoli);
    } 'Writing fasta from a genome runs ok';

    ok ((-s $fasta_fpath), "fasta file written for $obj_Ecoli.\n");
    print "First 10 lines of the FASTA file:\n";
    $annoutil->_print_fasta_gff(0, 10, $fasta_fpath);

    lives_ok {
        $fasta_fpath = $annoutil->_write_fasta_from_genome($obj_Echinacea);
    } 'Writing fasta from a genome runs ok';
    ok ((-s $fasta_fpath), "fasta file written for $obj_Echinacea.\n");
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


subtest '_save_genome_from_gff' => sub {
    ## repeat the portion from testing prodigal and
    ## _prodigal_then_glimmer3 in order to get the GFF
    my $md = 'meta';
    my $out_type = 'gff';
    my $trans = catfile($rast_genome_dir, 'protein_translation');
    my $nuc = catfile($rast_genome_dir, 'nucleotide_seq');
    my $out_file = catfile($rast_genome_dir, 'prodigal_output').'.'.$out_type;

    my $fa_input = $asmb_fasta;
    my $prd_gene_results;
    lives_ok {
        ($out_file, $prd_gene_results) = $annoutil->_prodigal_gene_call(
                               $fa_input, $trans, $nuc, $out_file, $out_type, $md);
    } 'Prodigal finished run.';
    ok( @{$prd_gene_results}, "Prodigal gene call on $fa_input returns result.");
    print "Prodigal gene call results3:\n".Dumper(@{$prd_gene_results}[0..10]);

    ## Check if the GFF file from Prodigal is tab delimited
    print "***********First 10 lines of prodigal gff file for $fa_input:\n";
    $annoutil->_print_fasta_gff(0, 10, $out_file);

    ## Test the _save_genome_from_gff function with Prodigal $out_file
    my $out_gn = 'prd_ed_Carsonella';
    my $input_asmb = $obj_asmb;
    my $mygn = {};
    lives_ok {
        $mygn = $annoutil->_save_genome_from_gff(
                   $ws, $out_gn, $input_asmb, $out_file);
    } "_save_genome run without errors on $input_asmb.\n";
    ok (exists $mygn->{genome_ref},
        "genome saved with genome_ref=$mygn->{genome_ref}");
    ok (exists $mygn->{genome_info}, 'genome saved with genome_info');
    is ($mygn->{genome_info}[1], $out_gn, 'saved genome name is correct');
    is ($mygn->{genome_info}[7], $ws, 'saved genome to the correct workspace');

    ## prodigal_then_glimmer3
    $fa_input = $asmb_fasta;
    my ($pNg_gff_file, $pNg_gene_results);
    my $out_file1 = catfile($rast_genome_dir, 'prodigal_output1').'.'.$out_type;
    lives_ok {
        ($pNg_gff_file, $pNg_gene_results) = $annoutil->_prodigal_then_glimmer3(
                               $fa_input, $trans, $nuc, $out_file1, $out_type, $md);
    } "_prodigal_then_glimmer3 finished run on $fa_input.";
    ok( @{$pNg_gene_results} > 0, "_prodigal_then_glimmer3 on $asmb_fasta returns result.");
    print "_prodigal_then_glimmer3 on $fa_input results:\n".Dumper(@{$pNg_gene_results}[0..10]);

    ## Check if the GFF file from Prodigal is tab delimited
    print "***********First 10 lines of prodigalNglimmer3 gff file for $fa_input:\n";
    $annoutil->_print_fasta_gff(0, 10, $pNg_gff_file);

    ## Test the _save_genome_from_gff function with prodigal_then_glimmer3 $gff_fpath
    $out_gn = 'pNg_Carsonella';
    $input_asmb = $obj_asmb;
    $mygn = {};
    lives_ok {
        $mygn = $annoutil->_save_genome_from_gff(
                    $ws, $out_gn, $input_asmb, $pNg_gff_file);
    } "_save_genome_from_gff run without errors on $input_asmb.\n";
    ok (exists $mygn->{genome_ref},
        "genome saved with genome_ref=$mygn->{genome_ref}");
    ok (exists $mygn->{genome_info}, 'genome saved with genome_info');
    is ($mygn->{genome_info}[1], $out_gn, 'saved genome name is correct');
    is ($mygn->{genome_info}[7], $ws, 'saved genome to the correct workspace');
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
        '_annotate_process_allInOne returns ERROR due to kmer data absence or other causes.';
    $parms = {
        "object_ref" => $obj_Echinacea,
        "output_genome_name" => "rasted_Echinace_prod",
        "output_workspace" => $ws
    };

    throws_ok {
        $rast_ret = $annoutil->rast_genome($parms);
    } qr/$allre/,
        '_annotate_process_allInOne returns ERROR due to kmer data absence or other causes.';
    $parms = {
        "object_ref" => $obj_asmb,
        "output_genome_name" => "rasted_assembly",
        "output_workspace" => $ws
    };

    lives_ok {
        $rast_ret = $annoutil->rast_genome($parms);
    } '$annoutil->rast_genome returns the original input because of empty features.';
    if(defined($rast_ret)) {
        ok (($rast_ret->{output_genome_ref} =~ m/[^\\w\\|._-]/), "rast_genome returns a VALID ref: $rast_ret->{output_genome_ref}");
    }
};


## testing Impl_rast_genome_assembly using obj ids from prod ONLY
subtest 'Impl_rast_genome_assembly' => sub {
    my $parms = {
        "object_ref" => $obj_Ecoli,
        "output_genome_name" => "rasted_genome",
        "output_workspace" => $ws
    };
    my $rast_ann;
    lives_ok {
        $rast_ann = $rast_impl->rast_genome_assembly($parms);
    } 'Impl rast_genome call returns normally.';

    $parms = {
        "object_ref" => $obj_asmb,
        "output_genome_name" => "rasted_assembly",
        "output_workspace" => $ws
    };
    lives_ok {
        $rast_ann = $rast_impl->rast_genome_assembly($parms);
    } 'Impl rast_genome call returns normally.';
    print "rast_genome_assembly returns:\n".Dumper($rast_ann);
};
=cut


=begin
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
    lives_ok {
        #$params->{input_genomes} = ["31020/5/1"]; # an appdev object
        $params->{input_genomes} = ["48109/9/1"]; # an prod object
        $params->{genome_text} = '';
        my $ret_ann4 = $rast_impl->annotate_genomes($params);
    } 'Should not throw error due to blank genome_text AND non-empty input_genoms';
};
=cut


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
    } "anno_utils rast_genomes_assemblies call on 1N1 returns normally.";

    lives_ok {
        $params->{output_workspace} = get_ws_name();
        #$params->{input_genomes} = ["31020/5/1"]; # array of an appdev object
        $params->{input_genomes} = [$obj_Echinacea, $obj_Ecoli]; # array of prod objects
        $params->{input_assemblies} = [];
        $params->{input_text} = '';
        my $ret_ann7 = $rast_impl->rast_genomes_assemblies($params);
    } "anno_utils rast_genomes_assemblies call on array of 2 genomes returns normally.";

    lives_ok {
        $params->{output_workspace} = get_ws_name();
        $params->{input_assemblies} = [$obj_asmb_refseq, $obj_asmb]; # array of prod objects
        $params->{input_genomes} = [];
        $params->{input_text} = '';
        my $ret_ann8 = $rast_impl->rast_genomes_assemblies($params);
    } "anno_utils rast_genomes_assemblies call on array of 2 assemblies returns normally.";

    lives_ok {
        $params->{output_workspace} = get_ws_name();
        $params->{input_genomes} = [$obj_Echinacea, $obj_Ecoli]; # array of prod objects
        $params->{input_assemblies} = [$obj_asmb_refseq, $obj_asmb]; # array of prod objects
        $params->{input_text} = '';
        my $ret_ann9 = $rast_impl->rast_genomes_assemblies($params);
    } "anno_utils rast_genomes_assemblies call on two arrays returns normally.";
};

=begin
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
=cut

=begin
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
        #print "Stats on $obj_Ecoli: \n".Dumper(\%ret_stats1);
    } $stats_ok;
    is(keys %ret_stats1, 2, "_generate_stats_from_gffContents on $obj_Ecoli should return non-empty.\n");

    my $gff_contents2 = [];
    my %ret_stats2;
    lives_ok {
        ($gff_contents2, $attr_delimiter) = $annoutil->_parse_gff($gff_path2, $attr_delimiter);
        %ret_stats2 = $annoutil->_generate_stats_from_gffContents($gff_contents2);
        #print "Stats on $obj7: \n".Dumper(\%ret_stats2);
    } $stats_ok;
    is(keys %ret_stats2, 2, "_generate_stats_from_gffContents on $obj_Ecoli_ann should return non-empty.\n");
    ok(exists($ret_stats2{gene_role_map}), '_generate_stats_from_gffContents stats contains gene_roles.');
    ok(exists($ret_stats2{function_roles}), '_generate_stats_from_gffContents stats contains function roles.');

    my %ftr_tab = $annoutil->_get_feature_function_lookup($test_ftrs);
    #print "\nFeature lookup:\n".Dumper(\%ftr_tab);
    my $ret_rpt = $annoutil->_generate_genome_report($obj_Ecoli, $obj_Ecoli_ann, $gff_contents1,
                                            $gff_contents2, \%ftr_tab);
    print "Report return: \n".Dumper($ret_rpt);
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
        #print "Stats on $obj_asmb_ann: \n".Dumper(\%ret_stats);
    } $stats_ok;
    is(keys %ret_stats, 2, "_generate_stats_from_gffContents on $obj_asmb_ann should return non-empty.\n");
    ok(exists($ret_stats{gene_role_map}), '_generate_stats_from_gffContents stats contains gene_roles.');
    ok(exists($ret_stats{function_roles}), '_generate_stats_from_gffContents stats contains function roles.');

    my %ftr_tab = $annoutil->_get_feature_function_lookup($test_ftrs);
    #print "\nFeature lookup:\n".Dumper(\%ftr_tab);
    my $ret_rpt = $annoutil->_generate_genome_report($obj_asmb, $obj_asmb_ann, [],
                                            $gff_contents, \%ftr_tab);
    print "Report return: \n".Dumper($ret_rpt);
    ok( exists($ret_rpt->{report_ref}), 'Report generation returns report_ref.');
    ok( exists($ret_rpt->{report_name}), 'Report generation returns report_name.');
    ok( exists($ret_rpt->{output_genome_ref}), 'Report generation returns output_gemome_ref.');
};
=cut

=begin
# testing _generate_stats_from_aa using obj ids from prod ONLY
subtest '_generate_stats_from_aa' => sub {
    my %ret_stats = $annoutil->_generate_stats_from_aa($obj4);
    #print "AMA stats return: \n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generation from AMA $obj4 returns result.");

    %ret_stats = $annoutil->_generate_stats_from_aa($obj5);
    #print "AMA stats return: \n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generation from AMA $obj5 returns result.");

    %ret_stats = $annoutil->_generate_stats_from_aa($obj7);
    #print "AMA stats return: \n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generation from AMA $obj7 returns result.");
};

# testing _generate_stats_from_gffContents using obj ids from prod ONLY
subtest '_generate_stats_from_gffContents' => sub {
    my $gff_path = $annoutil->_write_gff_from_ama($obj7);

    my ($gff_contents, $attr_delimiter) = ([], '=');
    ($gff_contents, $attr_delimiter) = $annoutil->_parse_gff($gff_path, $attr_delimiter);

    my %ret_stats = $annoutil->_generate_stats_from_gffContents($gff_contents);
    ok(keys %ret_stats, 'Statistics generation from gff_contents returns result.');
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
