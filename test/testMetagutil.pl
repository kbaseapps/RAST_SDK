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
use_ok "metag_utils";
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
my $mgutil = new metag_utils($config, $ctx);

my $scratch = $config->{'scratch'}; #'/kb/module/work/tmp';
my $rast_metag_dir = $mgutil->_create_rast_subdir($scratch, "metag_annotation_dir_");
my $rast_genome_dir = $mgutil->_create_rast_subdir($scratch, "genome_annotation_dir_");


sub genome_to_fasta {
    my($gn_ref) = @_;

    my $gfu = new installed_clients::GenomeFileUtilClient($call_back_url);

    my $fasta_result = $gfu->genome_proteins_to_fasta({"genome_ref" => $gn_ref});
    print "First 10 lines of the FASTA file from gfu->genome_proteins_to_fasta:\n";
    $mgutil->_print_fasta_gff(0, 10, $fasta_result->{file_path});
    return $fasta_result->{file_path};
}


sub generate_genome {
    my($ws, $gn_name, $gbff) = @_;
    my $gbff_path = catfile($rast_metag_dir, $gbff);

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

sub generate_metagenome {
    my($ws, $metag_name, $fasta, $gff) = @_;
    my $fasta_path = catfile($rast_metag_dir, $fasta_scrt);
    my $gff_path = catfile($rast_metag_dir, $gff_scrt);

    copy($fasta, $fasta_path) || croak "Copy file failed: $!\n";
    copy($gff, $gff_path) || croak "Copy file failed: $!\n";

    my $gfu = new installed_clients::GenomeFileUtilClient($call_back_url);
    my $mg = $gfu->fasta_gff_to_metagenome({
        "gff_file" => {'path' => $gff_path},
        "fasta_file" => {'path' => $fasta_path},
        "genome_name" => $metag_name,
        "workspace_name" => $ws,
        "generate_missing_genes" => 1
    });
    return $mg;
}

=begin
## global objects/variables for multiple subtests
my $ret_metag = generate_metagenome($ws, $out_name, $fasta1, $gff1);
print Dumper($ret_metag);
my $input_obj_ref = $ret_metag->{metagenome_ref};

my $input_fasta_file = catfile($rast_metag_dir, 'prodigal_input.fasta');
my $gff_filename = catfile($rast_metag_dir, 'genome.gff');
my ($fasta_contents, $gff_contents, $attr_delimiter) = ([], [], "=");

$input_fasta_file = $mgutil->_write_fasta_from_ama($input_obj_ref);
$gff_filename = $mgutil->_write_gff_from_ama($input_obj_ref);
# fetch protein sequences and gene IDs from fasta and gff files
$fasta_contents = $mgutil->_parse_fasta($input_fasta_file);
($gff_contents, $attr_delimiter) = $mgutil->_parse_gff(
                                       $gff_filename, $attr_delimiter);

my $gene_seqs = $mgutil->_extract_cds_sequences_from_fasta(
                    $fasta_contents, $gff_contents);
print "**********Gene sequences************\n" . Dumper($gene_seqs);

my $protein_seqs = $mgutil->_translate_gene_to_protein_sequences($gene_seqs);
print "**********Protein sequences************\n" . Dumper($protein_seqs);
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
my $obj4 = "55141/33/1";  # prod metag
my $obj5 = "55141/50/1";  # prod metag
my $obj6 = "55141/105/1";  # prod metag
my $obj7 = "55141/107/1";  # prod metag
my $obj8 = "55141/114/1";  # prod metag
my $obj9 = "55141/117/1";  # prod metag
my $obj10 = "55141/120/1";  # prod metag

my $asmb_fasta = $mgutil->_get_fasta_from_assembly($obj_asmb);

#####RAST_SDK Module test objects #####
my $obj2_1 = "63171/315/1";


#my $ecoli_fasta = genome_to_fasta($obj_Ecoli);
my $ecoli_fasta = $mgutil->_write_fasta_from_genome($obj_Ecoli);
unless (-s $ecoli_fasta) {
    print "**_write_fasta_from_genome on $obj_Ecoli ERROR: FASTA file is empty!\n";
}

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
subtest '_check_annotation_params' => sub {
    my $obj = '1234/56/7';

    my $missing_params = "Missing required parameters for annotating genome.\n";
    my $req1 = "'output_workspace' is required for running rast_genome.\n";
    my $req2 = "'object_ref' is required for running rast_genome.\n";

    throws_ok {
        $mgutil->_check_annotation_params()
    } qr/$missing_params/,
        '_check_annotation_params dies without params';

    throws_ok {
        $mgutil->_check_annotation_params( {} )
    } qr/$missing_params/,
        '_check_annotation_params dies with an empty hashref';

    throws_ok {
        my $p = {output_workspace => $ws,
                 output_genome_name => $outgn_name};
        print "input parameter=\n". Dumper($p);
        $mgutil->_check_annotation_params($p)
    } qr/$req2/,
        '_check_annotation_params dies with no object_ref';

    throws_ok {
        my $p = {output_workspace => $ws,
                 output_genome_name => $outgn_name,
                 object_ref => ''};
        print "input parameter=\n". Dumper($p);
        $mgutil->_check_annotation_params($p)
    } qr/$req2/,
        '_check_annotation_params dies with blank object_ref';

    throws_ok {
        $mgutil->_check_annotation_params(
            {object_ref => $obj,
             output_genome_name => $outgn_name})
    } qr/$req1/,
        '_check_annotation_params_metag dies with no outpout_workspace';

    throws_ok {
        $mgutil->_check_annotation_params(
            {workspace => $ws,
             output_genome_name => $outgn_name,
             obect_ref => $obj})
    } qr/$req1/,
        '_check_annotation_params dies with wrong workspace key';

    throws_ok {
        $mgutil->_check_annotation_params(
            {output_workspace => '',
             output_genome_name => $outgn_name,
             obect_ref => $obj})
    } qr/$req1/,
        '_check_annotation_params dies with blank workspace name';

    lives_ok {
        $mgutil->_check_annotation_params(
            {output_workspace => $ws,
             output_genome_name => $outgn_name,
             object_ref => 'abc/1/2'});
    } '_check_annotation_params object_ref check ok';

    lives_ok {
        $mgutil->_check_annotation_params(
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

    my $ret = $mgutil->_check_annotation_params(
            {output_workspace => $ws,
             output_genome_name => undef,
             object_ref => '456/1/2'});
    ok ($ret->{output_genome_name} eq $expected->{output_genome_name},
        'When undefined, '.$set_default_ok);

    $ret = $mgutil->_check_annotation_params(
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
        my $ret_parms1 = $mgutil->_check_bulk_annotation_params($params);
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
    my $ret_parms2 = $mgutil->_check_bulk_annotation_params($params);
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
    my $ret_parms3 = $mgutil->_check_bulk_annotation_params($params);
    ok (!@{$ret_parms3->{input_genomes}} && !@{$ret_parms3->{input_assemblies}},
        "No input genome or assembly was specified as input.");

    $params->{input_genomes} = ["48109/9/1"]; # array of a prod object
    $params->{input_text} = '';
    my $ret_parms4 = $mgutil->_check_bulk_annotation_params($params);
    ok ($ret_parms4->{input_genomes} eq $params->{input_genomes},
        "Input genome array is not empty.");

    $params->{input_genomes} = []; # array of a prod object
    $params->{input_text} = '48109/9/1;123/4/5';
    my $ret_parms5 = $mgutil->_check_bulk_annotation_params($params);
    ok ($ret_parms5->{input_text} eq $params->{input_text},
        "Input text is not empty.");
    ok (!@{$ret_parms5->{input_genomes}} && !@{$ret_parms5->{input_assemblies}},
        "No input genome or assembly was specified as input.");

    $params->{input_genomes} = "48109/9/1"; # non-array
    $params->{input_text} = '';
    my $ret_parms6 = $mgutil->_check_bulk_annotation_params($params);
    ok ($ret_parms6->{input_genomes}->[0] eq "48109/9/1",
        "Non array input genome converted into array");

    $params->{input_genomes} = [$obj_Ecoli]; # array of prod objects
    $params->{input_assemblies} = [$obj_asmb]; # array of prod objects
    $params->{input_text} = '';
    my $ret_parms7 = $mgutil->_check_bulk_annotation_params($params);
    ok (@{$ret_parms7->{input_genomes}} && @{$ret_parms7->{input_assemblies}},
        "Both input_genomes and input_assemblies arrays are not empty.");

    $params->{input_genomes} = [$obj_Echinacea, $obj_Ecoli]; # array of prod objects
    $params->{input_assemblies} = [];
    $params->{input_text} = '';
    my $ret_parms8 = $mgutil->_check_bulk_annotation_params($params);
    ok (@{$ret_parms8->{input_genomes}}==2 && @{$ret_parms8->{input_assemblies}}==0,
        "The input_genomes array has 2 elements while input_assemblies is empty.");

    $params->{input_assemblies} = [$obj_asmb_refseq, $obj_asmb]; # array of prod objects
    $params->{input_genomes} = [];
    $params->{input_text} = '';
    my $ret_parms9 = $mgutil->_check_bulk_annotation_params($params);
    ok (@{$ret_parms9->{input_genomes}}==0 && @{$ret_parms9->{input_assemblies}}==2,
        "The input_assemblies array has 2 elements while input_genomes is empty.");

    $params->{input_genomes} = [$obj_Echinacea, $obj_Ecoli]; # array of prod objects
    my $ret_parms10 = $mgutil->_check_bulk_annotation_params($params);
    ok (@{$ret_parms10->{input_genomes}}==2 && @{$ret_parms10->{input_assemblies}}==2,
        "Both input_genomes and input_assemblies arrays have 2 elements.");
};


subtest '_check_annotation_params_metag' => sub {
    my $obj = '1234/56/7';

    my $missing_params = "Missing required parameters for annotating metagenome.\n";
    my $req1 = "'output_workspace' is required for running rast_metagenome.\n";
    my $req2 = "'object_ref' is required for running rast_metagenome.\n";

    throws_ok {
        $mgutil->_check_annotation_params_metag()
    } qr/$missing_params/,
        '_check_annotation_params_metag dies without params';

    throws_ok {
        $mgutil->_check_annotation_params_metag( {} )
    } qr/$missing_params/,
        '_check_annotation_params_metag dies with an empty hashref';

    throws_ok {
        my $p = {output_workspace => $ws,
                 output_metagenome_name => $out_name};
        print "input parameter=\n". Dumper($p);
        $mgutil->_check_annotation_params_metag($p)
    } qr/$req2/,
        '_check_annotation_params_metag dies with no object_ref';

    throws_ok {
        my $p = {output_workspace => $ws,
                 output_metagenome_name => $out_name,
                 object_ref => ''};
        print "input parameter=\n". Dumper($p);
        $mgutil->_check_annotation_params_metag($p)
    } qr/$req2/,
        '_check_annotation_params_metag dies with blank object_ref';

    throws_ok {
        $mgutil->_check_annotation_params_metag(
            {object_ref => $obj,
             output_metagenome_name => $out_name})
    } qr/$req1/,
        '_check_annotation_params_metag dies with no outpout_workspace';

    throws_ok {
        $mgutil->_check_annotation_params_metag(
            {workspace => $ws,
             output_metagenome_name => $out_name,
             obect_ref => $obj})
    } qr/$req1/,
        '_check_annotation_params_metag dies with wrong workspace key';

    throws_ok {
        $mgutil->_check_annotation_params_metag(
            {output_workspace => '',
             output_metagenome_name => $out_name,
             obect_ref => $obj})
    } qr/$req1/,
        '_check_annotation_params_metag dies with blank workspace name';

    lives_ok {
        $mgutil->_check_annotation_params_metag(
            {output_workspace => $ws,
             output_metagenome_name => $out_name,
             object_ref => 'abc/1/2'});
    } '_check_annotation_params_metag object_ref check ok';

    lives_ok {
        $mgutil->_check_annotation_params_metag(
            {output_workspace => 'ab:c',
             output_metagenome_name => $out_name,
             object_ref => '456/1/2'});
    } '_check_annotation_params_metag workspace name check ok';

    # _check_annotation_params_metag passed
    my $expected = {
             output_workspace => $ws,
             output_metagenome_name => 'rast_annotated_metagenome',
             object_ref => '456/1/2'};
    my $set_default_ok = '_check_annotation_params_metag sets the default value for output_metagenome_name.';

    my $ret = $mgutil->_check_annotation_params_metag(
            {output_workspace => $ws,
             output_metagenome_name => undef,
             object_ref => '456/1/2'});
    ok ($ret->{output_metagenome_name} eq $expected->{output_metagenome_name},
        'When undefined, '.$set_default_ok);

    $ret = $mgutil->_check_annotation_params_metag(
            {output_workspace => $ws,
             output_metagenome_name => '',
             object_ref => '456/1/2'});
    ok ($ret->{output_metagenome_name} eq $expected->{output_metagenome_name},
        'When blank, '.$set_default_ok);
};

# test the html writing function using a small portion of a real AMA's stats data#
subtest '_write_html_from_stats' => sub {
    my %obj_stats = ('contig_count' => 123, 'id' => 'Test ID',
                     'num_features' => 456, 'gc_content' => 0.55);
    my %gff_stats = ('function_roles' => {
                                'FIG00500935: hypothetical protein' => {
                                                                         'gene_count' => 1,
                                                                         'gene_list' => '5785_1'
                                                                       },
                                'Respiratory nitrate reductase gamma chain (EC 1.7.99.4)' => {
                                                                                               'gene_count' => 3,
                                                                                               'gene_list' => '15370_3;8513_2;15_18'
                                                                                             },
                                'Twin-arginine translocation protein TatC' => {
                                                                                'gene_list' => '17470_2;13111_3;11624_3;11477_1;10224_3;7691_2;6222_4;4405_6;3274_2;2241_4;1632_1;718_3;575_5;256_4;249_3;15_58',
                                                                                'gene_count' => 16
                                                                              },
                                '1"-phosphate phophatase related protein' => {
                                                                        'gene_count' => 2,
                                                                        'gene_list' => '4599_1;627_2'
                                                                      }
                               },
            'gene_role_map' => {
                               '255_8' => 'Uptake hydrogenase small subunit precursor (EC 1.12.99.6)',
                               '13484_1' => 'Methyltransferase type 11',
                               '16059_2' => 'Methylated-DNA--protein-cysteine methyltransferase (EC 2.1.1.63)',
                               '11702_1' => 'Membrane-associated zinc metalloprotease',
                               '990_3' => 'Mobile element protein',
                               '11_44' => 'OsmC/Ohr family protein'
                               }
);

    my %subsys_info = $mgutil->_fetch_subsystem_info();
    my %ft_tab = $mgutil->_get_feature_function_lookup($test_ftrs);
    my @ret_html = $mgutil->_write_html_from_stats(\%obj_stats, \%gff_stats,
                                                   \%subsys_info, \%ft_tab, undef);
    ok(exists($ret_html[0]{path}), "html report written with file path returned.");
};

# test reading subsystem info
subtest '_fetch_subsystem_info' => sub {
    my $subsys_ok = 'subsystem info reading runs ok.\n';
    my %ret_subsysInfo = ();

    lives_ok {
        %ret_subsysInfo = $mgutil->_fetch_subsystem_info();
    } $subsys_ok;
    is(keys %ret_subsysInfo, 920, "_fetch_subsystem_info returns expected data.\n");
};

# test the prodigal functions
subtest '_build_prodigal_cmd' => sub {
    my $req = "An input FASTA/Genbank file is required for Prodigal to run.";
    my $set_default_ok = '_build_prodigal_cmd sets the default values ok.';
    my $outfile_default = 'prodigal_output.gff';
    my $outtype_default = 'gff';
    my $mode_default = 'meta';
    my $trans_default = 'protein_translation';
    my $nuc_default = 'nucleotide_seq';

    my $p_input = $fasta1;
    my ($v, $fpath, $f) = splitpath($p_input);

    my @exp_cmd_default = (
          $prodigal_cmd,
          '-i',
          $p_input,
          '-f',
          $outtype_default,
          '-o',
          catfile($fpath, $outfile_default),
          '-a',
          catfile($fpath, $trans_default),
          '-d',
          catfile($fpath, $nuc_default),
          '-c',
          '-g',
          11,
          '-q',
          '-p',
          $mode_default,
          '-m'
    );

    throws_ok {
        $mgutil->_build_prodigal_cmd()
    } qr/$req/,
        '_build_prodigal_cmd dies missing required fasta or gbk or gff file.';

    # set default values correctly
    my @ret_cmd = $mgutil->_build_prodigal_cmd($p_input);
    isa_ok( \@ret_cmd, 'ARRAY' );
    #print Dumper(\@ret_cmd);
    cmp_deeply(\@ret_cmd, \@exp_cmd_default, $set_default_ok);

    # with specified output file name and the rest default
    my $out_file = 'out_dir_path/your_outfile';
    @ret_cmd = $mgutil->_build_prodigal_cmd($p_input, '', '', $out_file, '', '');
    is($ret_cmd[0], $exp_cmd_default[0], 'Prodigal command set correctly');
    is($ret_cmd[2], $p_input, 'input file name set correctly');
    is($ret_cmd[6], $out_file, 'output file name set correctly');

    # with specified output type
    my $out_type = 'sco';
    @ret_cmd = $mgutil->_build_prodigal_cmd($p_input, '', '', $out_file, $out_type, '');
    is($ret_cmd[4], $out_type, 'output type set correctly');

    # with specified translation file name
    my $trans = 'out_dir_path/your_translation';
    @ret_cmd = $mgutil->_build_prodigal_cmd($p_input, $trans, '', $out_file, $out_type, '');
    is($ret_cmd[8], $trans, 'translation file set correctly');

    # with specified nucleotide sequence file name
    my $nuc = 'out_dir_path/your_nuc';
    @ret_cmd = $mgutil->_build_prodigal_cmd($p_input, $trans, $nuc, $out_file, $out_type, '');
    is($ret_cmd[10], $nuc, 'nucleotide sequence file set correctly');

    # with specified procedure mode
    my $md = 'single';
    @ret_cmd = $mgutil->_build_prodigal_cmd($p_input, $trans, $nuc, $out_file, $out_type, $md);
    is($ret_cmd[16], $md, 'nucleotide sequence file set correctly');

    # set default values for '' and undef inputs
    @ret_cmd = $mgutil->_build_prodigal_cmd($p_input, undef, '', undef, undef, '');
    isa_ok( \@ret_cmd, 'ARRAY' );
    cmp_deeply(\@ret_cmd, \@exp_cmd_default, $set_default_ok);
};

subtest '_run_prodigal' => sub {
    my $run_ok = '_run_prodigal_cmd runs ok.\n';
    my $run_err = 'ERROR Prodigal run failed';
    my $seq_too_short = 'Error:  Sequence must be 20000 characters';
    my $cannot_open = "Prodigal returned Error: can't open input file";
    my $prd_ret = 0;

    my $infile = '';
    my @p_cmd = (
          $prodigal_cmd,
          '-i',
          $infile
    );
    isnt($mgutil->_run_prodigal(@p_cmd), 0, $cannot_open);

    $infile = $fasta1;
    @p_cmd = (
          'wrong_path/prodigal',
          '-i',
          $infile
    );
    isnt($mgutil->_run_prodigal(@p_cmd), 0, "Can't exec wrong prodigal cmd.\n");

    @p_cmd = (
          $prodigal_cmd,
          '-i',
          $infile
    );
    lives_ok {
        $prd_ret = $mgutil->_run_prodigal(@p_cmd)
    } $run_ok;
    isnt($prd_ret, 0, "Prodigal returned:$seq_too_short\n");

    @p_cmd = (
        $prodigal_cmd,
	'-i',
	$infile,
	'-p',
	'meta'
    );
    lives_ok {
        $prd_ret = $mgutil->_run_prodigal(@p_cmd)
    } $run_ok;
    is($prd_ret, 0, "Prodigal runs ok with -p meta option.\n");

    ##---------A long file takes much longer time!!!!!---------##
    $infile = $fasta2;
    @p_cmd = (
          $prodigal_cmd,
          '-i',
          $infile
    );
    lives_ok {
        $prd_ret = $mgutil->_run_prodigal(@p_cmd)
    } $run_ok;
    is($prd_ret, 0, "_run_prodigal successfully returned.\n");
};


subtest '_parse_translation' => sub {
    my $trans_path = catfile($rast_metag_dir, 'trans_scrt');
    copy($trans_file, $trans_path) || croak "Copy file failed: $!\n";

    %trans_tab = $mgutil->_parse_translation($trans_path);
    ok( keys %trans_tab , "Prodigal translation parsing returns result.");
};

subtest '_parse_sco' => sub {
    $sco_tab = $mgutil->_parse_sco($ecoli_sco, %trans_tab);
    ok( @{$sco_tab} >0, "Prodigal SCO parsing returns result.");
};

subtest '_parseNwrite_gff' => sub {
    # testing get the gff from a genome using obj ids from prod ONLY
    my $gff_fpath;
    lives_ok {
        $gff_fpath = $mgutil->_write_gff_from_genome($obj_Echinacea);
    } 'Writing gff from a genome runs ok';

    ok((-e $gff_fpath), "GFF file created for $obj_Echinacea.\n");
    ok ((-s $gff_fpath), "GFF file written for $obj_Echinacea to $gff_fpath.\n");

    print "First 20 lines of the GFF file written from genome $obj_Echinacea:\n";
    $mgutil->_print_fasta_gff(0, 20, $gff_fpath);

    my $temp_gff_file = catfile($tmp_write_dir, 'test_gff_file.gff');
    copy $gff_fpath, $temp_gff_file;
    print "First 20 lines of the GFF file copied from genome $gff_fpath:\n";
    $mgutil->_print_fasta_gff(0, 20, $temp_gff_file);

    my $gff_contents1;
    lives_ok {
        ($gff_contents1, $attr_delimiter) = $mgutil->_parse_gff($gff_fpath, '=');
    } "Testing _parse_gff on $gff_fpath succeeded.";
    ok( @{$gff_contents1} >0, "Parsing GFF on $gff_fpath returns result.\n");

    my $gff_contents2;
    lives_ok {
        ($gff_contents2, $attr_delimiter) = $mgutil->_parse_gff($temp_gff_file, '=');
    } "Testing _parse_gff on $temp_gff_file succeeded.";
    ok( @{$gff_contents1} >0, "Parsing GFF on $temp_gff_file returns result.\n");

    is_deeply($gff_contents1, $gff_contents2, 'GFF data structures should be the same!');

    my $Echinacea_gff_file = "data/Echinacea_purpurea_1762.gff";
    my $gff_contents3;
    lives_ok {
        ($gff_contents3, $attr_delimiter) = $mgutil->_parse_gff($Echinacea_gff_file, '=');
    } "Testing _parse_gff on file $Echinacea_gff_file succeeded.";
    ok( @{$gff_contents3} >0, "Parsing GFF on $Echinacea_gff_file returns result.\n");
    print "Parsed ". scalar @{$gff_contents3}." GFF contents.\n";

    is_deeply($gff_contents1, $gff_contents3, 'GFF data structures should be the same 2!');

    my $test_gff_file_written;
    lives_ok {
        $test_gff_file_written = catfile($rast_genome_dir, 'test_written.gff');
        $mgutil->_write_gff($gff_contents1, $test_gff_file_written , '=');
    } "Writing the gff contents back to a gff file is ok";
    ok ( (-s $test_gff_file_written), "GFF file written to $test_gff_file_written.\n");
};


subtest '_parse_prodigal_results' => sub {
    my $prd_out_path = catfile($rast_metag_dir, 'prodigal_output.gff');
    my $trans_path = catfile($rast_metag_dir, 'protein_translation');
    copy($trans_file, $trans_path) || croak "Copy file failed: $!\n";

    # Prodigal generate a GFF output file
    my $out_type = 'gff';
    copy($ecoli_gff, $prd_out_path) || croak "Copy file failed: $!\n";
    my ($prd_results, %trans_tab) = $mgutil->_parse_prodigal_results(
                          $trans_path, $prd_out_path, $out_type);
    ok( @{$prd_results} >0, "Prodigal GFF parsing returns result.");
    ok( keys %trans_tab , "Prodigal GFF parsing returns translation table.");

    # Prodigal generate an SCO output file
    $out_type = 'sco';
    copy($ecoli_sco, $prd_out_path) || croak "Copy file failed: $!\n";
    ($prd_results, %trans_tab) = $mgutil->_parse_prodigal_results(
                          $trans_path, $prd_out_path, $out_type);
    ok( @{$prd_results} >0, "Prodigal SCO parsing returns result.");
    ok( keys %trans_tab , "Prodigal GFF parsing returns translation table.");
};

subtest '_prodigal_gene_call' => sub {
    my $p_input = $fasta1;
    my $md = 'meta';
    my $out_type = 'gff';
    my $gff_filename = catfile($rast_metag_dir, 'genome.gff');
    my $trans = catfile($rast_metag_dir, 'protein_translation');
    my $nuc = catfile($rast_metag_dir, 'nucleotide_seq');
    my $out_file = catfile($rast_metag_dir, 'prodigal_output').'.'.$out_type;

    my $prd_gene_results;
    lives_ok {
        ($out_file, $prd_gene_results) = $mgutil->_prodigal_gene_call(
                               $p_input, $trans, $nuc, $out_file, $out_type, $md);
    } 'Prodigal finished run 1.';
    ok( @{$prd_gene_results}=0, "Prodigal gene call on $p_input returns 0 result.");

    ## Check if the GFF file from Prodigal is tab delimited
    print "***********First 10 lines of prodigal gff file for $p_input:\n";
    $mgutil->_print_fasta_gff(0, 10, $out_file);

    $p_input = $ecoli_fasta;
    lives_ok {
        ($out_file, $prd_gene_results) = $mgutil->_prodigal_gene_call(
                               $p_input, $trans, $nuc, $out_file, $out_type, $md);
    } 'Prodigal finished run 2.';
    ok( @{$prd_gene_results}, "Prodigal gene call on $p_input returns result.");
    print "Prodigal gene call results2:\n".Dumper(@{$prd_gene_results}[0..10]);

    ## Check if the GFF file from Prodigal is tab delimited
    print "***********First 10 lines of prodigal gff file for $p_input:\n";
    $mgutil->_print_fasta_gff(0, 10, $out_file);

    $p_input = $fasta4;
    lives_ok {
        ($out_file, $prd_gene_results) = $mgutil->_prodigal_gene_call(
                               $p_input, $trans, $nuc, $out_file, $out_type, $md);
    } 'Prodigal finished run 3.';
    ok( @{$prd_gene_results}, "Prodigal gene call on $p_input returns result.");
    print "Prodigal gene call results3:\n".Dumper(@{$prd_gene_results}[0..10]);

    $p_input = $asmb_fasta;
    lives_ok {
        ($out_file, $prd_gene_results) = $mgutil->_prodigal_gene_call(
                               $p_input, $trans, $nuc, $out_file, $out_type, $md);
    } 'Prodigal finished run 4.';
    ok( @{$prd_gene_results}, "Prodigal gene call on $p_input returns result.");
    print "Prodigal gene call results3:\n".Dumper(@{$prd_gene_results}[0..10]);

    ## Check if the GFF file from Prodigal is tab delimited
    print "***********First 10 lines of prodigal gff file for $p_input:\n";
    $mgutil->_print_fasta_gff(0, 10, $out_file);
};

# test _glimmer3_gene_call
subtest '_glimmer3_gene_call' => sub {
    my $glimmer3_ok = "Glimmer3 gene call runs ok.";
    my $glimmer3_notOk = "ERROR";

    my $glimmer3_ret;
    throws_ok {
        $glimmer3_ret = $mgutil->_glimmer3_gene_call($fasta1);
    } qr/$glimmer3_notOk/,
        '_glimmer3_gene_call errors with contigs too short';

    lives_ok {
        $glimmer3_ret = $mgutil->_glimmer3_gene_call($ecoli_fasta);
    } $glimmer3_ok;

    lives_ok {
        $glimmer3_ret = $mgutil->_glimmer3_gene_call($fasta4);
    } $glimmer3_ok;
    ok( @{$glimmer3_ret} > 0, "_glimmer3_gene_call on $fasta4 returns gene call result.\n");
    print "Glimmer3 gene call results:\n". Dumper(@{$glimmer3_ret}[0..10]);
};

subtest '_prodigal_then_glimmer3' => sub {
    my $fa_input = $fasta4; # $ecoli_fasta;
    my $md = 'meta';
    my $out_type = 'gff';
    my $gff_filename = catfile($rast_metag_dir, 'genome.gff');
    my $trans = catfile($rast_metag_dir, 'protein_translation');
    my $nuc = catfile($rast_metag_dir, 'nucleotide_seq');
    my $out_file = catfile($rast_metag_dir, 'prodigal_output').'.'.$out_type;

    my ($pNg_gene_results, $pNg_gff_file);
    lives_ok {
        ($pNg_gff_file, $pNg_gene_results) = $mgutil->_prodigal_then_glimmer3(
                               $fa_input, $trans, $nuc, $out_file, $out_type, $md);
    } "_prodigal_then_glimmer3 finished run 1.";
    ok( @{$pNg_gene_results} > 0, "_prodigal_then_glimmer3 on $fa_input returns result.");
    print "_prodigal_then_glimmer3 on $fa_input results:\n".Dumper(@{$pNg_gene_results}[0..10]);

    ## Check if the GFF file from Prodigal is tab delimited
    # print "***********First 10 lines of prodigalNglimmer3 gff file for $fa_input:\n";
    # $mgutil->_print_fasta_gff(0, 10, $pNg_gff_file);

    $fa_input = $asmb_fasta;
    lives_ok {
        ($pNg_gff_file, $pNg_gene_results) = $mgutil->_prodigal_then_glimmer3(
                               $fa_input, $trans, $nuc, $out_file, $out_type, $md);
    } "_prodigal_then_glimmer3 finished run 2.";
    ok( @{$pNg_gene_results} > 0, "_prodigal_then_glimmer3 on $asmb_fasta returns result.");
    print "_prodigal_then_glimmer3 on $fa_input results:\n".Dumper(@{$pNg_gene_results}[0..10]);
};


subtest '_write_fasta_from_ama' => sub {
    my $fa_test1 = $mgutil->_write_fasta_from_ama($input_obj_ref);
    ok((-e $fa_test1), 'fasta file created');
    ok((-s $fa_test1), 'fasta file has data');
    # ok(compare($fa_test1, $fasta1) ==dd 0, 'fasta file written correctly');
};

subtest '_write_gff_from_ama' => sub {
    # test by using prod obj of type KBaseGenomeAnnotations.Assembly-5.0
    my $obj_wrong_type = "55141/119/1";
    throws_ok {
       $gff_test2 = $mgutil->_write_gff_from_ama($obj_wrong_type);
    } qr/ValueError/,
        '_write_gff_from_ama dies due to wrong object type.';

    # test using obj in the current test workspace
    my $gff_test1 = $mgutil->_write_gff_from_ama($input_obj_ref);
    ok((-e $gff_test1), 'gff file created');
    ok((-s $gff_test1), 'gff file has data');
};

subtest '_save_metagenome' => sub {
    my $req_params = "Missing required parameters for saving metagenome.\n";
    my $not_found = "file not found.\n";
    my $req1 = "Both 'output_workspace' and 'output_metagenome_name' are required.\n";
    my $req2 = "Both 'obj_ref' and 'gff_file' are required.\n";

    throws_ok {
        $mgutil->_save_metagenome()
    } qr/$req_params/,
        '_save_metagenome dies without params';

    throws_ok {
        $mgutil->_save_metagenome($ws, $out_name, undef, 'abc')
    } qr/$req_params/,
        '_save_metagenome dies with no fasta file';

    throws_ok {
        $mgutil->_save_metagenome($ws, $out_name, undef)
    } qr/$req_params/,
        '_save_metagenome dies with no gff file';

    throws_ok {
        $mgutil->_save_metagenome($ws, $out_name, 'a/b/c', 'def')
    } qr/$not_found/,
        '_save_metagenome dies because input file not found.';

    throws_ok {
        $mgutil->_save_metagenome($ws, $out_name, 'not_readable/4/5', $gff2)
    } qr/cannot be accessed/,
      '_save_metagenome dies because given workspace "not_readable" cannot be read';

    throws_ok {
        $mgutil->_save_metagenome($ws, $out_name, $fasta2,
                                  'data/nosuchfile.gff')
    } qr/$not_found/,
      '_save_metagenome dies because GFF file not found';

    my $gff_path = catfile($rast_metag_dir, $gff_scrt);
    copy($gff1, $gff_path) || croak "Copy file failed: $!\n";
    my $mymetag = {};
    lives_ok {
        $mymetag = $mgutil->_save_metagenome(
                       $ws, $out_name, $input_obj_ref, $gff_path)
    } '_save_metagenome run without errors on short_one.\n';
    ok (exists $mymetag->{metagenome_ref},
        "metagenome saved with metagenome_ref='$mymetag->{metagenome_ref}'");
    ok (exists $mymetag->{metagenome_info}, 'metagenome saved with metagenome_info');
    is ($mymetag->{metagenome_info}[1], $out_name, 'saved metagenome name is correct');
    is ($mymetag->{metagenome_info}[7], $ws, 'saved metagenome to the correct workspace');
};

# test by using prod/appdev obj id
subtest 'annotate_metagenome_prod' => sub {
    my $parms = {
        #"object_ref" => $obj1, # appdev obj
        "object_ref" => $obj10, # prod obj
        "output_metagenome_name" => "rasted_AMA",
        "output_workspace" => $ws
    };
    throws_ok {
        my $rast_ann = $rast_impl->annotate_metagenome($parms);
    } qr/ERROR calling rast run_pipeline/,
        'RAST annotate_metagenome call returns ERROR due to kmer data absence or other causes.';

};

subtest '_run_rast_genecalls' => sub {
    my $input_obj = $obj_Ecoli;
    my $inparams = {
        "object_ref" => $input_obj,
        "output_genome_name" => "ann_gn",
        "output_workspace" => $ws,
        "create_report" => 0
    };
    throws_ok {
        my $rast_ret = $mgutil->_run_rast_genecalls($inparams);
    } qr/ERROR calling rast run_pipeline/,
      '_run_rast_genecalls threw ERROR when local testing.';
};

=cut

=begin
## a CI object
my $ci_obj_id = '47032/4/8';

subtest 'mgutil_write_fasta_from_genome' => sub {
    # testing get the fasta from a genome using obj ids from prod ONLY
    my $fasta_fpath;
    lives_ok {
        $fasta_fpath = $mgutil->_write_fasta_from_genome($ci_obj_id);
    } 'Writing fasta from a genome runs ok';

    ok ((-s $fasta_fpath), "fasta file written for $ci_obj_id.\n");
    print "First 10 lines of the FASTA file:\n";
    $mgutil->_print_fasta_gff(0, 10, $fasta_fpath);
};


## testing rast-annotating genome functions
subtest '_write_fasta_from_genome' => sub {
    # testing get the fasta from a genome using obj ids from prod ONLY
    my $fasta_fpath;
    lives_ok {
        $fasta_fpath = $mgutil->_write_fasta_from_genome($obj_Ecoli);
    } 'Writing fasta from a genome runs ok';

    ok ((-s $fasta_fpath), "fasta file written for $obj_Ecoli.\n");
    print "First 10 lines of the FASTA file:\n";
    $mgutil->_print_fasta_gff(0, 10, $fasta_fpath);

    lives_ok {
        $fasta_fpath = $mgutil->_write_fasta_from_genome($obj_Echinacea);
    } 'Writing fasta from a genome runs ok';
    ok ((-s $fasta_fpath), "fasta file written for $obj_Echinacea.\n");
};

subtest '_write_gff_from_genome' => sub {
    # testing get the gff from a genome using obj ids from prod ONLY
    my $gff_fpath;
    lives_ok {
        $gff_fpath = $mgutil->_write_gff_from_genome($obj_Ecoli);
    } 'Writing gff from a genome runs ok';

    ok((-e $gff_fpath), "GFF file created for $obj_Ecoli.\n");
    ok ((-s $gff_fpath), "GFF file written for $obj_Ecoli.\n");

    print "First 10 lines of the GFF file:\n";
    $mgutil->_print_fasta_gff(0, 10, $gff_fpath);

    lives_ok {
        $gff_fpath = $mgutil->_write_gff_from_genome($obj_Echinacea);
    } 'Writing gff from a genome runs ok';

    ok((-e $gff_fpath), "GFF file created for $obj_Echinacea.\n");
    ok ((-s $gff_fpath), "GFF file written for $obj_Echinacea.\n");

    print "ALL lines of the GFF file:\n";
    $mgutil->_print_fasta_gff(0, 2000, $gff_fpath);
};


subtest '_save_genome' => sub {
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
        ($out_file, $prd_gene_results) = $mgutil->_prodigal_gene_call(
                               $fa_input, $trans, $nuc, $out_file, $out_type, $md);
    } 'Prodigal finished run.';
    ok( @{$prd_gene_results}, "Prodigal gene call on $fa_input returns result.");
    print "Prodigal gene call results3:\n".Dumper(@{$prd_gene_results}[0..10]);

    ## Check if the GFF file from Prodigal is tab delimited
    print "***********First 10 lines of prodigal gff file for $fa_input:\n";
    $mgutil->_print_fasta_gff(0, 10, $out_file);

    ## Test the _save_genome function with Prodigal $out_file
    my $out_gn = 'prd_ed_Carsonella';
    my $input_asmb = $obj_asmb;
    my $mygn = {};
    lives_ok {
        $mygn = $mgutil->_save_genome($ws, $out_gn, $input_asmb, $out_file);
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
        ($pNg_gff_file, $pNg_gene_results) = $mgutil->_prodigal_then_glimmer3(
                               $fa_input, $trans, $nuc, $out_file1, $out_type, $md);
    } "_prodigal_then_glimmer3 finished run on $fa_input.";
    ok( @{$pNg_gene_results} > 0, "_prodigal_then_glimmer3 on $asmb_fasta returns result.");
    print "_prodigal_then_glimmer3 on $fa_input results:\n".Dumper(@{$pNg_gene_results}[0..10]);

    ## Check if the GFF file from Prodigal is tab delimited
    print "***********First 10 lines of prodigalNglimmer3 gff file for $fa_input:\n";
    $mgutil->_print_fasta_gff(0, 10, $pNg_gff_file);
    ## Test the _save_genome function with prodigal_then_glimmer3 $gff_fpath
    $out_gn = 'pNg_Carsonella';
    $input_asmb = $obj_asmb;
    $mygn = {};
    lives_ok {
        $mygn = $mgutil->_save_genome($ws, $out_gn, $input_asmb, $pNg_gff_file);
    } "_save_genome run without errors on $input_asmb.\n";
    ok (exists $mygn->{genome_ref},
        "genome saved with genome_ref=$mygn->{genome_ref}");
    ok (exists $mygn->{genome_info}, 'genome saved with genome_info');
    is ($mygn->{genome_info}[1], $out_gn, 'saved genome name is correct');
    is ($mygn->{genome_info}[7], $ws, 'saved genome to the correct workspace');
};


subtest 'mgutil_rast_genome' => sub {
    # testing metag_utils rast_genome using obj ids from prod ONLY
    my $parms = {
        "object_ref" => $obj_Ecoli,
        "output_genome_name" => "rasted_ecoli_prod",
        "output_workspace" => $ws
    };
    my $rast_ref;
    throws_ok {
        $rast_ref = $mgutil->rast_genome($parms);
    } qr/ERROR calling rast run_pipeline/,
        'metag_utils rast_genome call returns ERROR due to kmer data absence or other causes.';
    if(defined($rast_ref)) {
        print "rast_genome returns: $rast_ref";
        ok (($rast_ref !~ m/[^\\w\\|._-]/), "rast_genome returned an INVALID ref: $rast_ref");
    }

    $parms = {
        "object_ref" => $obj_Echinacea,
        "output_genome_name" => "rasted_Echinace_prod",
        "output_workspace" => $ws
    };

    throws_ok {
        $rast_ref = $mgutil->rast_genome($parms);
    } qr/ERROR calling rast run_pipeline/,
        'metag_utils rast_genome call returns ERROR due to kmer data absence or other causes.';
    if(defined($rast_ref)) {
        print "rast_genome returns: $rast_ref";
        ok (($rast_ref !~ m/[^\\w\\|._-]/), "rast_genome returns an INVALID ref: $rast_ref");
    }

    $parms = {
        "object_ref" => $obj_asmb,
        "output_genome_name" => "rasted_assembly",
        "output_workspace" => $ws,
        "create_report" => 1
    };
    throws_ok {
        $rast_ref = $mgutil->rast_genome($parms);
    } qr/ERROR calling rast run_pipeline/,
        'metag_utils rast_genome call returns ERROR due to kmer data absence or other causes.';
};


subtest 'Impl_rast_genome_assembly' => sub {
    my $parms = {
        "object_ref" => $obj_Ecoli,
        "output_genome_name" => "rasted_genome",
        "output_workspace" => $ws
    };
    my $rast_ann;
    throws_ok {
        $rast_ann = $rast_impl->rast_genome_assembly($parms);
    } qr/ERROR calling rast run_pipeline/,
        'Impl rast_genome call returns ERROR due to kmer data absence or other causes.';

    $parms = {
        "object_ref" => $obj_asmb,
        "output_genome_name" => "rasted_assembly",
        "output_workspace" => $ws
    };
    throws_ok {
        $rast_ann = $rast_impl->rast_genome_assembly($parms);
    } qr/ERROR calling rast run_pipeline/,
        'Impl rast_genome call returns ERROR due to kmer data absence or other causes.';
};


## testing generate_genome_report using obj ids from prod ONLY
subtest 'generate_genome_report' => sub {
    my $stats_ok = 'stats generation runs ok.\n';

    my $gff_path = $mgutil->_write_gff_from_genome($obj_asmb_ann);

    my ($gff_contents, $attr_delimiter) = ([], '=');

    my %ret_stats;
    lives_ok {
        ($gff_contents, $attr_delimiter) = $mgutil->_parse_gff($gff_path, $attr_delimiter);
        %ret_stats = $mgutil->_generate_stats_from_gffContents($gff_contents);
        #print "Stats on $obj_asmb_ann: \n".Dumper(\%ret_stats);
    } $stats_ok;
    is(keys %ret_stats, 2, "_generate_stats_from_gffContents on $obj_asmb_ann should return non-empty.\n");
    ok(exists($ret_stats{gene_role_map}), '_generate_stats_from_gffContents stats contains gene_roles.');
    ok(exists($ret_stats{function_roles}), '_generate_stats_from_gffContents stats contains function roles.');

    my %ftr_tab = $mgutil->_get_feature_function_lookup($test_ftrs);
    #print "\nFeature lookup:\n".Dumper(\%ftr_tab);
    my $ret_rpt = $mgutil->_generate_genome_report($obj_asmb, $obj_asmb_ann, [],
                                            $gff_contents, \%ftr_tab);
    print "Report return: \n".Dumper($ret_rpt);
    ok( exists($ret_rpt->{report_ref}), 'Report generation returns report_ref.');
    ok( exists($ret_rpt->{report_name}), 'Report generation returns report_name.');
    ok( exists($ret_rpt->{output_genome_ref}), 'Report generation returns output_gemome_ref.');
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
    lives_ok {
        #$params->{input_genomes} = ["31020/5/1"]; # an appdev object
        $params->{input_genomes} = ["48109/9/1"]; # an prod object
        $params->{genome_text} = '';
        my $ret_ann4 = $rast_impl->annotate_genomes($params);
    } 'Should not throw error due to blank genome_text AND non-empty input_genoms';
};

subtest 'rast_genomes_assemblies' => sub {
    my $error_message = qr/ERROR:Missing required inputs/;
    my $error_mand = qr/Mandatory arguments missing/;

    my $params = {
        "output_GenomeSet_name" => "out_genomeSet"
    };

    throws_ok {
        $params->{output_workspace} = get_ws_name();
        $params->{input_genomes} = [$obj_Ecoli]; # array of prod objects
        $params->{input_assemblies} = [$obj_asmb]; # array of prod objects
        $params->{input_text} = '';
        my $ret_ann6 = $rast_impl->rast_genomes_assemblies($params);
    } qr/ERROR calling rast run_pipeline/,
        'metag_utils rast_genome call returns ERROR due to kmer data absence or other causes.';

    throws_ok {
        $params->{output_workspace} = get_ws_name();
        #$params->{input_genomes} = ["31020/5/1"]; # array of an appdev object
        $params->{input_genomes} = [$obj_Echinacea, $obj_Ecoli]; # array of prod objects
        $params->{input_assemblies} = [];
        $params->{input_text} = '';
        my $ret_ann7 = $rast_impl->rast_genomes_assemblies($params);
    } qr/ERROR calling rast run_pipeline/,
	'metag_utils rast_genome call returns ERROR due to kmer data absence or other causes.';

    throws_ok {
        $params->{output_workspace} = get_ws_name();
        $params->{input_assemblies} = [$obj_asmb_refseq, $obj_asmb]; # array of prod objects
	$params->{input_genomes} = [];
        $params->{input_text} = '';
        my $ret_ann8 = $rast_impl->rast_genomes_assemblies($params);
    } qr/ERROR calling rast run_pipeline/,
        'metag_utils rast_genome call returns ERROR due to kmer data absence or other causes.';

    throws_ok {
        $params->{output_workspace} = get_ws_name();
        $params->{input_genomes} = [$obj_Echinacea, $obj_Ecoli]; # array of prod objects
        $params->{input_assemblies} = [$obj_asmb_refseq, $obj_asmb]; # array of prod objects
        $params->{input_text} = '';
        my $ret_ann9 = $rast_impl->rast_genomes_assemblies($params);
    } qr/ERROR calling rast run_pipeline/,
        'metag_utils rast_genome call returns ERROR due to kmer data absence or other causes.';
};


#----- For checking the stats of a given obj id in prod ONLY-----#
my $stats_ok = 'stats generation runs ok.\n';

subtest '_generate_stats_from_aa & from_gffContents' => sub {
    my ($gff_contents, $attr_delimiter) = ([], '=');
    my $gff_path = '';

    # $obj8
    my %ret_stats = $mgutil->_generate_stats_from_aa($obj8);
    print "Stats from AMA on $obj8:\n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generated from AMA $obj8.");

    $gff_path = $mgutil->_write_gff_from_ama($obj8);
    ($gff_contents, $attr_delimiter) = $mgutil->_parse_gff($gff_path, $attr_delimiter);
    %ret_stats = $mgutil->_generate_stats_from_gffContents($gff_contents);
    #print "Stats from GFF on $obj8: \n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generated from gffContents on $obj8.");

    # $obj9
    %ret_stats = $mgutil->_generate_stats_from_aa($obj9);
    print "Stats from AMA on $obj9:\n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generated from AMA $obj9.");

    $gff_path = $mgutil->_write_gff_from_ama($obj9);
    ($gff_contents, $attr_delimiter) = $mgutil->_parse_gff($gff_path, $attr_delimiter);
    %ret_stats = $mgutil->_generate_stats_from_gffContents($gff_contents);
    #print "Stats from GFF on $obj9: \n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generated from gffContents on $obj9.");

    # $obj10
    %ret_stats = $mgutil->_generate_stats_from_aa($obj10);
    print "Stats from ama on $obj10:\n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generated from AMA $obj10.");

    $gff_path = $mgutil->_write_gff_from_ama($obj10);
    ($gff_contents, $attr_delimiter) = $mgutil->_parse_gff($gff_path, $attr_delimiter);
    %ret_stats = $mgutil->_generate_stats_from_gffContents($gff_contents);
    print "Stats from GFF on $obj10: \n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generated from gffContents on $obj10.");
};

# testing generate_metag_report using obj ids from prod ONLY
subtest 'generate_metag_report' => sub {
    my $stats_ok = 'stats generation runs ok.\n';

    my $gff_path1 = $mgutil->_write_gff_from_ama($obj6);

    my $gff_path2 = $mgutil->_write_gff_from_ama($obj7);

    my ($gff_contents1, $attr_delimiter) = ([], '=');
    my %ret_stats1;
    lives_ok {
        ($gff_contents1, $attr_delimiter) = $mgutil->_parse_gff($gff_path1, $attr_delimiter);
        %ret_stats1 = $mgutil->_generate_stats_from_gffContents($gff_contents1);
        #print "Stats on $obj6: \n".Dumper(\%ret_stats1);
    } $stats_ok;
    is(keys %ret_stats1, 0, "_generate_stats_from_gffContents on $obj6 should return empty.\n");

    my $gff_contents2 = [];
    my %ret_stats2;
    lives_ok {
        ($gff_contents2, $attr_delimiter) = $mgutil->_parse_gff($gff_path2, $attr_delimiter);
        %ret_stats2 = $mgutil->_generate_stats_from_gffContents($gff_contents2);
        #print "Stats on $obj7: \n".Dumper(\%ret_stats2);
    } $stats_ok;
    is(keys %ret_stats2, 2, "_generate_stats_from_gffContents on $obj7 should return non-empty.\n");
    ok(exists($ret_stats2{gene_role_map}), '_generate_stats_from_gffContents stats contains gene_roles.');
    ok(exists($ret_stats2{function_roles}), '_generate_stats_from_gffContents stats contains function roles.');

    my %ftr_tab = $mgutil->_get_feature_function_lookup($test_ftrs);
    #print "\nFeature lookup:\n".Dumper(\%ftr_tab);
    my $ret_rpt = $mgutil->generate_metag_report($obj6, $obj7, $gff_contents1,
                                            $gff_contents2, \%ftr_tab);
    print "Report return: \n".Dumper($ret_rpt);
    ok( exists($ret_rpt->{report_ref}), 'Report generation returns report_ref.');
    ok( exists($ret_rpt->{report_name}), 'Report generation returns report_name.');
    ok( exists($ret_rpt->{output_genome_ref}), 'Report generation returns output_gemome_ref.');
};


# testing rast_metagenome using obj ids from appdev ONLY
subtest 'rast_metagenome_appdev' => sub {
    # an appdev assembly
    my $parms = {
        "object_ref" => "37798/14/1",
        "output_metagenome_name" => "rasted_shortOne_appdev",
        "output_workspace" => $ws
    };
    my $rast_mg_ref = $mgutil->rast_metagenome($parms);
    print "rast_metagenome returns: $rast_mg_ref" if defined($rast_mg_ref);
    ok (($rast_mg_ref !~ m/[^\\w\\|._-]/), 'rast_metagenome returns an INVALID ref');

    # an appdev genome
    my $parms = {
        "object_ref" => "37798/15/1",
        "output_metagenome_name" => "rasted_shortOne_appdev",
        "output_workspace" => $ws
    };
    my $rast_mg_ref = $mgutil->rast_metagenome($parms);
    print "rast_metagenome returns: $rast_mg_ref" if defined($rast_mg_ref);
    ok (($rast_mg_ref !~ m/[^\\w\\|._-]/), 'rast_metagenome returns an INVALID ref');
    $parms = {
        object_ref => $ret_metag->{metagenome_ref},
        output_metagenome_name => 'rasted_metagenome',
        output_workspace => $ws
    };

    throws_ok {
        $rast_mg_ref = $mgutil->rast_metagenome($parms);
    } qr/**rast_metagenome ERROR/,
        'calling rast_metagenome dies file not found';

    throws_ok {
        $rast_mg_ref = $mgutil->rast_metagenome($parms);
        print "rast_metagenome returns: $rast_mg_ref" if defined($rast_mg_ref);
        if ($rast_mg_ref !~ m/[^\\w\\|._-]/) {
            croak "Invalid metagenome object reference:$rast_mg_ref.";
        }
    } qr/Invalid metagenome object reference/,
        'calling rast_metagenome fails to generate a valid metagenome';
};

# testing rast_metagenome using obj ids from prod ONLY
subtest 'rast_metagenome_prod' => sub {
    # a prod assembly
    my $parms = {
        "object_ref" => $obj2_1,
        "output_metagenome_name" => "rasted_ama",
        "output_workspace" => $ws
    };

    my $rast_mg_ref;
    throws_ok {
        $rast_mg_ref = $mgutil->rast_metagenome($parms);
    } qr/ERROR calling rast run_pipeline/,
        'mgutil rast_metagenome call returns ERROR due to kmer data absence or other causes.';

    # a prod metagenome assembly
    $parms = {
        "object_ref" => $obj4,
        "output_metagenome_name" => "rasted_obj4_prod",
        "output_workspace" => $ws
    };
    throws_ok {
        $rast_mg_ref = $mgutil->rast_metagenome($parms);
    } qr/ERROR calling rast run_pipeline/,
        'mgutil rast_metagenome call returns ERROR due to kmer data absence or other causes.';
};

subtest '_prepare_genome_4annotation' => sub {
    # a prod assembly
    my ($gff, $gn, $fa_file, $gff_file);

    throws_ok {
        $fa_file = $mgutil->_get_fasta_from_assembly($obj3);
        $gff_file = $mgutil->_write_gff_from_ama($obj3);
        ($gff, $gn) = $mgutil->_prepare_genome_4annotation($fa_file, $gff_file);
    } qr/ValueError: Object is not an AnnotatedMetagenomeAssembly/,
      "ValueError: Object is not an AnnotatedMetagenomeAssembly because $obj3 did not point to an AMA.";

    # a prod metagenome assembly
    lives_ok {
        $fa_file = $mgutil->_write_fasta_from_ama($obj4);
        $gff_file = $mgutil->_write_gff_from_ama($obj4);
        ($gff, $gn) = $mgutil->_prepare_genome_4annotation($fa_file, $gff_file);
    } 'mgutil->_prepare_genome_4annotation returns normally.';
    ok (@{$gff} > 0, "_prepare_genome_4annotation returns ". scalar @{$gff}." lines of GFF contents.");
    ok (@{$gn->{features}} > 0,
        "_prepare_genome_4annotation returns genome with ". scalar @{$gn->{features}}." features.");
};

# testing generate_metag_report using obj ids from prod ONLY
subtest '_generate_stats_from_aa' => sub {
    my %ret_stats = $mgutil->_generate_stats_from_aa($obj4);
    #print "AMA stats return: \n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generation from AMA $obj4 returns result.");

    %ret_stats = $mgutil->_generate_stats_from_aa($obj5);
    #print "AMA stats return: \n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generation from AMA $obj5 returns result.");

    %ret_stats = $mgutil->_generate_stats_from_aa($obj7);
    #print "AMA stats return: \n".Dumper(\%ret_stats);
    ok(keys %ret_stats, "Statistics generation from AMA $obj7 returns result.");
};

# testing _generate_stats_from_gffContents using obj ids from prod ONLY
subtest '_generate_stats_from_gffContents' => sub {
    my $gff_path = $mgutil->_write_gff_from_ama($obj7);

    my ($gff_contents, $attr_delimiter) = ([], '=');
    ($gff_contents, $attr_delimiter) = $mgutil->_parse_gff($gff_path, $attr_delimiter);

    my %ret_stats = $mgutil->_generate_stats_from_gffContents($gff_contents);
    ok(keys %ret_stats, 'Statistics generation from gff_contents returns result.');
};

# test by using appdev obj id
subtest 'annotate_metagenome' => sub {
    my $parms = {
        "object_ref" => "37798/7/1",
        "output_metagenome_name" => "rasted_shortOne_appdev",
        "output_workspace" => $ws
    };
    my $rast_ann = $rast_impl->annotate_metagenome($parms);
    print Dumper($rast_ann);

};

# metagenome saved successfully by using appdev obj ids
subtest '_save_metagenome' => sub {
    my $mymetag = {};
    lives_ok {
        $mymetag = $mgutil->_save_metagenome(
                       $ws, $out_name, $obj1, $gff1);
    } '__save_metagenome run without errors on short_one.\n';
    ok (exists $mymetag->{metagenome_ref},
        "metagenome saved with metagenome_ref='$mymetag->{metagenome_ref}'");
    ok (exists $mymetag->{metagenome_info}, 'metagenome saved with metagenome_info');
    is ($mymetag->{metagenome_info}[1], $out_name, 'saved metagenome name is correct');
    is ($mymetag->{metagenome_info}[7], $ws, 'saved metagenome to the correct workspace');
    
    lives_ok {
        $mymetag = $mgutil->_save_metagenome(
                       $ws, $out_name, $obj2, $gff2);
    } '_save_metagenome runs without errors on 59111.assembled.\n';
    ok (exists $mymetag->{metagenome_ref},
        "metagenome saved with metagenome_ref='$mymetag->{metagenome_ref}'");
    ok (exists $mymetag->{metagenome_info}, 'metagenome saved with metagenome_info');
    is ($mymetag->{metagenome_info}[1], $out_name, 'saved metagenome name is correct');
    is ($mymetag->{metagenome_info}[7], $ws, 'saved metagenome to the correct workspace');
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
