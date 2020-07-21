use Test::Most;
use RASTTestUtils;
use Data::Dumper::Concise;

my $DEBUG = 'N';

unless ( $DEBUG eq 'N' ) {
    ok 'skipping tests in genomesets_and_multiples.t';
    done_testing;
    exit 0;
}

print "Running $0\n";
print "PURPOSE:\n";
print "    1.  Test of Multiple Genomes and GenomeSets. \n";
print "        a. Two individual genomes and annotate.\n";
print "        b. Take the two individual genomes, make a set and annotate.\n";
print "        c. Take the genomeSet, add an individual genomes and annotate.\n";
print "    2.  The optional parameters are minimal because the sets are the test.\n";
print "    3.  In debug mode, using a genome reference in CI/prod.\n";
print
    "         Otherwise, load a genome in data dir. This takes more time but every user has access.\n";
print "         For this reason, the test isn't for a specific number of changes\n\n";

my $ws_client = RASTTestUtils::get_ws_client();
my $ws_name   = RASTTestUtils::get_ws_name();
my $su        = RASTTestUtils::get_setutils_client();

my $genome_obj_name1 = "Dactylopius_coccus";
my $genome_ref1      = "40046/2/1";
my $genome_obj_name2 = "Drosophila_melanogaster";
my $genome_ref2      = "40046/3/1";
my $genome_obj_name3 = "Nomada_ferruginata";
my $genome_ref3      = "40046/4/1";

# PROD
#my $genome_ref1 = "19217/183555/4";
#my $genome_ref2 = "19217/191347/4";
#my $genome_ref3 = "19217/165153/4";

my $tmp_obj;

$genome_obj_name1 = "Carsonella";
my $genome_gbff_name1 = "Carsonella.gbk";
( $tmp_obj, $genome_ref1 ) = RASTTestUtils::prepare_gbff( $genome_gbff_name1, $genome_obj_name1 );

$genome_obj_name2 = "Methanosarcina";
my $genome_gbff_name2 = "Methanosarcina_acetivorans_C2A.gbff";
( $tmp_obj, $genome_ref2 ) = RASTTestUtils::prepare_gbff( $genome_gbff_name2, $genome_obj_name2 );

my $genome_gbff_name3 = "GCF_000287295.1_ASM28729v1_genomic.gbff";
( $tmp_obj, $genome_ref3 ) = RASTTestUtils::prepare_gbff( $genome_gbff_name3, $genome_obj_name3 );

my $params = {
    "input_genomes"               => [ $genome_ref1, $genome_ref2 ],
    "call_features_tRNA_trnascan" => '1',
};

my ( $orig_genome1, $orig_funcs1 ) = RASTTestUtils::get_genome( $genome_ref1 );
my ( $orig_genome2, $orig_funcs2 ) = RASTTestUtils::get_genome( $genome_ref2 );
my ( $orig_genome3, $orig_funcs3 ) = RASTTestUtils::get_genome( $genome_ref3 );

subtest 'two individual genomes' => sub {
    $params->{ input_genomes } = [ $genome_ref1, $genome_ref2 ];
    my ( $genome_set_obj ) = RASTTestUtils::submit_set_annotation( 'multi_genomes',
        $params->{ input_genomes }, $params );

    my $data  = $ws_client->get_objects( [ { ref => $genome_set_obj } ] )->[ 0 ]{ refs };
    ok @$data == 2, "Correct number of genomes in output GenomeSet";
};

subtest "Build and Test a genome set" => sub {
    my $genome_set_name = 'new_genome_set';
    my $genome_set      = $su->KButil_Build_GenomeSet( {
        workspace_name => $ws_name,
        input_refs     => [ $genome_ref1, $genome_ref2 ],
        output_name    => $genome_set_name,
        desc           => 'GenomeSet Description'
    } );

    my $gs     = $ws_name . "/" . $genome_set_name;
    my $info   = $ws_client->get_objects( [ { ref => $gs } ] )->[ 0 ]{ info };
    my $newref = $info->[ 6 ] . "/" . $info->[ 0 ] . "/" . $info->[ 4 ];
    $params->{ input_genomes } = [ $newref ];

    my ( $genome_set_obj ) = RASTTestUtils::submit_set_annotation( $genome_set_name,
        $params->{ input_genomes }, $params );

    my $data = $ws_client->get_objects( [ { ref => $genome_set_obj } ] )->[ 0 ]{ refs };
    ok @$data == 2, "Correct number of genomes in output GenomeSet";
};

subtest 'build and test a genome set plus one' => sub {
    my $genome_set_name = 'genome_plus_set';
    my $genome_set      = $su->KButil_Build_GenomeSet( {
        workspace_name => $ws_name,
        input_refs     => [ $genome_ref1, $genome_ref2 ],
        output_name    => $genome_set_name,
        desc           => 'GenomeSet Description'
    } );

    my $gs     = $ws_name . "/" . $genome_set_name;
    my $info   = $ws_client->get_objects( [ { ref => $gs } ] )->[ 0 ]{ info };
    my $newref = $info->[ 6 ] . "/" . $info->[ 0 ] . "/" . $info->[ 4 ];
    $params->{ input_genomes } = [ $newref, $genome_ref3 ];

    my ( $genome_set_obj ) = RASTTestUtils::submit_set_annotation( $genome_set_name,
        $params->{ input_genomes }, $params );

    my $data = $ws_client->get_objects( [ { ref => $genome_set_obj } ] )->[ 0 ]{ refs };
    ok @$data == 2, "Correct number of genomes in output GenomeSet";
};

RASTTestUtils::clean_up();

done_testing;
