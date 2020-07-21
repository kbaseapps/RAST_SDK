use Test::Most;
use RASTTestUtils;
use Data::Dumper::Concise;
use feature qw( say );

my $DEBUG = 'N';

unless ( $DEBUG eq 'N' ) {

    # skip all the tests
    ok 'skipping tests in assemblySet.t';
    done_testing;
    exit 0;
}

say join "\n", (
    "PURPOSE:",
    "  1.  Test annotate Multiple Assemblies.",
    "      Test two assemblies, an assemblySet, an asseblySet plus singleton, and redundant assemblies",
    "  2.  Minimum test for tRNA",
    "  3.  In debug mode, using a genome reference in CI/prod.",
    "      Otherwise, load an assembly in data dir. This takes more time but every user has access.",
    "      For this reason, the tests aren't for a specific number of changes or names.",
);

my $ws_client     = RASTTestUtils::get_ws_client();
my $ws_name       = RASTTestUtils::get_ws_name();
my $su            = RASTTestUtils::get_setutils_client();

my $assembly_obj_name1 = "Dactylopius coccus";
my $assembly_ref1      = "40046/5/1";
$assembly_ref1 = "40619/11/1";
my $assembly_obj_name2 = "Drosophila melanogaster";
my $assembly_ref2      = "40046/6/1";
$assembly_ref2 = "40619/39/1";
my $assembly_obj_name3 = "Nomada ferruginata";
my $assembly_ref3      = "40046/7/1";
my $genome_set_name    = "New_GenomeSet";

if ( $DEBUG ne 'Y' ) {
    $assembly_obj_name1 = "bogus.fna";
    $assembly_ref1      = RASTTestUtils::prepare_assembly( $assembly_obj_name1 );

    $assembly_obj_name2 = "bogus2.fna";
    $assembly_ref2      = RASTTestUtils::prepare_assembly( $assembly_obj_name2 );

    $assembly_obj_name3 = "bogus2.fna";
    $assembly_ref3      = RASTTestUtils::prepare_assembly( $assembly_obj_name2 );
}

my $params = {
    "input_genomes"               => '',
    "call_features_tRNA_trnascan" => '1',
    "output_genome"               => $genome_set_name
};

#
#   Set up the needed AssemblySet
#
my $assembly_set_name = 'new_assembly_set';
my $assembly_set      = $su->KButil_Build_AssemblySet( {
    workspace_name => $ws_name,
    input_refs     => [ $assembly_ref1, $assembly_ref2 ],
    output_name    => $assembly_set_name,
    desc           => 'Test AssemblySet Description'
} );

#
#   TEST TWO ASSEMBLIES -
#
say "ASSEMBLYREF1 = $assembly_ref1 and ASSEMBLYREF2 = $assembly_ref2";

subtest 'Two assemblies' => sub {

    my $genome_refs = [ $assembly_ref1, $assembly_ref2 ];
    $params->{ input_genomes } = $genome_refs;
    my ( $genome_set_obj ) = RASTTestUtils::submit_set_annotation(
        $genome_set_name, $genome_refs, $params
    );

    my $data = $ws_client->get_objects( [ { ref => $genome_set_obj } ] )->[ 0 ]{ refs };
    ok @$data == 2, "Correct number of genomes in output set"
        or diag explain $data;
};

subtest 'AssemblySet with two Assemblies' => sub {

    my $gs     = $ws_name . "/" . $assembly_set_name;
    my $info   = $ws_client->get_objects( [ { ref => $gs } ] )->[ 0 ]{ info };
    my $newref = $info->[ 6 ] . "/" . $info->[ 0 ] . "/" . $info->[ 4 ];
    $params->{ input_genomes } = [ $newref ];

    my ( $genome_set_obj ) = RASTTestUtils::submit_set_annotation( $genome_set_name,
        $params->{ input_genomes }, $params );
    my $data = $ws_client->get_objects( [ { ref => $genome_set_obj } ] )->[ 0 ]{ refs };
    ok @$data == 2, "Correct number of genomes in output set"
        or diag explain $data;
};

subtest "AssemblySet plus single assembly" => sub {

    my $gs     = $ws_name . "/" . $assembly_set_name;
    my $info   = $ws_client->get_objects( [ { ref => $gs } ] )->[ 0 ]{ info };
    my $newref = $info->[ 6 ] . "/" . $info->[ 0 ] . "/" . $info->[ 4 ];
    $params->{ input_genomes } = [ $newref, $assembly_ref3 ];

    my ( $genome_set_obj ) = RASTTestUtils::submit_set_annotation(
        $genome_set_name, $params->{ input_genomes }, $params
    );
    my $data = $ws_client->get_objects( [ { ref => $genome_set_obj } ] )->[ 0 ]{ refs };
    ok @$data == 3,"Correct number of genomes in output set"
        or diag explain $data;

};

subtest 'Redundant assemblies' => sub {

    $params->{ input_genomes } = [ $assembly_ref1, $assembly_ref1 ];
    my ( $genome_set_obj )     = RASTTestUtils::submit_set_annotation(
        $genome_set_name, $params->{ input_genomes }, $params
    );
    my $data = $ws_client->get_objects( [ { ref => $genome_set_obj } ] )->[ 0 ]{ refs };
    ok @$data == 1, "Correct number of genomes in output set"
        or diag explain $data;
};

RASTTestUtils::clean_up();

done_testing;


