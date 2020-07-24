use Test::Most;
use RASTTestUtils;
use Data::Dumper::Concise;
use feature qw( say );
use File::Copy;

print "PURPOSE:\n";
print "    1.  Test annotation of one small assembly. \n";
print
    "    2.  Test that the saved genome isn't using defaults. Must be Archaea and genetic code 4\n";
print "    3.  Test passing a taxon ID vs. a scientific name.\n";
print "    4.  Test 2 failure modes for contacting the RE.\n";
print "\n";

my $ws_client = RASTTestUtils::get_ws_client();
my $ws_name   = RASTTestUtils::get_ws_name();

my $assembly_obj_name = "Acidilobus_sp._CIS.fna";
my $assembly_ref      = RASTTestUtils::prepare_assembly( $assembly_obj_name );
my $genome_obj_name   = 'Acidilobus_sp_CIS';

my $params = {
    "input_contigset"            => $assembly_obj_name,
    "scientific_name"            => 'Acidilobus sp 7',
    "domain"                     => 'A',
    "genetic_code"               => '4',
    "call_features_CDS_prodigal" => '1',
};

$params = RASTTestUtils::set_params( $genome_obj_name, $params );

# Test processing a single assembly and saving as a genome.

subtest 'Running RAST annotation' => sub {

    lives_ok {
        RASTTestUtils::make_impl_call( "RAST_SDK.annotate_genome", $params );
    } 'annotate_genome runs successfully';

    my $genome_ref = $ws_name . "/" . $genome_obj_name;
    my $genome_obj = $ws_client->get_objects( [ { ref => $genome_ref } ] )->[ 0 ]{ data };

    say "\n\nOUTPUT OBJECT DOMAIN = $genome_obj->{domain}";
    say "OUTPUT OBJECT G_CODE = $genome_obj->{genetic_code}\n";

    ok( defined( $genome_obj->{ features } ),       "Features array is present" );
    ok( scalar @{ $genome_obj->{ features } } gt 0, "Number of features" );
    ok( defined( $genome_obj->{ cdss } ),           "CDSs array is present" );
    ok( scalar @{ $genome_obj->{ cdss } } gt 0,     "Number of CDSs" );
    ok( defined( $genome_obj->{ mrnas } ),          "mRNAs array is present" );
    ok( scalar @{ $genome_obj->{ mrnas } } gt 0,    "Number of mRNAs" );
    ok( $genome_obj->{ scientific_name } eq "Acidilobus sp 7", "Sci name is correct" );
    ok( !defined( $genome_obj->{ taxon_assignments } ),
        "Taxon assignments is undefined" );
};

# Test processing a single assembly and saving as a genome while specifing a NCBI taxon ID.
subtest 'Running RAST annotation with Taxon ID' => sub {

    my $params_copy = { %$params };

    # this tax ID's species name changed in the 2018-12 NCBI dump and again in the 2019-02 dump
    # so it is a good test case for making sure the timestamp is passed to the RE correctly.
    # It depends on the RE containing 2018 NCBI data, which it currently does
    $params_copy->{ ncbi_taxon_id }                = 2448083;
    $params_copy->{ relation_engine_timestamp_ms } = 1545000000000;    # epoch ms

    lives_ok {
        RASTTestUtils::make_impl_call( "RAST_SDK.annotate_genome", $params_copy );
    } 'annotate_genomes runs successfully';
    my $genome_ref = $ws_name . "/" . $genome_obj_name;
    my $genome_obj = $ws_client->get_objects( [ { ref => $genome_ref } ] )->[ 0 ]{ data };

    say "\n\nOUTPUT OBJECT DOMAIN = $genome_obj->{domain}";
    say "OUTPUT OBJECT G_CODE = $genome_obj->{genetic_code}\n";

    ok( defined( $genome_obj->{ features } ),       "Features array is present" );
    ok( scalar @{ $genome_obj->{ features } } gt 0, "Number of features" );
    ok( defined( $genome_obj->{ cdss } ),           "CDSs array is present" );
    ok( scalar @{ $genome_obj->{ cdss } } gt 0,     "Number of CDSs" );
    ok( defined( $genome_obj->{ mrnas } ),          "mRNAs array is present" );
    ok( scalar @{ $genome_obj->{ mrnas } } gt 0,    "Number of mRNAs" );
    is $genome_obj->{ scientific_name }, "Nectria sp. MJH 2018c",
        "Sci name is correct";
    cmp_deeply $genome_obj->{ taxon_assignments },
        { 'ncbi' => '2448083' },
        "Taxon assignments are correct";

};

subtest 'RAST annotation into GenomeSet' => sub {

    # Test processing a single assembly and saving as a genome set.

    my $genome_set_name = "New_GenomeSet";
    my $new_params      = {
        "input_genomes"               => [ $assembly_ref ],
        "call_features_tRNA_trnascan" => '1',
        "output_genome"               => $genome_set_name
    };

    my ( $genome_set_obj, $more_params )
        = RASTTestUtils::submit_set_annotation( $genome_set_name,
        $new_params->{ input_genomes }, $new_params );
    my $data = $ws_client->get_objects( [ { ref => $genome_set_obj } ] )->[ 0 ]{ refs };
    my $number_genomes = scalar @{ $data };
    ok( $number_genomes == 1,
        "Input: One Assembly. Output: $number_genomes in output GenomeSet" );

    my $genome_obj
        = $ws_client->get_objects( [ { ref => $data->[ 0 ] } ] )->[ 0 ]{ data };
    ok( $genome_obj->{ scientific_name } eq "unknown taxon", "Sci name is correct" );
    ok( !defined( $genome_obj->{ taxon_assignments } ),
        "Taxon assignments is undefined" );

    # I have no idea what this test is supposed to prove - 19/10/30
    my $report     = "/kb/module/work/tmp/annotation_report.$genome_set_name";
    my $directory  = "/kb/module/test/report_output/";
    my $local_path = $directory . "annotation_report.$genome_set_name";

    unless ( mkdir $directory ) { die "Unable to create directory " . $directory; }

    copy $report, $local_path or die "copy failed: $!";

    ok -e $local_path, 'File found';

};

# Test processing a single assembly with a taxon ID and saving as a genome set.
subtest 'Running RAST annotation into GenomeSet with Taxon ID' => sub {
    my $genome_set_name = "New_GenomeSet2";
    my $new_params      = {
        "input_genomes"               => [ $assembly_ref ],
        "call_features_tRNA_trnascan" => '1',
        "output_genome"               => $genome_set_name,

        # this tax ID's species name changed in the 2018-12 NCBI dump and again in
        # the 2019-02 dump so it is a good test case for making sure the timestamp
        # is passed to the RE correctly. It depends on the RE containing 2018 NCBI
        # data, which it currently does
        "ncbi_taxon_id"                => 2448083,
        "relation_engine_timestamp_ms" => 1545000000000    # epoch ms
    };

    my ( $genome_set_obj, $more_params )
        = RASTTestUtils::submit_set_annotation( $genome_set_name,
        $new_params->{ input_genomes }, $new_params );
    my $data = $ws_client->get_objects( [ { ref => $genome_set_obj } ] )->[ 0 ]{ refs };
    my $number_genomes = scalar @{ $data };
    ok( $number_genomes == 1,
        "Input: One Assembly. Output: $number_genomes in output GenomeSet" );

    my $genome_obj
        = $ws_client->get_objects( [ { ref => $data->[ 0 ] } ] )->[ 0 ]{ data };
    is $genome_obj->{ scientific_name },
        "Nectria sp. MJH 2018c",
        "Sci name is correct";
    cmp_deeply(
        $genome_obj->{ taxon_assignments },
        { 'ncbi' => '2448083' },
        "Taxon assignments is correct"
    );

    # I have no idea what this test is supposed to prove - 19/10/30
    my $report     = "/kb/module/work/tmp/annotation_report.$genome_set_name";
    my $directory  = "/kb/module/test/report_output2/";
    my $local_path = $directory . "annotation_report.$genome_set_name";

    unless ( mkdir $directory ) { die "Unable to create directory " . $directory; }

    copy $report, $local_path or die "copy failed: $!";

    ok -e $local_path, 'File found';
};

# Test failing to contact the RE with bad input.
subtest 'Running RAST annotation fail with bad RE input' => sub {
    my $params_copy = { %$params };
    $params_copy->{ ncbi_taxon_id }                = 32;
    $params_copy->{ relation_engine_timestamp_ms } = 'Sept 19 2020';    # oops
      # the error here is a JSON string appended with a file and line number, and so can't be
      # decoded. Arrg.
    throws_ok {
        RASTTestUtils::make_impl_call( "RAST_SDK.annotate_genome", $params_copy );
    } qr/Error contacting Relation Engine: 400 BAD REQUEST/,
        'annotate_genome fails with appropriate message';
};

# Test sending a non-existant taxon to the RE.
subtest 'Running RAST annotation fail with bad tax ID' => sub {
    my $params_copy = { %$params };
    $params_copy->{ ncbi_taxon_id } = 1000000000000;    # pretty sure there aren't 1T taxa
    $params_copy->{ relation_engine_timestamp_ms } = 1572648527000;
    throws_ok {
        RASTTestUtils::make_impl_call( "RAST_SDK.annotate_genome", $params_copy );
    } qr/No result from Relation Engine for NCBI taxonomy ID 1000000000000/,
        'annotate_genome fails with appropriate message';
};

RASTTestUtils::clean_up();

done_testing;

