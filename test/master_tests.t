use Test::Most;
use Test::Compile;
use RASTTestUtils;
use Data::Dumper::Concise;
use feature qw( say );

#Use this control to turn various tests on and off
my $master_test_list = {
	"Make sure all code compiles" => 1,
	"Two assemblies" => 1,
	"AssemblySet with two Assemblies" => 1,
	"AssemblySet plus single assembly" => 1,
	"Redundant assemblies" => 1,
	"Running RAST annotation prodigal" => 1,
	'Running RAST annotation glimmer3' => 1,
	"Running RAST annotation both" => 1,
	"Running RAST annotation" => 1,
	"Running RAST annotation with Taxon ID" => 1,
	"RAST annotation into GenomeSet" => 1,
	"Running RAST annotation into GenomeSet with Taxon ID" => 1,
	"Running RAST annotation fail with bad RE input" => 1,
	"Running RAST annotation fail with bad tax ID" => 1
};

# Global parameters used in tests
my $ws_client = RASTTestUtils::get_ws_client();
my $su = RASTTestUtils::get_setutils_client();
my $tempws = 0;
my $test_data_workspace = "68255";
my $assembly_object_name = "GCA_000350285.1_OR1_genomic.fna";
my $assembly_object_ref = $test_data_workspace."/".$assembly_object_name;
my $genome_obj_name   = 'SR1_bacterium_MGEHA_gc25';
my $output_ws = "chenry:narrative_1663301769839";#Either ask chenry for write permission or reset this to RASTTestUtils::get_ws_id()
my $tester = Test::Compile->new();
RASTTestUtils::set_ws($output_ws);

#Tests start here:
say qq(PURPOSE:
    1. Make sure all code compiles
);
if ($master_test_list->{"Make sure all code compiles"} == 1) {
	# check all .pm and .pl files compile
	$tester->all_files_ok( '/kb/module/' );
	
	# also check .t and .psgi files
	my @all_ext_files = RASTTestUtils::files_with_ext( $tester, qr/\.(t|psgi)$/, '/kb/module/' );
	
	for ( @all_ext_files ) {
	    ok $tester->pl_file_compiles( $_ ), $_ . ' compiles';
	}
}

say qq(PURPOSE:
	1.  Test annotate Multiple Assemblies. Test two assemblies, an assemblySet, an asseblySet plus singleton, and redundant assemblies
	2.  Minimum test for tRNA
);
my $params = {
    "call_features_tRNA_trnascan" => '1',
    "output_genome"               => "New_GenomeSet"
};
#   TEST TWO ASSEMBLIES -
if ($master_test_list->{"Two assemblies"} == 1) {
	subtest 'Two assemblies' => sub {
	    $params->{ input_genomes } = [ $test_data_workspace."/bogus.fna", $test_data_workspace."/bogus2.fna" ];
	    my ( $genome_set_obj ) = RASTTestUtils::submit_set_annotation($params);
	    my $data = $ws_client->get_objects( [ { ref => $genome_set_obj } ] )->[ 0 ]{ refs };
	    ok @$data == 2, "Correct number of genomes in output set"
	        or diag explain $data;
	};
}

if ($master_test_list->{"AssemblySet with two Assemblies"} == 1) {
	subtest 'AssemblySet with two Assemblies' => sub {
	    $params->{ input_genomes } = [ $test_data_workspace . "/new_assembly_set" ];
	    my ( $genome_set_obj ) = RASTTestUtils::submit_set_annotation($params);
	    my $data = $ws_client->get_objects( [ { ref => $genome_set_obj } ] )->[ 0 ]{ refs };
	    ok @$data == 2, "Correct number of genomes in output set"
	        or diag explain $data;
	};
}

if ($master_test_list->{"AssemblySet plus single assembly"} == 1) {
	subtest "AssemblySet plus single assembly" => sub {
	    $params->{ input_genomes } = [ $test_data_workspace . "/new_assembly_set", $assembly_object_ref ];
	    my ( $genome_set_obj ) = RASTTestUtils::submit_set_annotation($params);
	    my $data = $ws_client->get_objects( [ { ref => $genome_set_obj } ] )->[ 0 ]{ refs };
	    ok @$data == 3,"Correct number of genomes in output set"
	        or diag explain $data;
	};
}

if ($master_test_list->{"Redundant assemblies"} == 1) {
	subtest 'Redundant assemblies' => sub {
	    $params->{ input_genomes } = [ $test_data_workspace."/bogus.fna", $test_data_workspace."/bogus.fna" ];
	    my ( $genome_set_obj )     = RASTTestUtils::submit_set_annotation($params);
	    my $data = $ws_client->get_objects( [ { ref => $genome_set_obj } ] )->[ 0 ]{ refs };
	    ok @$data == 1, "Correct number of genomes in output set"
	        or diag explain $data;
	};
}

say qq(PURPOSE:
    1.  Test annotation of genetic code 25 for Mircea Podar.
    2.  Test the Prodigal is using the -meta when the domain is 'Unknown'
);
if ($master_test_list->{"Running RAST annotation prodigal"} == 1) {
	subtest 'Running RAST annotation prodigal' => sub {
	    my $params = {
	        "input_contigset"            => $assembly_object_ref,
	        "scientific_name"            => 'candidate division SR1 bacterium MGEHA',
	        "domain"                     => 'U',
	        "workspace"                  => $output_ws,
	        "output_genome"              => $genome_obj_name.'_prodigal',
	        "genetic_code"               => '25',
	        "call_features_CDS_prodigal" => '1'
	    };
	
	    $params = RASTTestUtils::set_params( $genome_obj_name, $params );
	
	    RASTTestUtils::make_impl_call( "RAST_SDK.annotate_genome", $params );
	    my $genome_ref = $output_ws . "/" . $params->{output_genome};
	    my $genome_obj = $ws_client->get_objects( [ { ref => $genome_ref } ] )->[ 0 ]{ data };
	
	    say "OUTPUT OBJECT DOMAIN = $genome_obj->{domain}\n"
	    . "OUTPUT OBJECT G_CODE = $genome_obj->{genetic_code}\n"
	    . "OUTPUT SCIENTIFIC NAME = $genome_obj->{scientific_name}\n";
	
	    ok $genome_obj->{ features } && @{ $genome_obj->{ features } },
	        "Features array is present and contains more than zero items"
	        or diag explain $genome_obj;
	};
}

if ($master_test_list->{"Running RAST annotation prodigal"} == 1) {
	subtest 'Running RAST annotation glimmer3' => sub {
	    my $params = {
	        "input_contigset"            => $assembly_object_ref,
	        "scientific_name"            => 'candidate division SR1 bacterium MGEHA',
	        "domain"                     => 'U',
	        "workspace"                  => $output_ws,
	        "output_genome"              => $genome_obj_name.'_glimmer3',
	        "genetic_code"               => '25',
	        "call_features_CDS_glimmer3" => '1'
	    };
	
	    $params = RASTTestUtils::set_params( $genome_obj_name, $params );
	
	    ## Before fixing the GC25 scenario the following block should throw and catch an error.
	    #throws_ok {
	    #    RASTTestUtils::make_impl_call( "RAST_SDK.annotate_genome", $params );
	    #} qr/Error invoking method call_features_CDS_glimmer3/,
	    #    'RAST_SDK.annotate_genome dies with "Could not extract training ORFs from contig"';
	    #
	
	    ## Now it should do, instead of glimmer3, prodigal calling
	    RASTTestUtils::make_impl_call( "RAST_SDK.annotate_genome", $params );
	    my $genome_ref = $output_ws . "/" . $params->{output_genome};
	    my $genome_obj = $ws_client->get_objects( [ { ref => $genome_ref } ] )->[ 0 ]{ data };
	
	    ok $genome_obj->{ features } && @{ $genome_obj->{ features } },
	        "Glimmer3 gene call was replaced by prodigal and features array is present and contains more than zero items"
	        or diag explain $genome_obj;
	};
}

if ($master_test_list->{"Running RAST annotation both"} == 1) {
	subtest 'Running RAST annotation both' => sub {
	    my $params = {
	        "input_contigset"            => $assembly_object_ref,
	        "scientific_name"            => 'candidate division SR1 bacterium MGEHA',
	        "domain"                     => 'U',
	        "workspace"                  => $output_ws,
	        "output_genome"              => $genome_obj_name.'_both',
	        "genetic_code"               => '25',
	        "call_features_CDS_prodigal" => '1',
	        "call_features_CDS_glimmer3" => '1'
	    };
	
	    $params = RASTTestUtils::set_params( $genome_obj_name, $params );
	
	    RASTTestUtils::make_impl_call( "RAST_SDK.annotate_genome", $params );
	    my $genome_ref = $output_ws . "/" . $params->{output_genome};
	    my $genome_obj = $ws_client->get_objects( [ { ref => $genome_ref } ] )->[ 0 ]{ data };
	
	    ok $genome_obj->{ features } && @{ $genome_obj->{ features } },
	        "Glimmer3 gene call was skipped and features array is present and contains more than zero items"
	        or diag explain $genome_obj;
	
	};
}

say qq(PURPOSE:
	1.  Test annotation of one small assembly.
	2.  Test that the saved genome isn't using defaults. Must be Archaea and genetic code 4
	3.  Test passing a taxon ID vs. a scientific name.
	4.  Test 2 failure modes for contacting the RE.
);
$params = {
    "input_contigset"            => $test_data_workspace."/Acidilobus_sp._CIS.fna",
    "scientific_name"            => 'Acidilobus sp 7',
    "domain"                     => 'A',
    "genetic_code"               => '4',
    "call_features_CDS_prodigal" => '1',
};
$params = RASTTestUtils::set_params( 'Acidilobus_sp_CIS', $params );
if ($master_test_list->{"Running RAST annotation"} == 1) {
	subtest 'Running RAST annotation' => sub {
	
	    lives_ok {
	        RASTTestUtils::make_impl_call( "RAST_SDK.annotate_genome", $params );
	    } 'annotate_genome runs successfully';
	
	    my $genome_ref = $output_ws . "/Acidilobus_sp_CIS";
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
}

# Test processing a single assembly and saving as a genome while specifing a NCBI taxon ID.
if ($master_test_list->{"Running RAST annotation with Taxon ID"} == 1) {
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
	    my $genome_ref = $output_ws . "/Acidilobus_sp_CIS";
	    my $genome_obj = $ws_client->get_objects( [ { ref => $genome_ref } ] )->[ 0 ]{ data };
	
	    say "\n\nOUTPUT OBJECT DOMAIN = $genome_obj->{domain}";
	    say "OUTPUT OBJECT G_CODE = $genome_obj->{genetic_code}\n";
	
	    ok( defined( $genome_obj->{ features } ),       "Features array is present" );
	    ok( scalar @{ $genome_obj->{ features } } gt 0, "Number of features" );
	    ok( defined( $genome_obj->{ cdss } ),           "CDSs array is present" );
	    ok( scalar @{ $genome_obj->{ cdss } } gt 0,     "Number of CDSs" );
	    ok( defined( $genome_obj->{ mrnas } ),          "mRNAs array is present" );
	    ok( scalar @{ $genome_obj->{ mrnas } } gt 0,    "Number of mRNAs" );
	    is $genome_obj->{ scientific_name }, "Metarhizium sp. MJH 2018c",
	        "Sci name is correct";
	    cmp_deeply $genome_obj->{ taxon_assignments },
	        { 'ncbi' => '2448083' },
	        "Taxon assignments are correct";
	
	};
}

if ($master_test_list->{"RAST annotation into GenomeSet"} == 1) {
	subtest 'RAST annotation into GenomeSet' => sub {
	
	    # Test processing a single assembly and saving as a genome set.

	    my $genome_set_name = "New_GenomeSet";
	    my $new_params      = {
	        "input_genomes"               => [ $test_data_workspace."/bogus.fna" ],
	        "call_features_tRNA_trnascan" => '1',
	        "output_genome"               => $genome_set_name
	    };
	
	    my ( $genome_set_obj, $more_params )
	        = RASTTestUtils::submit_set_annotation($new_params );
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
}

# Test processing a single assembly with a taxon ID and saving as a genome set.
if ($master_test_list->{"Running RAST annotation into GenomeSet with Taxon ID"} == 1) {
	subtest 'Running RAST annotation into GenomeSet with Taxon ID' => sub {
	    my $genome_set_name = "New_GenomeSet2";
	    my $new_params      = {
	        "input_genomes"               => [ $test_data_workspace."/bogus.fna" ],
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
	        = RASTTestUtils::submit_set_annotation($new_params );
	    my $data = $ws_client->get_objects( [ { ref => $genome_set_obj } ] )->[ 0 ]{ refs };
	    my $number_genomes = scalar @{ $data };
	    ok( $number_genomes == 1,
	        "Input: One Assembly. Output: $number_genomes in output GenomeSet" );
	
	    my $genome_obj
	        = $ws_client->get_objects( [ { ref => $data->[ 0 ] } ] )->[ 0 ]{ data };
	    is $genome_obj->{ scientific_name },
	        "Metarhizium sp. MJH 2018c",
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
}

# Test failing to contact the RE with bad input.
if ($master_test_list->{"Running RAST annotation fail with bad RE input"} == 1) {
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
}

# Test sending a non-existant taxon to the RE.
if ($master_test_list->{"Running RAST annotation fail with bad tax ID"} == 1) {
	subtest 'Running RAST annotation fail with bad tax ID' => sub {
	    my $params_copy = { %$params };
	    $params_copy->{ ncbi_taxon_id } = 1000000000000;    # pretty sure there aren't 1T taxa
	    $params_copy->{ relation_engine_timestamp_ms } = 1572648527000;
	    throws_ok {
	        RASTTestUtils::make_impl_call( "RAST_SDK.annotate_genome", $params_copy );
	    } qr/No result from Relation Engine for NCBI taxonomy ID 1000000000000/,
	        'annotate_genome fails with appropriate message';
	};
}

if ($tempws == 1) {
	RASTTestUtils::clean_up();
}

done_testing;
