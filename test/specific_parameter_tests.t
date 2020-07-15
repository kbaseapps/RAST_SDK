use Test::Most;
use RASTTestUtils;
use Data::Dumper::Concise;
use feature qw( say );

say "    In debug mode, use a genome reference in CI/prod.";
say "    Otherwise, load a genome in data dir. This takes more time but every user has access.";
say "    For this reason, the tests aren't for a specific number of changes";

my $DEBUG = 'N';

my %list = (
    'carson' => 1,
    'clos'   => 1,
    'meth'   => 1,
    'shew'   => 1,
    'strepp' => 1,
    'streps' => 1
);

subtest 'Carsonella tests' => sub {
    print "PURPOSE:\n";
    print "    1.  Test annotation of one small genome. (Carsonella)\n";
    print "    2.  Test for rRNA, tRNA, and protein_similarity\n";

    my $genome_ref       = "15792/210698/1";
    my $genome_gbff_name = "Carsonella.gbk";
    my $genome_obj_name  = "Carsonella";

    my $params = {
        "input_genome"                 => $genome_ref,
        "call_features_rRNA_SEED"      => '1',
        "call_features_tRNA_trnascan"  => '1',
        "annotate_proteins_similarity" => '1',
    };

    if ( $list{ carson } ) {
        my %counts = annotate( $genome_obj_name, $genome_gbff_name, $genome_ref, $params );
        ok( $counts{ num_diff } > 0,
            "Number of new annotation (" . $counts{ num_diff } . ")" );
        ok( $counts{ num_new } > 0, "Number of new features (" . $counts{ num_new } . ")" );
    }
    else {
        ok 'skipping tests';
    }
};

subtest 'Clostridium botulinum tests' => sub {
    print "PURPOSE:\n";
    print
        "    1.  Clostridium botulinum is used as a test for selenocysteine-containing and for CRISPRS.\n";
    print "    2.  Both of these and one prophage should exist in this genome\n";

    my $genome_obj_name  = "Clostridium_botulinum";
    my $genome_gbff_name = "Clostridium_botulinum_310.gbff";
    my $genome_ref       = '15792/147849/1';

    my $params = {
        "input_genome"                      => $genome_ref,
        "call_selenoproteins"               => '1',
        "call_features_crispr"              => '1',
        "retain_old_anno_for_hypotheticals" => '1',
        "call_features_prophage_phispy"     => '1',
    };

    if ( $list{ clos } ) {
        my %counts = annotate( $genome_obj_name, $genome_gbff_name, $genome_ref, $params );
        ok( $counts{ seleno } == 0,
            "Number of Selenocysteine-containing genes = " . $counts{ seleno } . "\n" );
        ok( $counts{ crispr } == 0,
            "Number of features with CRISPR = " . $counts{ crispr } . "\n" );
        ok( $counts{ prophage } == 0,
            "Number of features with prophage = " . $counts{ prophage } . "\n" );
    }
    else {
        ok 'skipping tests';
    }

};

subtest 'Methanosarcina_acetivorans_C2A tests' => sub {
    print "PURPOSE:\n";
    print
        "    1.  Methanosarcina_acetivorans_C2A.gbff is used as a test for pyrrolysine-containing genes.\n";
    print "    2.  This should exist in this genome\n";

    my $genome_obj_name  = "Methanosarcina_acetivorans_C2A";
    my $genome_gbff_name = "Methanosarcina_acetivorans_C2A.gbff";
    my $genome_ref       = "15792/114497/2";

    my $params = {
        "input_genome"                  => $genome_ref,
        "call_pyrrolysoproteins"        => '1',
        "call_features_prophage_phispy" => '1',
    };
    if ( $list{ meth } ) {
        my %counts = annotate( $genome_obj_name, $genome_gbff_name, $genome_ref, $params );
        ok( $counts{ prophage } == 0,
            "Number of features with prophage = " . $counts{ prophage } . "\n" );
        ok( $counts{ pyrro } == 0,
            "Number of Pyrrolysine-containing genes = " . $counts{ pyrro } . "\n" );
    }
    else {
        ok 'skipping tests';
    }
};

subtest 'Shewanella spp. tests' => sub {
    print "PURPOSE:\n";
    print
        "    1.  Shewanella_sp._ANA-3 is used as a test for RAST repeats, kmer v2 and kmer v2.\n";
    print "    2.  Both of these should exist in this genome\n";
    print "    3.  In debug mode, using a genome reference in CI/prod.\n";

    my $genome_obj_name  = "Shewanella_sp._ANA-3";
    my $genome_gbff_name = "Shewanella_sp._ANA-3.gbff";
    my $genome_ref       = "15792/209997/1";

    my $params = {
        "input_genome"                     => $genome_ref,
        "call_features_repeat_region_SEED" => '1',
        "annotate_proteins_kmer_v2"        => '1',
        "kmer_v1_parameters"               => '1',
    };
    if ( $list{ shew } ) {
        my %counts = annotate( $genome_obj_name, $genome_gbff_name, $genome_ref, $params );
        ok( $counts{ num_diff } > 0,
            "Number of new annotation (" . $counts{ num_diff } . ")" );
        ok( $counts{ num_new } > 0, "Number of new features (" . $counts{ num_new } . ")" );
    }
    else {
        ok 'skipping tests';
    }
};

subtest 'Streptococcus pneumoniae tests' => sub {
    print "PURPOSE:\n";
    print
        "    1.  Streptococcus pneumoniae is a test for the strep-pneumoniae specific repeat\n";
    print "    2.  Streptococcus pneumoniae should also have SEED repeats\n";

    my $genome_obj_name  = "Streptococcus_pneumoniae_D39";
    my $genome_gbff_name = "Streptococcus_pneumoniae_D39.gbff";
    my $genome_ref       = "15792/103507/1";

    my $params = {
        "input_genome"                      => $genome_ref,
        "call_features_repeat_region_SEED"  => '1',
        "call_features_strep_pneumo_repeat" => '1',
        "retain_old_anno_for_hypotheticals" => '1',
    };
    if ( $list{ strepp } ) {
        my %counts = annotate( $genome_obj_name, $genome_gbff_name, $genome_ref, $params );
        ok( $counts{ repeat } > 0, "Number of repeats (" . $counts{ repeat } . ")" );
    }
    else {
        ok 'skipping tests';
    }
};

subtest 'Streptococcus_suis tests' => sub {
    print "PURPOSE:\n";
    print "    1.  Streptococcus suis is a test for the strep-suis specific repeat\n";
    print "    2.  This should exist in this genome\n";
    print "    3.  In debug mode, using a genome reference in CI/prod.\n";
    print
        "         Otherwise, load a genome in data dir. This takes more time but every user has access.\n";
    print "         For this reason, the test isn't for a specific number of changes\n\n";

    my $genome_obj_name  = "Streptococcus_suis";
    my $genome_gbff_name = "Streptococcus_suis.gbff";
    my $genome_ref       = '15792/60950/1';

    my $params = {
        "input_genome"                    => $genome_ref,
        "call_features_strep_suis_repeat" => '1',
    };
    if ( $list{ streps } ) {
        my %counts = annotate( $genome_obj_name, $genome_gbff_name, $genome_ref, $params );
        ok( $counts{ repeat } > 0, "Number of repeats (" . $counts{ repeat } . ")" );
    }
    else {
        ok 'skipping tests';
    }
};

sub annotate {
    my ( $genome_obj_name, $genome_gbff_name, $genome_ref, $params ) = @_;

    unless ( $DEBUG eq 'Y' ) {
        my $tmp_genome_obj;
        ( $tmp_genome_obj, $genome_ref ) = RASTTestUtils::prepare_gbff( $genome_gbff_name, $genome_obj_name );
        $params->{ input_genome } = $genome_ref;
    }

    my $num_func_in  = 0;
    my $num_func_out = 0;
    my %counts       = (
        'num_new'  => 0,
        'num_diff' => 0,
        'seleno'   => 0,
        'crispr'   => 0,
        'prophage' => 0,
        'pyrro'    => 0,
        'repeat'   => 0
    );

    my ( $orig_genome ) = RASTTestUtils::get_genome( $genome_ref );
    my ( $orig_funcs )  = RASTTestUtils::summarize_input( $orig_genome );

    print "number of input features = "
        . scalar @{ $orig_genome->{ features } } . "\n";
    print "number of input non-coding features = "
        . scalar @{ $orig_genome->{ non_coding_features } } . "\n";
    $num_func_in = scalar @{ $orig_genome->{ features } }
        + @{ $orig_genome->{ non_coding_features } };

    print( "######## Running RAST annotation ########\n" );
    my ( $genome_obj, $params ) = RASTTestUtils::submit_annotation( $genome_obj_name, $genome_ref, $params );

    print "number of returned features = "
        . scalar @{ $genome_obj->{ features } } . "\n";
    print "number of returned non-coding features = "
        . scalar @{ $genome_obj->{ non_coding_features } } . "\n";
    $num_func_out = scalar @{ $genome_obj->{ features } }
        + @{ $genome_obj->{ non_coding_features } };

    for my $ftr ( @{ $genome_obj->{ features } } ) {
        my $func = $ftr->{ function } // "";
        my $orig_func = $orig_funcs->{ $ftr->{ id } } // '';
        $counts{ num_diff }++ if $func ne $orig_func;

        if ( $ftr->{ ontology_terms } ) {

            #print Dumper $ftr->{ontology_terms}->{SSO};
            for ( keys %{ $ftr->{ ontology_terms }->{ SSO } } ) {
                $counts{ seleno }++ if ( $_ =~ /SSO:000009304/ );
                $counts{ pyrro }++  if ( $_ =~ /SSO:000009291/ );
            }
        }
    }

    for my $ftr ( @{ $genome_obj->{ non_coding_features } } ) {
        my $func = $ftr->{ function } // "";
        my $orig_func = $orig_funcs->{ $ftr->{ id } } // '';
        $counts{ num_diff }++ if $func ne $orig_func;
        if ( $ftr->{ function } && $ftr->{ type } =~ /crispr/ ) {
            $counts{ crispr }++;
        }
        elsif ( $ftr->{ function } && $ftr->{ function } =~ /phiSpy/ ) {
            $counts{ prophage }++;
        }
        if ( $ftr->{ type } =~ /repeat/ ) {
            $counts{ repeat }++;
        }
    }

    print "**** Number of features post-annotation = $num_func_out\n";
    $counts{ num_new } = $num_func_out - $num_func_in;

    print "Summary for $genome_obj_name\n";
    print "Number of features with changed function: " . $counts{ num_diff } . "\n";
    return ( %counts );
}

RASTTestUtils::clean_up();

done_testing;
