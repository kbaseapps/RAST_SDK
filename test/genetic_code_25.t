use Test::Most;
use RASTTestUtils;
use Data::Dumper::Concise;
use feature qw( say );

say qq(PURPOSE:
    1.  Test annotation of genetic code 25 for Mircea Podar.
    2.  Test the Prodigal is using the -meta when the domain is 'Unknown'
);

my $ws_client     = RASTTestUtils::get_ws_client();
my $ws_name       = RASTTestUtils::get_ws_name();

subtest 'Running RAST annotation' => sub {

    my $assembly_obj_name = "GCA_000350285.1_OR1_genomic.fna";
    my $assembly_ref      = RASTTestUtils::prepare_assembly( $assembly_obj_name );
    my $genome_obj_name   = 'SR1_bacterium_MGEHA';
    my $genome_set_name   = "New_GenomeSet";

    my $params = {
        "input_contigset"            => $assembly_obj_name,
        "scientific_name"            => 'candidate division SR1 bacterium MGEHA',
        "domain"                     => 'U',
        "genetic_code"               => '25',
        "call_features_CDS_prodigal" => '1',
    };

    $params = RASTTestUtils::set_params( $genome_obj_name, $params );

    RASTTestUtils::make_impl_call( "RAST_SDK.annotate_genome", $params );
    my $genome_ref = $ws_name . "/" . $genome_obj_name;
    my $genome_obj = $ws_client->get_objects( [ { ref => $genome_ref } ] )->[ 0 ]{ data };

    say "OUTPUT OBJECT DOMAIN = $genome_obj->{domain}\n"
    . "OUTPUT OBJECT G_CODE = $genome_obj->{genetic_code}\n"
    . "OUTPUT SCIENTIFIC NAME = $genome_obj->{scientific_name}\n";

    ok $genome_obj->{ features } && @{ $genome_obj->{ features } },
        "Features array is present and contains more than zero items"
        or diag explain $genome_obj;

};

RASTTestUtils::clean_up();

done_testing;
