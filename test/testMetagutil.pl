use Test::Most;
use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catfile splitpath);
use File::Copy;
use Carp qw(croak);

use installed_clients::WorkspaceClient;
use installed_clients::GenomeFileUtilClient;


use_ok "metag_utils";

use testRASTutil;

my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('RAST_SDK');
my $ws_url = $config->{"workspace-url"};
my $ws = undef;
my $ws_client = new installed_clients::WorkspaceClient($ws_url,token => $token);
my $call_back_url = $ENV{ SDK_CALLBACK_URL };


$ws = get_ws_name();
my $out_name = 'annotated_metag';
my $scratch = metag_utils::_create_metag_dir();

sub generate_metagenome {
    my($ws, $metag_name, $fasta, $gff) = @_;
    my $fasta_path = catfile($scratch, "fasta_file.fa");
    my $gff_path = catfile($scratch, "gff_file.gff");

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

subtest '_save_metagenome' => sub {
    my $req_params = "Missing required parameters for saving metagenome.\n";
    my $not_found = "file not found.\n";
    my $gff_not_found = "GFF file not found.\n";
    my $req1 = "Both 'output_workspace' and 'output_metagenome_name' are required.\n";
    my $req2 = "Both 'fasta_file' and 'gff_file' are required.\n";

    throws_ok {
        metag_utils::_save_metagenome()
    } qr/$req_params/,
        '_save_metagenome dies without params';

    throws_ok {
        metag_utils::_save_metagenome($ws, $out_name, undef, 'abc')
    } qr/$req_params/,
        '_save_metagenome dies with no fasta file';

    throws_ok {
        metag_utils::_save_metagenome($ws, $out_name, undef)
    } qr/$req_params/,
        '_save_metagenome dies with no gff file';

    throws_ok {
        metag_utils::_save_metagenome($ws, $out_name, 'abc', 'def')
    } qr/$not_found/,
        '_save_metagenome dies because input file not found.';

    throws_ok {
        metag_utils::_save_metagenome($ws, $out_name,
             'data/metag_test/nosuchfile.fna',
             'data/metag_test/59111.assembled.gff')
    } qr/$not_found/,
      '_save_metagenome dies because fasta file not found';

    throws_ok {
        metag_utils::_save_metagenome($ws, $out_name,
             'data/metag_test/59111.assembled.fna',
             'data/metag_test/nosuchfile.gff')
    } qr/$not_found/,
      '_save_metagenome dies because GFF file not found';

    # metagenome saved successfully
    my $mymetag = {};
    lives_ok {
        $mymetag = metag_utils::_save_metagenome(
             $ws, $out_name,
             'data/short_one.fa',
             'data/short_one.gff',
             $scratch)
    } '__save_metagenome run without errors on short_one.\n';
    
    my $save_ok = '_save_metagenome runs ok with given parameters';
    # cmp_deeply($mymetag1->{metagenome_ref}, $exp_metag, $save_ok);

    lives_ok {
        $mymetag = metag_utils::_save_metagenome(
             $ws, $out_name,
             'data/metag_test/59111.assembled.fna',
             'data/metag_test/59111.assembled.gff',
             $scratch)
    } '__save_metagenome run without errors on 59111.assembled.\n';
    # print Dumper($mymetag);
    
    my $ret_metag = generate_metagenome(
             $ws, $out_name,
             'data/short_one.fa',
             'data/short_one.gff');

    print Dumper($ret_metag);
    # cmp_deeply($ret_metag, $exp_metag, $save_ok);

};


=begin
subtest '_check_annotation_params' => sub {
    my $obj = '1234/56/7';

    my $missing_params = "Missing required parameters for annotating metagenome.\n";
    my $req1 = "'output_workspace' is required for running rast_metagenome.\n";
    my $invald1 = "Invalid workspace name:";
    my $req2 = "'object_ref' is required for running rast_metagenome.\n";
    my $invald2 = "Invalid workspace object reference:";

    throws_ok {
        metag_utils::_check_annotation_params()
    } qr/$missing_params/,
        '_check_annotation_params dies without params';

    throws_ok {
        metag_utils::_check_annotation_params( {} )
    } qr/$missing_params/,
        '_check_annotation_params dies with an empty hashref';

    throws_ok {
        my $p = {output_workspace => $ws,
                 output_metagenome_name => $out_name};
        print "input parameter=\n". Dumper($p);
        metag_utils::_check_annotation_params($p)
    } qr/$req2/,
        '_check_annotation_params dies with no object_ref';

    throws_ok {
        my $p = {output_workspace => $ws,
                 output_metagenome_name => $out_name,
                 object_ref => ''};
        print "input parameter=\n". Dumper($p);
        metag_utils::_check_annotation_params($p)
    } qr/$req2/,
        '_check_annotation_params dies with blank object_ref';

    throws_ok {
        metag_utils::_check_annotation_params(
            {object_ref => $obj,
             output_metagenome_name => $out_name})
    } qr/$req1/,
        '_check_annotation_params dies with no outpout_workspace';

    throws_ok {
        metag_utils::_check_annotation_params(
            {workspace => $ws,
             output_metagenome_name => $out_name,
             obect_ref => $obj})
    } qr/$req1/,
        '_check_annotation_params dies with wrong workspace key';

    throws_ok {
        metag_utils::_check_annotation_params(
            {output_workspace => '',
             output_metagenome_name => $out_name,
             obect_ref => $obj})
    } qr/$req1/,
        '_check_annotation_params dies with blank workspace name';

    throws_ok {
        metag_utils::_check_annotation_params(
            {output_workspace => $ws,
             output_metagenome_name => $out_name,
             object_ref => 'abc/1/2'})
    } qr/$invald2/,
      '_check_annotation_params dies because of invalid object reference format';

    throws_ok {
        metag_utils::_check_annotation_params(
            {output_workspace => 'ab:c',
             output_metagenome_name => $out_name,
             object_ref => '456/1/2'})
    } qr/$invald1/,
      '_check_annotation_params dies because of invalid workspace name format';

    # _check_annotation_params passed
    my $expected = {
             output_workspace => $ws,
             output_metagenome_name => 'rast_annotated_metagenome',
             object_ref => '456/1/2'};
    my $set_default_ok = '_check_annotation_params sets the default value for output_metagenome_name.';

    my $ret = metag_utils::_check_annotation_params(
            {output_workspace => $ws,
             output_metagenome_name => undef,
             object_ref => '456/1/2'});
    cmp_deeply($ret, $expected, 'When undefined, '.$set_default_ok);

    $ret = metag_utils::_check_annotation_params(
            {output_workspace => $ws,
             output_metagenome_name => '',
             object_ref => '456/1/2'});
    cmp_deeply($ret, $expected, 'When blank, '.$set_default_ok);

};

=cut

done_testing();

=begin
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
=cut
