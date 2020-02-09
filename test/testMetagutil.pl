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

use RAST_SDK::RAST_SDKImpl;
use_ok "metag_utils";
use testRASTutil;

local $| = 1;

## global variables
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('RAST_SDK');
my $auth_token = Bio::KBase::AuthToken->new(
        token => $token, ignore_authrc => 1, auth_svc=>$config->{'auth-service-url'});
my $ws_url = $config->{"workspace-url"};
my $ws = undef;
my $ws_client = new installed_clients::WorkspaceClient($ws_url,token => $token);
my $call_back_url = $ENV{ SDK_CALLBACK_URL };

$ws = get_ws_name();
my $out_name = 'annotated_metag';
my $fasta1 = 'data/short_one.fa';
my $gff1 = 'data/short_one.gff';
my $fasta2 = 'data/metag_test/59111.assembled.fna';
my $gff2 = 'data/metag_test/59111.assembled.gff';
my $fasta_scrt = 'fasta_file.fa';
my $gff_scrt = 'gff_file.gff';

my $rast_impl = new RAST_SDK::RAST_SDKImpl();

my $ctx = LocalCallContext->new($token, $auth_token->user_id);
my $mgutil = new metag_utils($config, $ctx);

my $scratch = $config->{'scratch'}; #'/kb/module/work/tmp';
my $rast_dir = $mgutil->_create_metag_dir($scratch);


sub generate_metagenome {
    my($ws, $metag_name, $fasta, $gff) = @_;
    my $fasta_path = catfile($rast_dir, $fasta_scrt);
    my $gff_path = catfile($rast_dir, $gff_scrt);

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

## global objects/variables for multiple subtests
my $ret_metag = generate_metagenome($ws, $out_name, $fasta1, $gff1);
print Dumper($ret_metag);
my $input_obj_ref = $ret_metag->{genome_ref};

my $input_fasta_file = catfile($rast_dir, 'prodigal_input.fasta');
my $gff_filename = catfile($rast_dir, 'genome.gff');
my ($fasta_contents, $gff_contents, $attr_delimiter) = ([], [], "=");

$input_fasta_file = $mgutil->_write_fasta_from_metagenome(
		    $input_fasta_file, $input_obj_ref);
$gff_filename = $mgutil->_write_gff_from_metagenome(
	        $gff_filename, $input_obj_ref);
# fetch protein sequences and gene IDs from fasta and gff files
$fasta_contents = $mgutil->_parse_fasta($input_fasta_file);
($gff_contents, $attr_delimiter) = $mgutil->_parse_gff(
				  $gff_filename, $attr_delimiter);

my $gene_seqs = $mgutil->_extract_cds_sequences_from_fasta(
                    $fasta_contents, $gff_contents);
print "**********Gene sequences************\n" . Dumper($gene_seqs);
my $protein_seqs = $mgutil->_translate_gene_to_protein_sequences($gene_seqs);
print "**********Protein sequences************\n" . Dumper($protein_seqs);



##-----------------Test Blocks--------------------##

my $ecoli_gff = 'data/metag_test/ecoli_out.gff';
my $ecoli_sco = 'data/metag_test/ecoli_out.sco';
my $trans_file = 'data/metag_test/translationfile';
my %trans_tab;
my $sco_tab = [];

subtest '_parse_translation' => sub {
    my $trans_path = catfile($rast_dir, 'trans_scrt');
    copy($trans_file, $trans_path) || croak "Copy file failed: $!\n";

    %trans_tab = $mgutil->_parse_translation($trans_path);
    ok( keys %trans_tab , 'Prodigal translation parsing returns result.');
};

subtest '_parse_sco' => sub {
    $sco_tab = $mgutil->_parse_sco($ecoli_sco, %trans_tab);
    ok( @{$sco_tab} >0, 'Prodigal SCO parsing returns result.');
};

subtest '_parse_gff' => sub {
    my ($gff_contents, $attr_delimiter) = $mgutil->_parse_gff(
                                              $ecoli_gff, '=');
    ok( @{$gff_contents} >0, 'Parsing GFF returns result.');
};

subtest '_parse_prodigal_results' => sub {
    my $prd_out_path = catfile($rast_dir, 'prd_out');
    my $trans_path = catfile($rast_dir, 'trans_scrt');
    copy($trans_file, $trans_path) || croak "Copy file failed: $!\n";

    # Prodigal generate a GFF output file
    my $out_type = 'gff';
    copy($ecoli_gff, $prd_out_path) || croak "Copy file failed: $!\n";
    my $prd_results = $mgutil->_parse_prodigal_results(
                          $trans_path, $prd_out_path, $out_type);
    ok( @{$prd_results} >0, 'Prodigal GFF parsing returns result.');

    # Prodigal generate an SCO output file
    $out_type = 'sco';
    copy($ecoli_sco, $prd_out_path) || croak "Copy file failed: $!\n";
    $prd_results = $mgutil->_parse_prodigal_results(
                          $trans_path, $prd_out_path, $out_type);
    ok( @{$prd_results} >0, 'Prodigal SCO parsing returns result.');

};

subtest '_run_rast' => sub {
    my $inputgenome = {
        features => []
    };
    my $i=1;
    foreach my $gene (sort keys %$protein_seqs){
        push(@{$inputgenome->{features}},{
             id => "peg".$i,
             protein_translation => $protein_seqs->{$gene}
        });
        $i++;
    }

    my $rast_ret = $mgutil->_run_rast($inputgenome);
    isnt($rast_ret, undef, 'RAST run_pipeline call returns results');

};

subtest 'rast_metagenome' => sub {
    my $parms = {
        object_ref => $ret_metag->{metagenome_ref},
        output_metagenome_name => 'rasted_metagenome',
        output_workspace => $ws
    };
    my $rast_mg_ref;

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

    $parms = {
        "object_ref" => "37798/7/1",
        "output_metagenome_name" => "rasted_shortOne_appdev",
        "output_workspace" => $ws
    };
    $rast_mg_ref = $mgutil->rast_metagenome($parms);
    print "rast_metagenome returns: $rast_mg_ref" if defined($rast_mg_ref);
    ok (($rast_mg_ref !~ m/[^\\w\\|._-]/), 'rast_metagenome returns an INVALID ref');
};

subtest 'annotate_metagenome' => sub {
    my $parms = {
        "object_ref" => "37798/7/1",
        "output_metagenome_name" => "rasted_shortOne_appdev",
        "output_workspace" => $ws
    };
    my $rast_ann = $rast_impl->annotate_metagenome($parms);
    print Dumper($rast_ann);

};


=begin
subtest '_write_fasta_from_metagenome' => sub {
    my $fa_test1 = catfile($rast_dir, 'test1.fasta');
    $fa_test1 = $mgutil->_write_fasta_from_metagenome(
                    $fa_test1, $input_obj_ref);

    ok((-e $fa_test1), 'fasta file created');
    ok((-s $fa_test1), 'fasta file has data');
    # ok(compare($fa_test1, $fasta1) == 0, 'fasta file written correctly');
};

subtest '_write_gff_from_metagenome' => sub {
    my $gff_test1 = catfile($rast_dir, 'test1.gff');
    $gff_test1 = $mgutil->_write_gff_from_metagenome(
		   $gff_test1, $input_obj_ref);

    ok((-e $gff_test1), 'gff file created');
    ok((-s $gff_test1), 'gff file has data');
    # ok(compare($gff_test1, $gff1) == 0, 'GFF file written correctly');

};
=cut
=begin
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
    isnt($prd_ret, 0, 'Prodigal returned: '.$seq_too_short.'\n');

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
    is($prd_ret, 0, 'Prodigal runs ok with -p meta option.\n');

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

    is($prd_ret, 0, '_run_prodigal successfully returned.\n');

};


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

subtest '_save_metagenome' => sub {
    my $req_params = "Missing required parameters for saving metagenome.\n";
    my $not_found = "file not found.\n";
    my $req1 = "Both 'output_workspace' and 'output_metagenome_name' are required.\n";
    my $req2 = "Both 'fasta_file' and 'gff_file' are required.\n";

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
        $mgutil->_save_metagenome($ws, $out_name, 'abc', 'def')
    } qr/$not_found/,
        '_save_metagenome dies because input file not found.';

    throws_ok {
        $mgutil->_save_metagenome($ws, $out_name,
                                      'data/nosuchfile.fna', $gff2)
    } qr/$not_found/,
      '_save_metagenome dies because fasta file not found';

    throws_ok {
        $mgutil->_save_metagenome($ws, $out_name,
                                      $fasta2, 'data/nosuchfile.gff')
    } qr/$not_found/,
      '_save_metagenome dies because GFF file not found';

    # metagenome saved successfully
    my $mymetag = {};
    lives_ok {
        $mymetag = $mgutil->_save_metagenome(
             $ws, $out_name, $fasta1, $gff1, $rast_dir)
    } '__save_metagenome run without errors on short_one.\n';
    ok (exists $mymetag->{metagenome_ref},
        "metagenome saved with metagenome_ref='$mymetag->{metagenome_ref}'");
    ok (exists $mymetag->{metagenome_info}, 'metagenome saved with metagenome_info');
    is ($mymetag->{metagenome_info}[1], $out_name, 'saved metagenome name is correct');
    is ($mymetag->{metagenome_info}[7], $ws, 'saved metagenome to the correct workspace');
    
    lives_ok {
        $mymetag = $mgutil->_save_metagenome(
             $ws, $out_name, $fasta2, $gff2, $rast_dir)
    } '_save_metagenome runs without errors on 59111.assembled.\n';
    ok (exists $mymetag->{metagenome_ref},
        "metagenome saved with metagenome_ref='$mymetag->{metagenome_ref}'");
    ok (exists $mymetag->{metagenome_info}, 'metagenome saved with metagenome_info');
    is ($mymetag->{metagenome_info}[1], $out_name, 'saved metagenome name is correct');
    is ($mymetag->{metagenome_info}[7], $ws, 'saved metagenome to the correct workspace');

};

subtest '_check_annotation_params' => sub {
    my $obj = '1234/56/7';

    my $missing_params = "Missing required parameters for annotating metagenome.\n";
    my $req1 = "'output_workspace' is required for running rast_metagenome.\n";
    my $req2 = "'object_ref' is required for running rast_metagenome.\n";

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
                 output_metagenome_name => $out_name};
        print "input parameter=\n". Dumper($p);
        $mgutil->_check_annotation_params($p)
    } qr/$req2/,
        '_check_annotation_params dies with no object_ref';

    throws_ok {
        my $p = {output_workspace => $ws,
                 output_metagenome_name => $out_name,
                 object_ref => ''};
        print "input parameter=\n". Dumper($p);
        $mgutil->_check_annotation_params($p)
    } qr/$req2/,
        '_check_annotation_params dies with blank object_ref';

    throws_ok {
        $mgutil->_check_annotation_params(
            {object_ref => $obj,
             output_metagenome_name => $out_name})
    } qr/$req1/,
        '_check_annotation_params dies with no outpout_workspace';

    throws_ok {
        $mgutil->_check_annotation_params(
            {workspace => $ws,
             output_metagenome_name => $out_name,
             obect_ref => $obj})
    } qr/$req1/,
        '_check_annotation_params dies with wrong workspace key';

    throws_ok {
        $mgutil->_check_annotation_params(
            {output_workspace => '',
             output_metagenome_name => $out_name,
             obect_ref => $obj})
    } qr/$req1/,
        '_check_annotation_params dies with blank workspace name';

    lives_ok {
        $mgutil->_check_annotation_params(
            {output_workspace => $ws,
             output_metagenome_name => $out_name,
             object_ref => 'abc/1/2'})
    } '_check_annotation_params object_ref check ok';

    lives_ok {
        $mgutil->_check_annotation_params(
            {output_workspace => 'ab:c',
             output_metagenome_name => $out_name,
             object_ref => '456/1/2'})
    } '_check_annotation_params workspace name check ok';

    # _check_annotation_params passed
    my $expected = {
             output_workspace => $ws,
             output_metagenome_name => 'rast_annotated_metagenome',
             object_ref => '456/1/2'};
    my $set_default_ok = '_check_annotation_params sets the default value for output_metagenome_name.';

    my $ret = $mgutil->_check_annotation_params(
            {output_workspace => $ws,
             output_metagenome_name => undef,
             object_ref => '456/1/2'});
    cmp_deeply($ret, $expected, 'When undefined, '.$set_default_ok);

    $ret = $mgutil->_check_annotation_params(
            {output_workspace => $ws,
             output_metagenome_name => '',
             object_ref => '456/1/2'});
    cmp_deeply($ret, $expected, 'When blank, '.$set_default_ok);

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

{
     package LocalCallContext;
     use strict;
     sub new {
         my($class,$token,$user) = @_;
         my $self = {
             token => $token,
             user_id => $user
         };
         return bless $self, $class;
     }
     sub user_id {
         my($self) = @_;
         return $self->{user_id};
     }
     sub token {
         my($self) = @_;
         return $self->{token};
     }
     sub provenance {
         my($self) = @_;
         return [{'service' => 'RAST_SDK', 'method' => 'please_never_use_it_in_production', 'method_params' => []}];
     }
     sub authenticated {
         return 1;
     }
     sub log_debug {
         my($self,$msg) = @_;
         print STDERR $msg."\n";
     }
     sub log_info {
         my($self,$msg) = @_;
         print STDERR $msg."\n";
     }
     sub method {
         my($self) = @_;
         return "TEST_METHOD";
     }
 }
