use Test::Most;
use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catfile splitpath);
use File::Copy;
use Carp qw(croak);

use installed_clients::WorkspaceClient;
use installed_clients::GenomeFileUtilClient;


use_ok "metag_utils";
use testRASTutil;


## global variables
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
my $fasta1 = 'data/short_one.fa';
my $gff1 = 'data/short_one.gff';
my $fasta2 = 'data/metag_test/59111.assembled.fna';
my $gff2 = 'data/metag_test/59111.assembled.gff';
my $fasta_shk = 'fasta_file.fa';
my $gff_shk = 'gff_file.gff';



sub generate_metagenome {
    my($ws, $metag_name, $fasta, $gff) = @_;
    my $fasta_path = catfile($scratch, $fasta_shk);
    my $gff_path = catfile($scratch, $gff_shk);

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


##-----------------Test Blocks--------------------##

subtest '_run_rast' => sub {
    my $ret_metag = generate_metagenome($ws, $out_name, $fasta1, $gff1);
    print Dumper($ret_metag);
    my $inputgenome = {
        features => []
    };

    my $input_obj_ref = $ret_metag->{metagenome_ref};
    my $input_fasta_file = catfile($metag_dir, 'prodigal_input.fasta');
    my $gff_filename = catfile($scratch, 'genome.gff');
    my ($fasta_contents, $gff_contents, $attr_delimiter) = ([], [], "=");
    
    $input_fasta_file = metag_utils::_write_fasta_from_metagenome(
                           $input_fasta_file, $input_obj_ref);
    $gff_filename = metag_utils::_write_gff_from_metagenome(
                       $gff_filename, $input_obj_ref);

    # fetch protein sequences and gene IDs from fasta and gff files
    $fasta_contents = metag_utils::_parse_fasta($input_fasta_file);
    ($gff_contents, $attr_delimiter) = metag_utils::_parse_gff(
                                          $gff_filename, $attr_delimiter);

    my $gene_seqs = metag_utils::_extract_cds_sequences_from_fasta(
                    $fasta_contents, $gff_contents);
    my $protein_seqs = metag_utils::_translate_gene_to_protein_sequences($gene_seqs);

    my %gene_id_index=();
    my $i=1;
    foreach my $gene (sort keys %$protein_seqs){
        push(@{$inputgenome->{features}},{
             id => "peg".$i,
                protein_translation => $protein_seqs->{$gene}
        });
        $gene_id_index{$i}=$gene;
        $i++;
    }

    my $rast_ret = metag_utils::_run_rast($inputgenome);
    print Dumper($rast_ret);
};

=begin
subtest 'rast_metagenome' => sub {
    my $input_params = {
        object_ref => $ret_metag->{metagenome_ref},
        output_metagenome_name => 'rast_metagenome',
        output_workspace => $ws
    };
 
    my $rast_mg = metag_utils::rast_metagenome($input_params);


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

    isnt(metag_utils::_run_prodigal(@p_cmd), 0, $cannot_open);

    $infile = $fasta1;
    @p_cmd = (
          'wrong_path/prodigal',
          '-i',
          $infile
    );

    isnt(metag_utils::_run_prodigal(@p_cmd), 0, "Can't exec wrong prodigal cmd.\n");

    @p_cmd = (
          $prodigal_cmd,
          '-i',
          $infile
    );
    lives_ok {
        $prd_ret = metag_utils::_run_prodigal(@p_cmd)
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
        $prd_ret = metag_utils::_run_prodigal(@p_cmd)
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
        $prd_ret = metag_utils::_run_prodigal(@p_cmd)
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
        metag_utils::_build_prodigal_cmd()
    } qr/$req/,
        '_build_prodigal_cmd dies missing required fasta or gbk or gff file.';

    # set default values correctly
    my @ret_cmd = metag_utils::_build_prodigal_cmd($p_input);
    isa_ok( \@ret_cmd, 'ARRAY' );
    #print Dumper(\@ret_cmd);
    cmp_deeply(\@ret_cmd, \@exp_cmd_default, $set_default_ok);

    # with specified output file name and the rest default
    my $out_file = 'out_dir_path/your_outfile';
    @ret_cmd = metag_utils::_build_prodigal_cmd($p_input, '', '', $out_file, '', '');
    is($ret_cmd[0], $exp_cmd_default[0], 'Prodigal command set correctly');
    is($ret_cmd[2], $p_input, 'input file name set correctly');
    is($ret_cmd[6], $out_file, 'output file name set correctly');

    # with specified output type
    my $out_type = 'sco';
    @ret_cmd = metag_utils::_build_prodigal_cmd($p_input, '', '', $out_file, $out_type, '');
    is($ret_cmd[4], $out_type, 'output type set correctly');

    # with specified translation file name
    my $trans = 'out_dir_path/your_translation';
    @ret_cmd = metag_utils::_build_prodigal_cmd($p_input, $trans, '', $out_file, $out_type, '');
    is($ret_cmd[8], $trans, 'translation file set correctly');

    # with specified nucleotide sequence file name
    my $nuc = 'out_dir_path/your_nuc';
    @ret_cmd = metag_utils::_build_prodigal_cmd($p_input, $trans, $nuc, $out_file, $out_type, '');
    is($ret_cmd[10], $nuc, 'nucleotide sequence file set correctly');

    # with specified procedure mode
    my $md = 'single';
    @ret_cmd = metag_utils::_build_prodigal_cmd($p_input, $trans, $nuc, $out_file, $out_type, $md);
    is($ret_cmd[16], $md, 'nucleotide sequence file set correctly');

    # set default values for '' and undef inputs
    @ret_cmd = metag_utils::_build_prodigal_cmd($p_input, undef, '', undef, undef, '');
    isa_ok( \@ret_cmd, 'ARRAY' );
    cmp_deeply(\@ret_cmd, \@exp_cmd_default, $set_default_ok);
};

subtest '_save_metagenome' => sub {
    my $req_params = "Missing required parameters for saving metagenome.\n";
    my $not_found = "file not found.\n";
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
                                      'data/nosuchfile.fna', $gff2)
    } qr/$not_found/,
      '_save_metagenome dies because fasta file not found';

    throws_ok {
        metag_utils::_save_metagenome($ws, $out_name,
                                      $fasta2, 'data/nosuchfile.gff')
    } qr/$not_found/,
      '_save_metagenome dies because GFF file not found';

    # metagenome saved successfully
    my $mymetag = {};
    lives_ok {
        $mymetag = metag_utils::_save_metagenome(
             $ws, $out_name, $fasta1, $gff1, $scratch)
    } '__save_metagenome run without errors on short_one.\n';
    ok (exists $mymetag->{metagenome_ref},
        "metagenome saved with metagenome_ref='$mymetag->{metagenome_ref}'");
    ok (exists $mymetag->{metagenome_info}, 'metagenome saved with metagenome_info');
    is ($mymetag->{metagenome_info}[1], $out_name, 'saved metagenome name is correct');
    is ($mymetag->{metagenome_info}[7], $ws, 'saved metagenome to the correct workspace');
    
    lives_ok {
        $mymetag = metag_utils::_save_metagenome(
             $ws, $out_name, $fasta2, $gff2, $scratch)
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
