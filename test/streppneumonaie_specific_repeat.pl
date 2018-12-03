use strict;
use Data::Dumper;
use Test::More;
use Test::Exception;
use Config::Simple;
use Time::HiRes qw(time);
use Workspace::WorkspaceClient;
use JSON;
use File::Copy;
use AssemblyUtil::AssemblyUtilClient;
use GenomeFileUtil::GenomeFileUtilClient;
use Storable qw(dclone);
use Bio::KBase::kbaseenv;

my $purpose = "Streptococcus pneumoniae is a test for the strep-pneumoniae specific repeat\n";
$purpose   .= "Streptococcus pneumoniae should also have SEED repeats\n";
$purpose   .= "This should exist in this genome\n";

local $| = 1;
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('RAST_SDK');
my $ws_url = $config->{"workspace-url"};
my $ws_name = undef;
my $ws_client = new Workspace::WorkspaceClient($ws_url,token => $token);
my $call_back_url = $ENV{ SDK_CALLBACK_URL };
my $au = new AssemblyUtil::AssemblyUtilClient($call_back_url);
my $gfu = new GenomeFileUtil::GenomeFileUtilClient($call_back_url);

sub get_ws_name {
    if (!defined($ws_name)) {
        my $suffix = int(time * 1000);
        $ws_name = 'test_RAST_SDK_' . $suffix;
        $ws_client->create_workspace({workspace => $ws_name});
    }
    return $ws_name;
}

sub make_impl_call {
    my ($method, $params) = @_;
    # Prepare json file input
    my $input = {
        method => $method,
        params => $params,
        version => "1.1",
        id => "1"
    };
    open my $fh, ">", "/kb/module/work/input.json";
    print $fh encode_json($input);
    close $fh;
    my $output_path = "/kb/module/work/output.json";
    if (-e $output_path) {
        unlink($output_path);
    }
    # Run run_async.sh
    system("sh", "/kb/module/scripts/run_async.sh");
    # Load json file with output
    unless (-e $output_path) {
        die "Output file of JSON-RPC call (in CLI mode) wasn't created";
    }
    open my $fh2, "<", $output_path;
    my $output_json = <$fh2>;
    close $fh2;
    my $json = JSON->new;
    my $output = $json->decode($output_json);
    if (defined($output->{error})) {
        die encode_json($output->{error});
    }
    my $ret_obj = $output->{result}->[0];
    return $ret_obj;
}

sub prepare_gbff {
    my($genome_gbff_name,$genome_obj_name) = @_;
	my $genome_gbff_path = "/kb/module/test/data/$genome_gbff_name";
    my $temp_path = "/kb/module/work/tmp/$genome_gbff_name";
    copy $genome_gbff_path, $temp_path;

	print "***** Loading and Saving the Genome Object *****\n";  
    my $genome_obj = $gfu->genbank_to_genome({
        workspace_name => get_ws_name(),
        genome_name => $genome_obj_name,
        file => {"path" => $temp_path},
		source => "Streptococcus_pneumoniae_D39"
    });
#	print "Genome Object Keys\n";
#	foreach my $key (keys(%$genome_obj))
#	{
#		print "\t$key ".ref($genome_bj->{$key})." \n";
#	}
	my $genome_ref = $genome_obj->{'genome_ref'};
	return ($genome_obj, $genome_ref);
}

sub reannotate_genome {
    my($genome_obj_name, $genome_upa) = @_;
    my $params={"input_genome"=>$genome_upa,
             "call_features_rRNA_SEED"=>'0',
             "call_features_tRNA_trnascan"=>'0',
             "call_selenoproteins"=>'0',
             "call_pyrrolysoproteins"=>'0',
             "call_features_repeat_region_SEED"=>'1',
             "call_features_insertion_sequences"=>'0',
             "call_features_strep_suis_repeat"=>'0',
             "call_features_strep_pneumo_repeat"=>'1',
             "call_features_crispr"=>'0',
             "call_features_CDS_glimmer3"=>'0',
             "call_features_CDS_prodigal"=>'0',
             "annotate_proteins_kmer_v2"=>'0',
             "kmer_v1_parameters"=>'0',
             "annotate_proteins_similarity"=>'0',
             "retain_old_anno_for_hypotheticals"=>'1',
             "resolve_overlapping_features"=>'0',
             "call_features_prophage_phispy"=>'0',
             "output_genome"=>$genome_obj_name,
             "workspace"=>get_ws_name()
           };
    return make_impl_call("RAST_SDK.annotate_genome", $params);
}

my $diff_count = 0;
my $num_func = 0;
my $num_func_out = 0;
my $repeat = 0;
lives_ok {
    my $genome_obj_name = "genome.1";
    my $genome_gbff_name = "Streptococcus_pneumoniae_D39.gbff";
    my ($tmp_genome_obj, $genome_ref) = prepare_gbff($genome_gbff_name,$genome_obj_name);
	my $orig_genome = $ws_client->get_objects([{ref=>$genome_ref}])->[0]->{data};

    if (defined($orig_genome->{taxon_ref})) {
        $orig_genome->{taxon_ref} = "ReferenceTaxons/unknown_taxon" ;
    }
    if (defined($orig_genome->{genbank_handle_ref})) {
        delete $orig_genome->{genbank_handle_ref};
    }

	my $orig_funcs = {};
	my %types;
	for (my $i=0; $i < scalar @{$orig_genome->{features}}; $i++) {
		my $ftr = $orig_genome->{features}->[$i];
		my $func = $ftr->{function};
		if (not defined($func)) {
			$func = "";
		}
		if (not defined($ftr->{type})) {
			$ftr->{type} = "other";
		}
		if ($ftr->{type} eq "gene" and not(defined($ftr->{protein_translation}))) {
			if (lc($ftr->{function}) =~ /ribosomal/) {
				$ftr->{type} = 'rRNA';
			} elsif (lc($ftr->{function}) =~ /trna/) {
				$ftr->{type} = 'tRNA';
			} else {
				$ftr->{type} = 'other non-coding';
			}
		}
		if (exists $types{$ftr->{type}}) {
			$types{$ftr->{type}} ++;
		} else {
			$types{$ftr->{type}} = 1;
		}
		$num_func++;
        $orig_funcs->{$ftr->{id}} = $func;
    }

#	print "**** Number of input features = $num_func\n";
#	foreach my $key (keys %types)
#	{
#		print "\t$key = $types{$key}\n";
#	}

	$genome_obj_name .= "_new";
    print("$purpose\nRunning RAST annotation\n");
    my $ret = reannotate_genome($genome_obj_name, $genome_ref);
	my $report_ref = $ret->{report_ref};
    my $report_obj = $ws_client->get_objects([{ref=>$report_ref}])->[0]->{data};
	my $report_text = $report_obj->{direct_html};
    print "\n$purpose\n\nReport: " . $report_text . "\n\n";
    my $genome_ref = get_ws_name() . "/" . $genome_obj_name ;
    my $genome_obj = $ws_client->get_objects([{ref=>$genome_ref}])->[0]->{data};

#	print "Genome Object Keys\n";
#	foreach my $key (keys(%$genome_obj))
#	{
#		print "\t$key ".ref($genome_obj->{$key})." \n";
#	}


	print "number of features = ".scalar  @{$genome_obj->{features}}."\n";
	print "number of non-coding features = ".scalar  @{$genome_obj->{non_coding_features}}."\n";

    %types = ();
    for (my $i=0; $i < @{$genome_obj->{features}}; $i++) {
        my $ftr = $genome_obj->{features}->[$i];
        my $func = $ftr->{function};
        if (not defined($func)) {
            $func = "";
        }
        my $orig_func = $orig_funcs->{$ftr->{id}};
        if (not($func eq $orig_func)) {
            $diff_count++;
        }
		$types{$ftr->{type}} += 1;
		$num_func_out++;
    }
    for (my $i=0; $i < @{$genome_obj->{non_coding_features}}; $i++) {
        my $ftr = $genome_obj->{non_coding_features}->[$i];
        my $func = $ftr->{non_coding_function};
        if (not defined($func)) {
            $func = "";
        }
        my $orig_func = $orig_funcs->{$ftr->{id}};
        if (not defined($orig_func)) {
            $func = "";
        }
        if (not($func eq $orig_func)) {
            $diff_count++;
        }
		$types{$ftr->{type}} += 1;
		$num_func_out++;
		if ($ftr->{type} =~ /repeat/) {
#			print Dumper $ftr->{ontology_terms};
#			print "Found a NON-coding repeat $ftr->{id}\n";
			$repeat++;
		}
    }

	print "Output Feature Types\n";
	foreach my $key (keys %types)
	{
		print "\t$key = $types{$key}\n";
	}

	print "**** Number of features post-annotation = $num_func_out\n";
} "Pipeline Runs";
print "Number of features with changed function: " . $diff_count . "\n";
ok($repeat > 0, "Number of repeats (" . $repeat . ")");
done_testing(2);

my $err = undef;
if ($@) {
    $err = $@;
}
eval {
    if (defined($ws_name)) {
        $ws_client->delete_workspace({workspace => $ws_name});
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


