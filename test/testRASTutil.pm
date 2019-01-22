use strict;
use warnings;
 
#use Data::Dumper;
#use Test::More;
#use Test::Exception;
#use Config::Simple;
#use Time::HiRes qw(time);
use Workspace::WorkspaceClient;
use JSON;
#use File::Copy;
#use AssemblyUtil::AssemblyUtilClient;
use GenomeFileUtil::GenomeFileUtilClient;
#use Storable qw(dclone);
#use Bio::KBase::kbaseenv;

local $| = 1;
my $token = $ENV{'KB_AUTH_TOKEN'};
my $config_file = $ENV{'KB_DEPLOYMENT_CONFIG'};
my $config = new Config::Simple($config_file)->get_block('RAST_SDK');
my $ws_url = $config->{"workspace-url"};
my $ws_name = undef;
my $ws_client = new Workspace::WorkspaceClient($ws_url,token => $token);
my $call_back_url = $ENV{ SDK_CALLBACK_URL };
#my $au = new AssemblyUtil::AssemblyUtilClient($call_back_url);
my $gfu = new GenomeFileUtil::GenomeFileUtilClient($call_back_url);

sub get_ws_name {
    if (!defined($ws_name)) {
        my $suffix = int(time * 1000);
        $ws_name = 'test_RAST_SDK_' . $suffix;
        $ws_client->create_workspace({workspace => $ws_name});
    }
    return $ws_name;
}

#--------------------------------------------
#	Call the RAST annotation Impl file
#	Create a JSON object of the parameters and submit
#--------------------------------------------
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
    return ($ret_obj);
}

#--------------------------------------------
#	Copy Genbank file to workspace and
#	use GFU to create the genome object and genome reference
#--------------------------------------------
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
		source => "Refseq"
    });
#	print "Genome Object Keys\n";
#	foreach my $key (keys(%$genome_obj))
#	{
#		print "\t$key ".ref($genome_obj->{$key})." \n";
#	}
	my $genome_ref = $genome_obj->{'genome_ref'};
	return ($genome_obj, $genome_ref);
}
 
sub get_and_prep {
	my $genome_ref = shift;
	my $orig_genome = $ws_client->get_objects([{ref=>$genome_ref}])->[0]->{data};

    if (defined($orig_genome->{taxon_ref})) {
        $orig_genome->{taxon_ref} = "ReferenceTaxons/unknown_taxon" ;
    }
    if (defined($orig_genome->{genbank_handle_ref})) {
        delete $orig_genome->{genbank_handle_ref};
    }

	my $orig_funcs = {};
	my %types;

	if (! exists $orig_genome->{features}) {
		$orig_genome->{features} = [];
	}
	if (! exists $orig_genome->{non_coding_features}) {
		$orig_genome->{non_coding_features} = [];
	}

	for (my $i=0; $i < scalar @{$orig_genome->{features}}; $i++) {
		my $ftr = $orig_genome->{features}->[$i];
		if (defined($ftr->{functions}) && scalar @{$ftr->{functions}} > 0){
					$ftr->{function} = join("; ", @{$ftr->{functions}});
		}
        my $func      = defined($ftr->{function}) ? $ftr->{function} : "";
        $orig_funcs->{$ftr->{id}} = $func;
		if (not defined($ftr->{type})) {
           if (defined($ftr->{protein_translation})) {
              $ftr->{type} = "gene";
           } else {
              $ftr->{type} = "other";
           }
		}
		if ($ftr->{type} eq "gene" and not(defined($ftr->{protein_translation}))) {
			if (lc($func) =~ /ribosomal/) {
				$ftr->{type} = 'rRNA';
			}
			elsif (lc($func) =~ /trna/) {
				$ftr->{type} = 'tRNA';
			} else {
				$ftr->{type} = 'other non-coding';
			}
		}
		%types = &count_types($ftr->{type},%types);
    }
	for (my $i=0; $i < scalar @{$orig_genome->{non_coding_features}}; $i++) {
		my $ftr = $orig_genome->{non_coding_features}->[$i];
		if (defined($ftr->{functions}) && scalar @{$ftr->{functions}} > 0){
					$ftr->{function} = join("; ", @{$ftr->{functions}});
		}
        my $func      = defined($ftr->{function}) ? $ftr->{function} : "";
        $orig_funcs->{$ftr->{id}} = $func;
		if (not defined($ftr->{type})) {
           if (defined($ftr->{protein_translation})) {
              $ftr->{type} = "gene";
           } else {
              $ftr->{type} = "other";
           }
		}
		%types = &count_types($ftr->{type},%types);
    }
	return ($orig_genome,$orig_funcs);

}

sub count_types {
	my ($feat_type,%types) = @_;

	if (exists $types{$feat_type}) {
		$types{$feat_type} ++;
	} else {
		$types{$feat_type} = 1;
	}
	return %types;
}

sub report_types {
	my %types = @_;
	foreach my $key (keys %types) 	{
		print "\t$key = $types{$key}\n";
	}
}

sub submit_annotation {
	my ($genome_obj_name, $genome_ref) = @_;

    my ($ret,$params) = reannotate_genome($genome_obj_name, $genome_ref);
	my $report_ref = $ret->{report_ref};
    my $report_obj = $ws_client->get_objects([{ref=>$report_ref}])->[0]->{data};
	my $report_text = $report_obj->{direct_html};
    print "\nReport: " . $report_text . "\n\n";
    $genome_ref = get_ws_name() . "/" . $genome_obj_name ;
    my $genome_obj = $ws_client->get_objects([{ref=>$genome_ref}])->[0]->{data};

#	print "Genome Object Keys\n";
#	foreach my $key (keys(%$genome_obj))
#	{
#		print "\t$key ".ref($genome_obj->{$key})." \n";
#	}
	if (! exists $genome_obj->{features}) {
		$genome_obj->{features} = [];
	}
	if (! exists $genome_obj->{non_coding_features}) {
		$genome_obj->{non_coding_features} = [];
	}
	return ($genome_obj,$params);
}

1;
