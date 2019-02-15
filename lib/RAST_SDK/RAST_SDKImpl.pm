package RAST_SDK::RAST_SDKImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = '0.0.16';
our $GIT_URL = 'https://github.com/landml/RAST_SDK';
our $GIT_COMMIT_HASH = '152bb0b02d601d031316b04ed4e817f978a049db';

=head1 NAME

RAST_SDK

=head1 DESCRIPTION

The SDK version of the KBaase Genome Annotation Service.
This wraps genome_annotation which is based off of the SEED annotations.

=cut

#BEGIN_HEADER
use Bio::KBase::AuthToken;
use Bio::KBase::utilities;
use Bio::KBase::kbaseenv;
use Config::IniFiles;
use warnings;
use JSON::XS;
use Digest::MD5;
use Data::Dumper;
use Getopt::Long;
use Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl;
use Bio::KBase::GenomeAnnotation::Service;


#Initialization function for call
sub util_initialize_call {
	my ($self,$params,$ctx) = @_;
	Bio::KBase::kbaseenv::initialize_call($ctx);
	$Bio::KBase::GenomeAnnotation::Service::CallContext = $ctx;
	Bio::KBase::kbaseenv::ac_client({refresh => 1});
	Bio::KBase::kbaseenv::ga_client({refresh => 1});
	Bio::KBase::kbaseenv::gfu_client({refresh => 1});
	Bio::KBase::kbaseenv::su_client({refresh => 1});
	return $params;
}

sub util_version {
	my ($self) = @_;
	return $VERSION;
}

sub util_log {
	my($self,$message) = @_;
	print $message."\n";
}

sub util_get_genome {
	my ($self,$workspace,$genomeid) = @_;
	my $ref = Bio::KBase::kbaseenv::buildref($workspace,$genomeid);
	my $output = Bio::KBase::kbaseenv::ga_client()->get_genome_v1({
		genomes => [{
			"ref" => $ref
		}],
		ignore_errors => 1,
		no_data => 0,
		no_metadata => 0,
		downgrade => 0
	});
	my $genome = $output->{genomes}->[0];
	$genome->{data}->{'_reference'} = $genome->{info}->[6]."/".$genome->{info}->[0]."/".$genome->{info}->[4];
	print("Genome $genome->{data}->{'_reference'} downloaded\n");
	return $genome->{data};
}

sub max_contigs {
	my $max_contigs = 10000;
	print ("Setting maximum contigs to $max_contigs\n");
	return $max_contigs;
}

sub util_get_contigs {
	my ($self,$workspace,$objid) = @_;
	my $ref = Bio::KBase::kbaseenv::buildref($workspace,$objid);
	my $info = Bio::KBase::kbaseenv::get_object_info([{ref=>$ref}],0);
	my $obj;
	if ($info->[0]->[2] =~ /Assembly/) {
		my $output = Bio::KBase::kbaseenv::ac_client()->get_assembly_as_fasta({
			"ref" => $ref
		});
		my $fasta = "";
		open(my $fh, "<", $output->{path}) || die "Could not find file:".$output->{path};
		while (my $line = <$fh>) {
			$fasta .= $line;
		}
		close($fh);
		$obj = {
			_reference => $info->[0]->[6]."/".$info->[0]->[0]."/".$info->[0]->[4],
			id => $objid,
			name => $objid,
			source_id => $objid,
			source => "KBase",
			type => "SingleGenome",
			contigs => []
		};
		$fasta =~ s/\>([^\n]+)\n/>$1\|\|\|/g;
		$fasta =~ s/\n//g;
		my $array = [split(/\>/,$fasta)];
		my $max_contigs = max_contigs();
		for (my $i=0; $i < @{$array}; $i++) {
			if (@{$obj->{contigs}} > $max_contigs){
				Bio::KBase::Exceptions::ArgumentValidationError->throw(error => 'too many contigs', 
										method_name => 'util_get_contigs');
			}
			if (length($array->[$i]) > 0) {
				my $subarray = [split(/\|\|\|/,$array->[$i])];
				if (@{$subarray} == 2) {
				    my $description = "unknown";
				    my $id = $subarray->[0];
				    if( $subarray->[0] =~ /^([^\s]+)\s(.+)$/) {
				    	$id = $1;
				    	$description = $2;
				    }
				    my $contigobject = {
						id => $id,
						name => $id,
						"length" => length($subarray->[1]),
						md5 => Digest::MD5::md5_hex($subarray->[1]),
						sequence => $subarray->[1],
						description => $description
					};
					$contigobject->{name} = $id;
					$contigobject->{description} = $description;
					push(@{$obj->{contigs}},$contigobject);
	 			}
			}
		}
		my $sortedarray = [sort { $a->{sequence} cmp $b->{sequence} } @{$obj->{contigs}}];
		my $str = "";
		for (my $i=0; $i < @{$sortedarray}; $i++) {
			if (length($str) > 0) {
				$str .= ";";
			}
			$str .= $sortedarray->[$i]->{sequence};
		}
		print("Assembly $obj->{_reference} Downloaded\n");
		$obj->{md5} = Digest::MD5::md5_hex($str);
		$obj->{_kbasetype} = "Assembly";
	} else {
		$obj = Bio::KBase::kbaseenv::get_objects([{
			ref=>Bio::KBase::kbaseenv::buildref($workspace,$objid)
		}]);
		$obj = $obj->[0]->{data};
		$obj->{_kbasetype} = "ContigSet";
		$obj->{_reference} = $info->[0]->[6]."/".$info->[0]->[0]."/".$info->[0]->[4];
		print("Contigset $obj->{_reference}  Downloaded\n");
	}
	my $totallength = 0;
	my $gclength = 0;
	for (my $i=0; $i < @{$obj->{contigs}}; $i++) {
		my $newseq = $obj->{contigs}->[$i]->{sequence};
		$totallength += length($newseq);
		$newseq =~ s/[atAT]//g;
		$gclength += length($newseq);	
	}
	$obj->{_gc} = int(1000*$gclength/$totallength+0.5);
	$obj->{_gc} = $obj->{_gc}/1000;
	return $obj;
}

sub annotate_process {
	my ($self,$parameters) = @_;	
	my $oldfunchash = {};
	my $oldtype     = {};
	#Creating default genome object
	my $inputgenome = {
  		id => $parameters->{output_genome},
  		genetic_code => $parameters->{genetic_code},
  		scientific_name => $parameters->{scientific_name},
  		domain => $parameters->{domain},
  		contigs => [],
  		features => []
  	};
  	my $contigobj;
  	my $message = "";
	my %types = ();

	if (defined($parameters->{input_genome})) {
		$inputgenome = $self->util_get_genome($parameters->{workspace},$parameters->{input_genome});
		for (my $i=0; $i < @{$inputgenome->{features}}; $i++) {
			my $ftr = $inputgenome->{features}->[$i];
			if (!defined($ftr->{type}) || $ftr->{type} lt '     ') {
				if (defined($ftr->{protein_translation})) {
					$ftr->{type} = 'gene';
				} else {
					$ftr->{type} = 'other'; 
				}
			}

			# Reset functions in protein features to "hypothetical protein" to make them available
			# for re-annotation in RAST service (otherwise these features will be skipped).
			if (lc($ftr->{type}) eq "cds" || lc($ftr->{type}) eq "peg" ||
					($ftr->{type} eq "gene" and defined($ftr->{protein_translation}))) {
				if (defined($ftr->{functions})){
					$ftr->{function} = join("; ", @{$ftr->{functions}});
				}
				$oldfunchash->{$ftr->{id}} = $ftr->{function};
				$ftr->{function} = "hypothetical protein";
			}
			elsif ($ftr->{type} eq "gene") {
				$ftr->{type} = 'Non-coding '.$ftr->{type};
			}
			$oldtype->{$ftr->{id}} = $ftr->{type};
			#
			#	Count the input feature types
			#
			if (exists $types{$ftr->{type}}) {
				$types{$ftr->{type}} += 1;
			} else {
				$types{$ftr->{type}} = 1;
			}
		}
		if (exists $inputgenome->{non_coding_features}) {
			for (my $i=0; $i < @{$inputgenome->{non_coding_features}}; $i++) {
				my $ftr = $inputgenome->{non_coding_features}->[$i];
				if (!defined($ftr->{type})) {
					$ftr->{type} = "Non-coding";
				} 
				$oldtype->{$ftr->{id}} = $ftr->{type};
				#
				#	Count the input feature types
				#
				if (exists $types{"Non-coding ".$ftr->{type}}) {
					$types{"Non-coding ".$ftr->{type}} += 1;
				} else {
					$types{"Non-coding ".$ftr->{type}} = 1;
				}
			}
		} else {
			$inputgenome->{non_coding_features} = [];
		}
		
		my $contigref;
		if($inputgenome->{domain} !~ /Eukaryota|Plant/){
		    if (defined($inputgenome->{contigset_ref})) {
			$contigref = $inputgenome->{contigset_ref};
		    } elsif (defined($inputgenome->{assembly_ref})) {
			$contigref = $inputgenome->{assembly_ref};
		    }
			$contigobj = $self->util_get_contigs(undef,$inputgenome->{_reference}.";".$contigref)
		}
		$parameters->{genetic_code} = $inputgenome->{genetic_code};
		$parameters->{domain} = $inputgenome->{domain};
		$parameters->{scientific_name} = $inputgenome->{scientific_name};
	} elsif (defined($parameters->{input_contigset})) {
		$parameters->{genetic_code} = $inputgenome->{genetic_code};
		$parameters->{domain} = $inputgenome->{domain};
		$parameters->{scientific_name} = $inputgenome->{scientific_name};
		$contigobj = $self->util_get_contigs($parameters->{workspace},$parameters->{input_contigset});	
	} else {
		Bio::KBase::utilities::error("Neither contigs nor genome specified!");
	}

	my $tax_domain =  (exists $inputgenome->{domain} && $inputgenome->{domain} =~ m/^([ABV])/o) ? $inputgenome->{domain} : 'U';
	if ($tax_domain eq 'U' ) {
		$message .= "Some RAST tools will not run unless the taxonomic domain is Archaea, Bacteria, or Virus. \nThese tools include: call selenoproteins, call pyrroysoproteins, call crisprs, and call prophage phispy features.\nYou may not get the results you were expecting with your current domain of $inputgenome->{domain}.\n";
	}

	if (defined($contigobj)) {
		my $count = 0;
		my $size = 0;
		if (defined($contigobj->{contigs})) {
			$inputgenome->{contigs} = $contigobj->{contigs};
			for (my $i=0; $i < @{$inputgenome->{contigs}}; $i++) {
				$count++;
				$size += length($inputgenome->{contigs}->[$i]->{sequence});
				$inputgenome->{contigs}->[$i]->{dna} = $inputgenome->{contigs}->[$i]->{sequence};
				delete $inputgenome->{contigs}->[$i]->{sequence};
			}
		}
	
		if ($contigobj->{_kbasetype} eq "ContigSet") {
			$inputgenome->{contigset_ref} = $contigobj->{_reference};
		} else {
			$inputgenome->{assembly_ref} = $contigobj->{_reference};
		}
		if (defined($parameters->{input_contigset})) {
			$message .= "The RAST algorithm was applied to annotating a genome sequence comprised of ".$count." contigs containing ".$size." nucleotides. \nNo initial gene calls were provided.\n";
		} else {
			$message .= "The RAST algorithm was applied to annotating an existing genome: ".$parameters->{scientific_name}.". \nThe sequence for this genome is comprised of ".$count." contigs containing ".$size." nucleotides. \nThe input genome has ".@{$inputgenome->{features}}." existing coding features and ".@{$inputgenome->{non_coding_features}}." existing non-coding features.\n";
			$message .= "NOTE: Older input genomes did not properly separate coding and non-coding features.\n" if (@{$inputgenome->{non_coding_features}} == 0);
		}		
	} else {
		if($inputgenome->{domain} !~ /Eukaryota|Plant/){
		    $message .= "The RAST algorithm was applied to annotating an existing genome: ".$parameters->{scientific_name}.". \nNo DNA sequence was provided for this genome, therefore new genes cannot be called. \nWe can only functionally annotate the ".@{$inputgenome->{features}}." existing features.\n";
		} else {
		    $message .= "The RAST algorithm was applied to functionally annotate ".@{$inputgenome->{features}}." coding features  and ".@{$inputgenome->{non_coding_features}}." existing non-coding features in an existing genome: ".$parameters->{scientific_name}.".\n";
			$message .= "NOTE: Older input genomes did not properly separate coding and non-coding features.\n" if (@{$inputgenome->{non_coding_features}} == 0);
		}
	}
	if  (%types) {
		$message .= "Input genome has the following feature types:\n";
		for my $key (sort keys(%types)) {
			$message .= sprintf("\t%-30s %5d \n", $key, $types{$key});
		}
	}
	
  	my $gaserv = Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl->new();
  	my $workflow = {stages => []};
	my $extragenecalls = "";
	if (defined($parameters->{call_features_rRNA_SEED}) && $parameters->{call_features_rRNA_SEED} == 1)	{
		if (length($extragenecalls) == 0) {
			$extragenecalls = "A scan was conducted for the following additional feature types: ";
		} else {
			$extragenecalls .= "; ";
		}
		$extragenecalls .= "rRNA";
		push(@{$workflow->{stages}},{name => "call_features_rRNA_SEED"});
		if (!defined($contigobj)) {
			Bio::KBase::utilities::error("Cannot call genes on genome with no contigs!");
		}
	}
	if (defined($parameters->{call_features_tRNA_trnascan}) && $parameters->{call_features_tRNA_trnascan} == 1)	{
		if (length($extragenecalls) == 0) {
			$extragenecalls = "A scan was conducted for the following additional feature types: ";
		} else {
			$extragenecalls .= "; ";
		}
		$extragenecalls .= "tRNA";
		push(@{$workflow->{stages}},{name => "call_features_tRNA_trnascan"});
		if (!defined($contigobj)) {
			Bio::KBase::utilities::error("Cannot call genes on genome with no contigs!");
		}
	}
	if (defined($parameters->{call_selenoproteins}) && $parameters->{call_selenoproteins} == 1)	{
		if ($tax_domain ne 'U' ) {
			if (length($extragenecalls) == 0) {
				$extragenecalls = "A scan was conducted for the following additional feature types: ";
			} else {
				$extragenecalls .= "; ";
			}
			$extragenecalls .= "selenoproteins";
			push(@{$workflow->{stages}},{name => "call_selenoproteins"});
			if (!defined($contigobj)) {
				Bio::KBase::utilities::error("Cannot call genes on genome with no contigs!");
			}
		} else {
			$message .= "Did not call selenoproteins because the domain is $parameters->{domain}\n\n";
		}
	}
	if (defined($parameters->{call_pyrrolysoproteins}) && $parameters->{call_pyrrolysoproteins} == 1)	{
		if ($tax_domain ne 'U' ) {
			if (length($extragenecalls) == 0) {
				$extragenecalls = "A scan was conducted for the following additional feature types: ";
			} else {
			$extragenecalls .= "; ";
			}
			$extragenecalls .= "pyrrolysoproteins";
			push(@{$workflow->{stages}},{name => "call_pyrrolysoproteins"});
			if (!defined($contigobj)) {
				Bio::KBase::utilities::error("Cannot call genes on genome with no contigs!");
			} 
		} else {
			$message .= "Did not call pyrrolysoproteins because the domain is $parameters->{domain}\n\n";	
		}
	}
	if (defined($parameters->{call_features_repeat_region_SEED}) && $parameters->{call_features_repeat_region_SEED} == 1)	{
		if (length($extragenecalls) == 0) {
			$extragenecalls = "A scan was conducted for the following additional feature types: ";
		} else {
			$extragenecalls .= "; ";
		}
		$extragenecalls .= "repeat regions";
		push(@{$workflow->{stages}},{
			name => "call_features_repeat_region_SEED",
			"repeat_region_SEED_parameters" => {
							"min_identity" => "95",
							"min_length" => "100"
					 }
		});
		if (!defined($contigobj)) {
			Bio::KBase::utilities::error("Cannot call genes on genome with no contigs!");
		}
	}
	if (defined($parameters->{call_features_strep_suis_repeat}) && $parameters->{call_features_strep_suis_repeat} == 1 && $parameters->{scientific_name} =~ /^Streptococcus\s/)	{
		if (length($extragenecalls) == 0) {
			$extragenecalls = "A scan was conducted for the following additional feature types: ";
		} else {
			$extragenecalls .= "; ";
		}
		$extragenecalls .= "strep suis repeats";
		push(@{$workflow->{stages}},{name => "call_features_strep_suis_repeat"});
		if (!defined($contigobj)) {
			Bio::KBase::utilities::error("Cannot call genes on genome with no contigs!");
		}
	}
	if (defined($parameters->{call_features_strep_pneumo_repeat}) && $parameters->{call_features_strep_pneumo_repeat} == 1 && $parameters->{scientific_name} =~ /^Streptococcus\s/)	{
		if (length($extragenecalls) == 0) {
			$extragenecalls = "A scan was conducted for the following additional feature types: ";
		} else {
			$extragenecalls .= "; ";
		}
		$extragenecalls .= "strep pneumonia repeats";
		push(@{$workflow->{stages}},{name => "call_features_strep_pneumo_repeat"});
		if (!defined($contigobj)) {
			Bio::KBase::utilities::error("Cannot call genes on genome with no contigs!");
		}
	}
	if (defined($parameters->{call_features_crispr}) && $parameters->{call_features_crispr} == 1)	{
		if ($tax_domain ne 'U' ) {
			if (length($extragenecalls) == 0) {
				$extragenecalls = "A scan was conducted for the following additional feature types: ";
			} else {
				$extragenecalls .= "; ";
			}
			$extragenecalls .= "crispr";
			push(@{$workflow->{stages}},{name => "call_features_crispr"});
			if (!defined($contigobj)) {
				Bio::KBase::utilities::error("Cannot call genes on genome with no contigs!");
			}
		} else {
			$message .= "Did not call crisprs because the domain is $parameters->{domain}\n\n";	
		}
	}
	$extragenecalls .= ".\n" if (length($extragenecalls) > 0);

	my $genecalls = "";
	if (defined($parameters->{call_features_CDS_glimmer3}) && $parameters->{call_features_CDS_glimmer3} == 1)	{
		if (@{$inputgenome->{features}} > 0) {
#			$inputgenome->{features} = [];
			$message .= "The existing gene features were cleared due to selection of gene calling with Glimmer3 or Prodigal.\n";
		}
		if (length($genecalls) == 0) {
			$genecalls = "Standard features were called using: ";
		} else {
			$genecalls .= "; ";
		}
		$genecalls .= "glimmer3";
		push(@{$workflow->{stages}},{
			name => "call_features_CDS_glimmer3",
			"glimmer3_parameters" => {
							"min_training_len" => "2000"
					 }
		});
		if (!defined($contigobj)) {
			Bio::KBase::utilities::error("Cannot train and call glimmer genes on a genome with no contigs > 2000 nt!\n");
		}
	}
	if (defined($parameters->{call_features_CDS_prodigal}) && $parameters->{call_features_CDS_prodigal} == 1)	{
		if ($tax_domain ne 'U' ) {
			if (@{$inputgenome->{features}} > 0) {
#				$inputgenome->{features} = [];
				$message .= "The existing gene features were cleared due to selection of gene calling with Glimmer3 or Prodigal.\n";
			}
			if (length($genecalls) == 0) {
				$genecalls = "Standard gene features were called using: ";
			} else {
				$genecalls .= "; ";
			}
			$genecalls .= "prodigal";
			push(@{$workflow->{stages}},{name => "call_features_CDS_prodigal"});
			if (!defined($contigobj)) {
				Bio::KBase::utilities::error("Cannot call genes on genome with no contigs!\n");
			}
		} else {
			$message .= "Did not predict prodigal genes because the domain is $parameters->{domain}\n\n";	
		}
	}
	$genecalls .= ".\n" if (length($genecalls) > 0);

#	if (defined($parameters->{call_features_CDS_genemark}) && $parameters->{call_features_CDS_genemark} == 1)	{
#		if (@{$inputgenome->{features}} > 0) {
##			$inputgenome->{features} = [];
#			$message .= " Existing gene features were cleared due to selection of gene calling with Glimmer3, Prodigal, or Genmark.";
#		}
#		if (length($genecalls) == 0) {
#			$genecalls = "Standard gene features were called using: ";
#		} else {
#			$genecalls .= "; ";
#		}
#		$genecalls .= "genemark";
#		push(@{$workflow->{stages}},{name => "call_features_CDS_genemark"});
#		if (!defined($contigobj)) {
#			Bio::KBase::utilities::error("Cannot call genes on genome with no contigs!");
#		}
#	}
	my $v1flag = 0;
	my $simflag = 0;
	my $annomessage = "";
	if (defined($parameters->{annotate_proteins_kmer_v2}) && $parameters->{annotate_proteins_kmer_v2} == 1)	{
		if (length($annomessage) == 0) {
			$annomessage = "The genome features were functionally annotated using the following algorithm(s): ";
		}
		$annomessage .= "Kmers V2";
		$v1flag = 1;
		$simflag = 1;
		push(@{$workflow->{stages}},{
			name => "annotate_proteins_kmer_v2",
			"kmer_v2_parameters" => {
							"min_hits" => "5",
							"annotate_hypothetical_only" => 1
					 }
		});
	}
	if (defined($parameters->{kmer_v1_parameters}) && $parameters->{kmer_v1_parameters} == 1)	{
		$simflag = 1;
		if (length($annomessage) == 0) {
			$annomessage = "The genome features were functionally annotated using the following algorithm(s): ";
		} else {
			$annomessage .= "; ";
		}
		$annomessage .= "Kmers V1";
		push(@{$workflow->{stages}},{
			name => "annotate_proteins_kmer_v1",
			 "kmer_v1_parameters" => {
							"dataset_name" => "Release70",
							"annotate_hypothetical_only" => $v1flag
					 }
		});
	}
	if (defined($parameters->{annotate_proteins_similarity}) && $parameters->{annotate_proteins_similarity} == 1)	{
		if (length($annomessage) == 0) {
			$annomessage = "The genome features were functionally annotated using the following algorithm(s): ";
		} else {
			$annomessage .= "; ";
		}
		$annomessage .= "protein similarity";
		push(@{$workflow->{stages}},{
			name => "annotate_proteins_similarity",
			"similarity_parameters" => {
							"annotate_hypothetical_only" => $simflag
					 }
		});
	}
	if (defined($parameters->{resolve_overlapping_features}) && $parameters->{resolve_overlapping_features} == 1)	{
		push(@{$workflow->{stages}},{
			name => "resolve_overlapping_features",
			"resolve_overlapping_features_parameters" => {}
		});
	}
	if (defined($parameters->{call_features_prophage_phispy}) && $parameters->{call_features_prophage_phispy} == 1)	{
		if ($tax_domain ne 'U' ) {
			push(@{$workflow->{stages}},{name => "call_features_prophage_phispy"});
		} else {
			$message .= "Did not call call features prophage phispy because the domain is $parameters->{domain}\n\n";	
		}
	}
	$annomessage .= ".\n" if (length($annomessage) > 0);

	if (length($genecalls) > 0) {
		push(@{$workflow->{stages}},{name => "renumber_features"});
		if (@{$inputgenome->{features}} > 0) {
			my $replace = [];
			for (my $i=0; $i< scalar @{$inputgenome->{features}}; $i++) {
				my $ftr = $inputgenome->{features}->[$i];
				if (!defined($ftr->{protein_translation}) || $ftr->{type} =~ /pseudo/) {
					push(@{$replace}, @{$inputgenome->{features}}->[$i]);
				} 
			}
			$inputgenome->{features} = $replace;
		}
		$message .= $genecalls;
	}
	if (length($extragenecalls) > 0) {
		$message .= $extragenecalls;
	}
	if (length($annomessage) > 0) {
		$message .= $annomessage;
	}

	my $genome = $inputgenome;
	my $genehash = {};
	my $num_coding = 0;
	if (defined($genome->{features})) {
		for (my $i=0; $i < @{$genome->{features}}; $i++) {
			# Caching feature functions for future comparison against new functions
			# defined by RAST service in order to calculate number of updated features.
			# If function is not set we treat it as empty string to avoid perl warning.
			my $func = $genome->{features}->[$i]->{function};
			if (not defined($func)) {
				$func = "";
			}
			$genehash->{$genome->{features}->[$i]->{id}}->{$func} = 1;
		}
	}
	if (defined($genome->{non_coding_features})) {
		for (my $i=0; $i < @{$genome->{non_coding_features}}; $i++) {
			$genehash->{$genome->{non_coding_features}->[$i]->{id}} = 1;
		}
	}
	if (defined($inputgenome->{features})) {
		## Checking if it's recent old genome containing both CDSs and genes in "features" array
		if ((not defined($inputgenome->{cdss})) && (not defined($inputgenome->{mrnas}))) {
			my $genes = [];
			my $non_genes = [];
			for (my $i=0; $i < @{$inputgenome->{features}}; $i++) {
				my $feature = $inputgenome->{features}->[$i];
				if ($feature->{type} eq "gene" and not(defined($feature->{protein_translation}))) {
					# gene without protein translation
					push(@{$genes}, $feature);
				} else {
					push(@{$non_genes}, $feature);
				}
			}
			if ((scalar @{$genes} > 0) && (scalar @{$non_genes} > 0)) {
				# removing genes from features
				$inputgenome->{features} = $non_genes;
			}
		}
		## Temporary switching feature type from 'gene' to 'CDS':
		for (my $i=0; $i < @{$inputgenome->{features}}; $i++) {
			my $feature = $inputgenome->{features}->[$i];
			if ($feature->{type} eq "gene" and defined($feature->{protein_translation})) {
				$feature->{type} = "CDS";
			}
		}
		# When RAST annotates the genome it becomes incompatible with the new
		# spec. Removing this attribute triggers an "upgrade" to the genome
		# that fixes these problems when saving with GFU
		delete $inputgenome->{feature_counts};
		if (defined($inputgenome->{ontology_events})){
			my $ont_event = {
				 "id" => "SSO",
				  "method" => Bio::KBase::utilities::method(),
				  "method_version" => $self->util_version(),
				  "ontology_ref" => "KBaseOntology/seed_subsystem_ontology",
				  "timestamp" => Bio::KBase::utilities::timestamp()
			};
			push(@{$inputgenome->{ontology_events}}, $ont_event);
		}
	}
	#----------------------
	# Runs, the annotation, comment out if you dont have the reference files
	#---------------------
	$genome = $gaserv->run_pipeline($inputgenome, $workflow);

	delete $genome->{contigs};
	delete $genome->{feature_creation_event};
	delete $genome->{analysis_events};
	$genome->{genetic_code} = $genome->{genetic_code}+0;
	$genome->{id} = $parameters->{output_genome};
	if (!defined($genome->{source})) {
		$genome->{source} = "KBase";
		$genome->{source_id} = $parameters->{output_genome};
	}
	if (defined($inputgenome->{gc_content})) {
		$genome->{gc_content} = $inputgenome->{gc_content};
	}
	if (defined($genome->{gc})) {
		$genome->{gc_content} = $genome->{gc}+0;
		delete $genome->{gc};
	}
	if (!defined($genome->{gc_content})) {
		$genome->{gc_content} = 0.5;
	}
	if (not defined($genome->{non_coding_features})) {
		$genome->{non_coding_features} = [];
	}
	my @splice_list = ();
	if (defined($genome->{features})) {
		for (my $i=0; $i < scalar @{$genome->{features}}; $i++) {
			my $ftr = $genome->{features}->[$i];
			if (defined($ftr->{aliases}) && scalar @{$ftr->{aliases}} > 0)  {
				if (ref($ftr->{aliases}->[0]) !~ /ARRAY/) {
				# Found some pseudogenes that have wrong structure for aliases
					my $tmp = [];
					foreach my $key (@{$ftr->{aliases}})  {
						my @ary = ('alias',$key);
						push(@{$tmp},\@ary); 
					}
					$ftr->{aliases} = $tmp;
				}
			 }
			if (defined ($ftr->{type}) && $ftr->{type} ne 'gene' && $ftr->{type} ne 'CDS') {
				push(@splice_list,$i);
            }
		}
	}
	my $count = 0;

#	
#	Move non-coding features from features to non_coding_features
#	They can have more than one location and they need md5 and dna_sequence_length
#
#	print "Array size =".scalar @{$genome->{features}}."\n";
#	print "Number to splice out = ".scalar @splice_list."\n";
	foreach my $key (reverse @splice_list) {
		if ($key =~ /\D/ ) {
			print "INVALID:size=".scalar @{$genome->{features}}."\n";
#			print Dumper $key;
		}
		else {
			my $ftr = $genome->{features}->[$key];
			if (defined($ftr->{location})) {
				$ftr->{dna_sequence_length} = 0;
				$ftr->{md5} = '';
				for (my $i=0; $i < scalar @{$ftr->{location}}; $i++) {
					$ftr->{location}->[$i]->[1]  = $ftr->{location}->[$i]->[1]+0;
					$ftr->{location}->[$i]->[3]  = $ftr->{location}->[$i]->[3]+0;
					$ftr->{dna_sequence_length} += $ftr->{location}->[$i]->[3]+0;
				}
			}
			delete  $ftr->{feature_creation_event};
			my $non = splice(@{$genome->{features}},$key,1);
			push(@{$genome->{non_coding_features}}, $non);

			$count++;
		}

	}
	
	if (defined($contigobj) && defined($contigobj->{contigs}) && scalar(@{$contigobj->{contigs}})>0 ) {
		$genome->{num_contigs} = @{$contigobj->{contigs}};
		$genome->{md5} = $contigobj->{md5};
	}
	#Getting the seed ontology dictionary
	#Sam: I need to build the PlantSEED ontology and use it here
	my $output = Bio::KBase::kbaseenv::get_objects([{
		workspace => "KBaseOntology",
		name => "seed_subsystem_ontology"
	}]);
	#Building a hash of standardized seed function strings
	my $funchash = {};
	foreach my $term (keys(%{$output->[0]->{data}->{term_hash}})) {
		my $rolename = lc($output->[0]->{data}->{term_hash}->{$term}->{name});
		$rolename =~ s/[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+//g;
		$rolename =~ s/\s//g;
		$rolename =~ s/\#.*$//g;
		$funchash->{$rolename} = $output->[0]->{data}->{term_hash}->{$term};
	}
	my $newftrs = 0;
	my $newncfs = 0;
	my $update_cdss = 'N';
	my $proteins = 0;
	my $others = 0;
	my $seedfunctions = 0;
	my $genomefunchash;
	my $seedfunchash;
	my $advancedmessage = '';
	%types = ();
	if (defined($genome->{features})) {
		for (my $i=0; $i < @{$genome->{features}}; $i++) {
			my $ftr = $genome->{features}->[$i];
			if (defined($genehash) && !defined($genehash->{$ftr->{id}})) {
				# Let's count number of features with functions updated by RAST service.
				# If function is not set we treat it as empty string to avoid perl warning.
				$newftrs++;
				my $func = $ftr->{function};
				if (not defined($func)) {
					$func = "";
				}
				if ($ftr->{type} eq 'CDS' && !defined($genome->{cdss}->[$i])) {
					#Some of the optional gene finders don't respect new genome format
					#They may not have the cdss
					$update_cdss = 'Y';
				}
			} else {
				$num_coding++;
			}
			if (!defined($ftr->{type}) && $ftr->{id} =~ m/(\w+)\.\d+$/) {
				$ftr->{type} = $1;
			}
			if (defined($ftr->{protein_translation})) {
				$ftr->{protein_translation_length} = length($ftr->{protein_translation})+0;
				$ftr->{md5} = Digest::MD5::md5_hex($ftr->{protein_translation});
				$proteins++;
			} else {
				$others++;
			}
			if (defined($ftr->{dna_sequence})) {
				$ftr->{dna_sequence_length} = length($ftr->{dna_sequence})+0;
			}
			if (defined($ftr->{quality}->{weighted_hit_count})) {
				$ftr->{quality}->{weighted_hit_count} = $ftr->{quality}->{weighted_hit_count}+0;
			}
			if (defined($ftr->{quality}->{hit_count})) {
				$ftr->{quality}->{hit_count} = $ftr->{quality}->{hit_count}+0;
			}
			if (defined($ftr->{annotations})) {
				delete $ftr->{annotations};
			}
			if (defined($ftr->{location})) {
				$ftr->{location}->[0]->[1] = $ftr->{location}->[0]->[1]+0;
				$ftr->{location}->[0]->[3] = $ftr->{location}->[0]->[3]+0;
			}
			delete $ftr->{feature_creation_event};
			if (defined($ftr->{function}) && length($ftr->{function}) > 0 && defined($ftr->{protein_translation})) {
				for (my $j=0; $j < @{$genome->{features}}; $j++) {
					if ($i ne $j) {
						if ($ftr->{location}->[0]->[0] eq $genome->{features}->[$j]->{location}->[0]->[0]) {
							if ($ftr->{location}->[0]->[1] == $genome->{features}->[$j]->{location}->[0]->[1]) {
								if ($ftr->{location}->[0]->[2] eq $genome->{features}->[$j]->{location}->[0]->[2]) {
									if ($ftr->{location}->[0]->[3] == $genome->{features}->[$j]->{location}->[0]->[3]) {
										if ($ftr->{type} ne $genome->{features}->[$j]->{type}) {
											$genome->{features}->[$j]->{function} = $ftr->{function};
										}
									}
								}
							}
						}
					}
				}
			}
		}
		for (my $i=0; $i < @{$genome->{features}}; $i++) {
			my $ftr = $genome->{features}->[$i];
			if (defined($oldfunchash->{$ftr->{id}}) && (!defined($ftr->{function}) || $ftr->{function} =~ /hypothetical\sprotein/)) {
				if (defined($parameters->{retain_old_anno_for_hypotheticals}) && $parameters->{retain_old_anno_for_hypotheticals} == 1)	{
					$ftr->{function} = $oldfunchash->{$ftr->{id}};
				}
			}
			if (defined($ftr->{function}) && length($ftr->{function}) > 0) {
				my $function = $ftr->{function};
				my $array = [split(/\#/,$function)];
				$function = shift(@{$array});
				$function =~ s/\s+$//;
				$array = [split(/\s*;\s+|\s+[\@\/]\s+/,$function)];
				my $marked = 0;
				for (my $j=0;$j < @{$array}; $j++) {
					my $rolename = lc($array->[$j]);
					$rolename =~ s/[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]+//g;
					$rolename =~ s/\s//g;
					$rolename =~ s/\#.*$//g;
					$genomefunchash->{$rolename} = 1;
					if (defined($funchash->{$rolename})) {
						$seedfunchash->{$rolename} = 1;
						if ($marked == 0) {
							$seedfunctions++;
							$marked = 1;
						}
						my $ont_term = $ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}};
						if (!defined($ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}})) {
							$ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}} = {
								 evidence => [],
								 id => $funchash->{$rolename}->{id},
								 term_name => $funchash->{$rolename}->{name},
								 ontology_ref => "KBaseOntology/seed_subsystem_ontology",
								 term_lineage => [],
							};
						}
						my $found = 0;
						if (ref($ont_term) eq 'ARRAY'){
							push(@{$ont_term}, $#{$inputgenome->{ontology_events}});
						} else {
							for (my $k = 0; $k < @{$ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}}; $k++) {
								if ($ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}->[$k]->{method} eq Bio::KBase::utilities::method()) {
									$ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}->[$k]->{timestamp} = Bio::KBase::utilities::timestamp();
									$ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}->[$k]->{method_version} = $self->util_version();
									$found = 1;
									last;
								}
							}
							if ($found == 0) {
								push(
									@{$ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}},
									{
										method         => Bio::KBase::utilities::method(),
										method_version => $self->util_version(),
										timestamp      => Bio::KBase::utilities::timestamp()
									});
							}
						}
						if (exists ($ftr->{ontology_terms}->{SSO})) {
							foreach my $sso (keys($ftr->{ontology_terms}->{SSO})) {
								if ($sso =~ /SSO:000009304/) {
									$advancedmessage .= "Found selenocysteine-containing gene $ftr->{ontology_terms}->{SSO}->{$sso}\n";
								}
								elsif ($sso =~ /SSO:000009291/) {
									$advancedmessage .= "Found pyrrolysine-containing gene $ftr->{ontology_terms}->{SSO}->{$sso}\n";
								}
							}
						}	
					}
				}
			}
		}

		## Rolling protein features back from 'CDS' to 'gene':
		for (my $i=0; $i < @{$genome->{features}}; $i++) {
			my $ftr = $genome->{features}->[$i];
			if ($ftr->{type} eq "CDS") {
				$ftr->{type} = "gene"
			}
			if (exists $oldtype->{$ftr->{id}} && $oldtype->{$ftr->{id}} =~ /gene/) {
				$ftr->{type} = $oldtype->{$ftr->{id}};					
			}
			if (defined($ftr->{location})) {
				$ftr->{location}->[0]->[1] = $ftr->{location}->[0]->[1]+0;
				$ftr->{location}->[0]->[3] = $ftr->{location}->[0]->[3]+0;
			}
#			$ftr->{type} = 'gene';
			my $type = '';
			if ( defined$ftr->{type}) {
				$type = 'Coding '.$ftr->{type};
			} else {
				$type = 'Coding ';
			}
			#	Count the output feature types
			if (exists $types{$type}) {
				$types{$type} += 1;
			} else {
				$types{$type} = 1;
			}
		}
		if ((not defined($genome->{cdss})) || (not defined($genome->{mrnas})) || $update_cdss eq 'Y') {
			## Reconstructing new feature arrays ('cdss' and 'mrnas') if they are not present:
			my $cdss = [];
			$genome->{cdss} = $cdss;
			my $mrnas = [];
			$genome->{mrnas} = $mrnas;
			for (my $i=0; $i < @{$genome->{features}}; $i++) {
				my $feature = $genome->{features}->[$i];
				if (defined($feature->{protein_translation})) {
					my $gene_id = $feature->{id};
					my $cds_id = $gene_id . "_CDS";
					my $mrna_id = $gene_id . "_mRNA";
					my $location = $feature->{location};
					my $md5 = $feature->{md5};
					my $function = $feature->{function};
					if (not defined($function)) {
						$function = "";
					}
					my $ontology_terms = $feature->{ontology_terms};
					if (not defined($ontology_terms)) {
						$ontology_terms = {};
					}
					my $protein_translation = $feature->{protein_translation};
					my $protein_translation_length = length($protein_translation);
					my $aliases = [];
					if (defined($feature->{aliases})) {
						$aliases = $feature->{aliases};
						# Found some pseudogenes that have wrong structure for aliases
						if (ref($aliases) !~ /ARRAY/) {
							$aliases = [$feature->{aliases}];
						}
					}
					push(@{$cdss}, {
						id => $cds_id,
						location => $location,
						md5 => $md5,
						parent_gene => $gene_id,
						parent_mrna => $mrna_id,
						function => $function,
						ontology_terms => $ontology_terms,
						protein_translation => $protein_translation,
						protein_translation_length => $protein_translation_length,
						aliases => $aliases
					});
					push(@{$mrnas}, {
						id => $mrna_id,
						location => $location,
						md5 => $md5,
						parent_gene => $gene_id,
						cds => $cds_id
					});
					$feature->{cdss} = [$cds_id];
					$feature->{mrnas} = [$mrna_id];
				}
			}
		}
	}
	my $num_non_coding = 0;
	if (defined($genome->{non_coding_features})) {
		for (my $i=0; $i < @{$genome->{non_coding_features}}; $i++) {
			my $ftr = $genome->{non_coding_features}->[$i];
			if (defined($genehash) && !defined($genehash->{$ftr->{id}})) {
				# Let's count number of non_coding_features with functions updated by RAST service.
				# If function is not set we treat it as empty string to avoid perl warning.
				$newncfs++;
				$newftrs++;
				my $func = $ftr->{function};
				if (not defined($func)) {
					$func = "";
				}

			} else {
				$num_non_coding++; 
			}
			my $type = '';
			if ( defined$ftr->{type} && $ftr->{type} =~ /coding/) {
				$type = $ftr->{type} ;
			} elsif (defined$ftr->{type}) {
				$type = 'Non-coding '.$ftr->{type} ;
			} else {
				$type = 'Non-coding ';
			}
			#	Count the output feature types
			if (exists $types{$type}) {
				$types{$type} += 1;
			} else {
				$types{$type} = 1;
			}
		}
	}
	if (defined($inputgenome)) {
		$message .= "In addition to the remaining original $num_coding coding features and $num_non_coding non-coding features, ".$newftrs." new features were called, of which $newncfs are non-coding.\n";
		if  (%types) {
			$message .= "Output genome has the following feature types:\n";
			for my $key (sort keys(%types)) {
				$message .= sprintf("\t%-30s %5d \n", $key, $types{$key});
			}
		}
	}

	$message .= "Overall, the genes have ".keys(%{$genomefunchash})." distinct functions. \nThe genes include ".$seedfunctions." genes with a SEED annotation ontology across ".keys(%{$seedfunchash})." distinct SEED functions.\n";
	$message .= "The number of distinct functions can exceed the number of genes because some genes have multiple functions.\n";
	print($message);
	if (!defined($genome->{assembly_ref})) {
		delete $genome->{assembly_ref};
	}
	if (defined($contigobj)) {
		$genome->{gc_content} = $contigobj->{_gc};
		if ($contigobj->{_kbasetype} eq "ContigSet") {
			$genome->{contigset_ref} = $contigobj->{_reference};
		} else {
			$genome->{assembly_ref} = $contigobj->{_reference};
		}
	}
#	print "SEND OFF FOR SAVING\n";
#	print "***** Number of features=".scalar  @{$genome->{features}}."\n";
#	print "***** Number of non_coding_features=".scalar  @{$genome->{non_coding_features}}."\n";
#	print "***** Number of cdss=    ".scalar  @{$genome->{cdss}}."\n";
#	print "***** Number of mrnas=   ".scalar  @{$genome->{mrnas}}."\n";
	my $gaout = Bio::KBase::kbaseenv::gfu_client()->save_one_genome({
		workspace => $parameters->{workspace},
        name => $parameters->{output_genome},
        data => $genome,
        provenance => [{
			"time" => DateTime->now()->datetime()."+0000",
			service_ver => $self->util_version(),
			service => "RAST_SDK",
			method => Bio::KBase::utilities::method(),
			method_params => [$parameters],
			input_ws_objects => [],
			resolved_ws_objects => [],
			intermediate_incoming => [],
			intermediate_outgoing => []
		}],
        hidden => 0
	});
	Bio::KBase::kbaseenv::add_object_created({
		"ref" => $gaout->{info}->[6]."/".$gaout->{info}->[0]."/".$gaout->{info}->[4],
		"description" => "Annotated genome"
	});
	Bio::KBase::utilities::print_report_message({
		message => "<pre>".$message."</pre>",
		append => 0,
		html => 0
	});
	return ({"ref" => $gaout->{info}->[6]."/".$gaout->{info}->[0]."/".$gaout->{info}->[4]},$message);
}
#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR
	Bio::KBase::utilities::read_config({
		service => "RAST_SDK",
		mandatory => ['workspace-url']
	});
    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}

=head1 METHODS



=head2 annotate_genome

  $return = $obj->annotate_genome($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a RAST_SDK.AnnotateGenomeParams
$return is a RAST_SDK.AnnotateGenomeResults
AnnotateGenomeParams is a reference to a hash where the following keys are defined:
	workspace has a value which is a string
	input_genome has a value which is a RAST_SDK.genome_id
	input_contigset has a value which is a RAST_SDK.contigset_id
	genetic_code has a value which is an int
	domain has a value which is a string
	scientific_name has a value which is a string
	output_genome has a value which is a string
	call_features_rRNA_SEED has a value which is a RAST_SDK.bool
	call_features_tRNA_trnascan has a value which is a RAST_SDK.bool
	call_selenoproteins has a value which is a RAST_SDK.bool
	call_pyrrolysoproteins has a value which is a RAST_SDK.bool
	call_features_repeat_region_SEED has a value which is a RAST_SDK.bool
	call_features_insertion_sequences has a value which is a RAST_SDK.bool
	call_features_strep_suis_repeat has a value which is a RAST_SDK.bool
	call_features_strep_pneumo_repeat has a value which is a RAST_SDK.bool
	call_features_crispr has a value which is a RAST_SDK.bool
	call_features_CDS_glimmer3 has a value which is a RAST_SDK.bool
	call_features_CDS_prodigal has a value which is a RAST_SDK.bool
	call_features_CDS_genemark has a value which is a RAST_SDK.bool
	annotate_proteins_kmer_v2 has a value which is a RAST_SDK.bool
	kmer_v1_parameters has a value which is a RAST_SDK.bool
	annotate_proteins_similarity has a value which is a RAST_SDK.bool
	resolve_overlapping_features has a value which is a RAST_SDK.bool
	call_features_prophage_phispy has a value which is a RAST_SDK.bool
	retain_old_anno_for_hypotheticals has a value which is a RAST_SDK.bool
genome_id is a string
contigset_id is a string
bool is an int
AnnotateGenomeResults is a reference to a hash where the following keys are defined:
	workspace has a value which is a RAST_SDK.workspace_name
	id has a value which is a string
	report_name has a value which is a string
	report_ref has a value which is a string
workspace_name is a string

</pre>

=end html

=begin text

$params is a RAST_SDK.AnnotateGenomeParams
$return is a RAST_SDK.AnnotateGenomeResults
AnnotateGenomeParams is a reference to a hash where the following keys are defined:
	workspace has a value which is a string
	input_genome has a value which is a RAST_SDK.genome_id
	input_contigset has a value which is a RAST_SDK.contigset_id
	genetic_code has a value which is an int
	domain has a value which is a string
	scientific_name has a value which is a string
	output_genome has a value which is a string
	call_features_rRNA_SEED has a value which is a RAST_SDK.bool
	call_features_tRNA_trnascan has a value which is a RAST_SDK.bool
	call_selenoproteins has a value which is a RAST_SDK.bool
	call_pyrrolysoproteins has a value which is a RAST_SDK.bool
	call_features_repeat_region_SEED has a value which is a RAST_SDK.bool
	call_features_insertion_sequences has a value which is a RAST_SDK.bool
	call_features_strep_suis_repeat has a value which is a RAST_SDK.bool
	call_features_strep_pneumo_repeat has a value which is a RAST_SDK.bool
	call_features_crispr has a value which is a RAST_SDK.bool
	call_features_CDS_glimmer3 has a value which is a RAST_SDK.bool
	call_features_CDS_prodigal has a value which is a RAST_SDK.bool
	call_features_CDS_genemark has a value which is a RAST_SDK.bool
	annotate_proteins_kmer_v2 has a value which is a RAST_SDK.bool
	kmer_v1_parameters has a value which is a RAST_SDK.bool
	annotate_proteins_similarity has a value which is a RAST_SDK.bool
	resolve_overlapping_features has a value which is a RAST_SDK.bool
	call_features_prophage_phispy has a value which is a RAST_SDK.bool
	retain_old_anno_for_hypotheticals has a value which is a RAST_SDK.bool
genome_id is a string
contigset_id is a string
bool is an int
AnnotateGenomeResults is a reference to a hash where the following keys are defined:
	workspace has a value which is a RAST_SDK.workspace_name
	id has a value which is a string
	report_name has a value which is a string
	report_ref has a value which is a string
workspace_name is a string


=end text



=item Description

annotate genome
params - a param hash that includes the workspace id and options

=back

=cut

sub annotate_genome
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_genome');
    }

    my $ctx = $RAST_SDK::RAST_SDKServer::CallContext;
    my($return);
    #BEGIN annotate_genome
    $self->util_initialize_call($params,$ctx);
    $params = Bio::KBase::utilities::args($params,["workspace","output_genome"],{
	    input_genome => undef,
	    input_contigset => undef,
	    genetic_code => 11,
	    domain => "Bacteria",
	    scientific_name => "Unknown species",
	    call_features_rRNA_SEED => 1,
	    call_features_tRNA_trnascan => 1,
	    call_selenoproteins => 1,
	    call_pyrrolysoproteins => 1,
	    call_features_repeat_region_SEED => 1,
	    call_features_strep_suis_repeat => 1,
	    call_features_strep_pneumo_repeat => 1,
	    call_features_crispr => 1,
	    call_features_CDS_glimmer3 => 1,
	    call_features_CDS_prodigal => 1,
	    call_features_CDS_genemark => 1,
	    annotate_proteins_kmer_v2 => 1,
	    kmer_v1_parameters => 1,
	    annotate_proteins_similarity => 1,
	    resolve_overlapping_features => 1,
	    call_features_prophage_phispy => 1,
	    retain_old_anno_for_hypotheticals => 1
	});
    my $output = $self->annotate_process($params);
    my $reportout = Bio::KBase::kbaseenv::create_report({
    	workspace_name => $params->{workspace},
    	report_object_name => $params->{output_genome}.".report",
    });
	$return = {
    	workspace => $params->{workspace},
    	id => $params->{output_genome},
    	report_ref => $reportout->{"ref"},
    	report_name => $params->{output_genome}.".report",
    	ws_report_id => $params->{output_genome}.".report"
    };
    Bio::KBase::utilities::close_debug();
    #END annotate_genome
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_genome:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_genome');
    }
    return($return);
}




=head2 annotate_genomes

  $return = $obj->annotate_genomes($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a RAST_SDK.AnnotateGenomesParams
$return is a RAST_SDK.AnnotateGenomesResults
AnnotateGenomesParams is a reference to a hash where the following keys are defined:
	workspace has a value which is a string
	genomes has a value which is a reference to a list where each element is a RAST_SDK.GenomeParams
	call_features_rRNA_SEED has a value which is a RAST_SDK.bool
	call_features_tRNA_trnascan has a value which is a RAST_SDK.bool
	call_selenoproteins has a value which is a RAST_SDK.bool
	call_pyrrolysoproteins has a value which is a RAST_SDK.bool
	call_features_repeat_region_SEED has a value which is a RAST_SDK.bool
	call_features_insertion_sequences has a value which is a RAST_SDK.bool
	call_features_strep_suis_repeat has a value which is a RAST_SDK.bool
	call_features_strep_pneumo_repeat has a value which is a RAST_SDK.bool
	call_features_crispr has a value which is a RAST_SDK.bool
	call_features_CDS_glimmer3 has a value which is a RAST_SDK.bool
	call_features_CDS_prodigal has a value which is a RAST_SDK.bool
	call_features_CDS_genemark has a value which is a RAST_SDK.bool
	annotate_proteins_kmer_v2 has a value which is a RAST_SDK.bool
	kmer_v1_parameters has a value which is a RAST_SDK.bool
	annotate_proteins_similarity has a value which is a RAST_SDK.bool
	resolve_overlapping_features has a value which is a RAST_SDK.bool
	call_features_prophage_phispy has a value which is a RAST_SDK.bool
	retain_old_anno_for_hypotheticals has a value which is a RAST_SDK.bool
GenomeParams is a reference to a hash where the following keys are defined:
	input_contigset has a value which is a RAST_SDK.contigset_id
	input_genome has a value which is a RAST_SDK.genome_id
	output_genome has a value which is a RAST_SDK.genome_id
	genetic_code has a value which is an int
	domain has a value which is a string
	scientific_name has a value which is a string
contigset_id is a string
genome_id is a string
bool is an int
AnnotateGenomesResults is a reference to a hash where the following keys are defined:
	workspace has a value which is a RAST_SDK.workspace_name
	report_name has a value which is a string
	report_ref has a value which is a string
workspace_name is a string

</pre>

=end html

=begin text

$params is a RAST_SDK.AnnotateGenomesParams
$return is a RAST_SDK.AnnotateGenomesResults
AnnotateGenomesParams is a reference to a hash where the following keys are defined:
	workspace has a value which is a string
	genomes has a value which is a reference to a list where each element is a RAST_SDK.GenomeParams
	call_features_rRNA_SEED has a value which is a RAST_SDK.bool
	call_features_tRNA_trnascan has a value which is a RAST_SDK.bool
	call_selenoproteins has a value which is a RAST_SDK.bool
	call_pyrrolysoproteins has a value which is a RAST_SDK.bool
	call_features_repeat_region_SEED has a value which is a RAST_SDK.bool
	call_features_insertion_sequences has a value which is a RAST_SDK.bool
	call_features_strep_suis_repeat has a value which is a RAST_SDK.bool
	call_features_strep_pneumo_repeat has a value which is a RAST_SDK.bool
	call_features_crispr has a value which is a RAST_SDK.bool
	call_features_CDS_glimmer3 has a value which is a RAST_SDK.bool
	call_features_CDS_prodigal has a value which is a RAST_SDK.bool
	call_features_CDS_genemark has a value which is a RAST_SDK.bool
	annotate_proteins_kmer_v2 has a value which is a RAST_SDK.bool
	kmer_v1_parameters has a value which is a RAST_SDK.bool
	annotate_proteins_similarity has a value which is a RAST_SDK.bool
	resolve_overlapping_features has a value which is a RAST_SDK.bool
	call_features_prophage_phispy has a value which is a RAST_SDK.bool
	retain_old_anno_for_hypotheticals has a value which is a RAST_SDK.bool
GenomeParams is a reference to a hash where the following keys are defined:
	input_contigset has a value which is a RAST_SDK.contigset_id
	input_genome has a value which is a RAST_SDK.genome_id
	output_genome has a value which is a RAST_SDK.genome_id
	genetic_code has a value which is an int
	domain has a value which is a string
	scientific_name has a value which is a string
contigset_id is a string
genome_id is a string
bool is an int
AnnotateGenomesResults is a reference to a hash where the following keys are defined:
	workspace has a value which is a RAST_SDK.workspace_name
	report_name has a value which is a string
	report_ref has a value which is a string
workspace_name is a string


=end text



=item Description

annotate genomes
params - a param hash that includes the workspace id and options

=back

=cut

sub annotate_genomes
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_genomes:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_genomes');
    }

    my $ctx = $RAST_SDK::RAST_SDKServer::CallContext;
    my($return);
    #BEGIN annotate_genomes
    $self->util_initialize_call($params,$ctx);
    #$params = Bio::KBase::utilities::args($params,["workspace"],{
    $params = Bio::KBase::utilities::args($params,["workspace","output_genome"],{
    	input_genomes => [],
    	genome_text => undef,
	    call_features_rRNA_SEED => 1,
	    call_features_tRNA_trnascan => 1,
	    call_selenoproteins => 1,
	    call_pyrrolysoproteins => 1,
	    call_features_repeat_region_SEED => 1,
	    call_features_strep_suis_repeat => 1,
	    call_features_strep_pneumo_repeat => 1,
	    call_features_crispr => 1,
	    call_features_CDS_glimmer3 => 1,
	    call_features_CDS_prodigal => 1,
	    call_features_CDS_genemark => 1,
	    annotate_proteins_kmer_v2 => 1,
	    kmer_v1_parameters => 1,
	    annotate_proteins_similarity => 1,
	    resolve_overlapping_features => 1,
	    call_features_prophage_phispy => 1,
	    retain_old_anno_for_hypotheticals => 1
	});
	my $htmlmessage = "";
	my $genomes = $params->{input_genomes};

	my $obj_type;
	if (ref $genomes eq 'ARRAY') {
		my $replace_genomes = [];
		foreach my $ref (@$genomes) {
	 		my $info = Bio::KBase::kbaseenv::get_object_info([{ref=>$ref}],0);
			$obj_type = $info->[0]->[2];
			if ($obj_type =~ /KBaseSearch\.GenomeSet/ || $obj_type =~ /KBaseSets\.AssemblySet/) {
				my $obj = Bio::KBase::kbaseenv::get_objects([{
					ref=>Bio::KBase::kbaseenv::buildref($params->{workspace},$ref)}])->[0]->{'refs'};

				if (ref($obj) eq 'HASH') {
					foreach my $key (keys %$obj) {
						push(@$replace_genomes,$key);
					}
				} elsif (ref($obj) eq 'ARRAY') {
					foreach my $key (@$obj) {
						push(@$replace_genomes,$key);
					}
					
				}
			} else {
				push(@$replace_genomes,$ref);
			}
		}
		$genomes = $replace_genomes;
	}
	
	if (defined($params->{genome_text})) {
		my $new_genome_list = [split(/[\n;\|]+/,$params->{genome_text})];
		for (my $i=0; $i < @{$new_genome_list}; $i++) {
			push(@{$genomes},$new_genome_list->[$i]);
		}
	}

	my $output_genomes = [];
	for (my $i=0; $i < @{$genomes}; $i++) {
		my $obj_type = '';
		my $input = $genomes->[$i];
		if ($input =~ m/\//) {
			my $array = [split(/\//,$input)];
			my $info = Bio::KBase::kbaseenv::get_object_info([{ref=>
				Bio::KBase::kbaseenv::buildref($array->[0],$array->[1],$array->[2])}
			],0);
			$input = $info->[0]->[1];
			$obj_type =  $info->[0]->[2];
		} else {
			my $info = Bio::KBase::kbaseenv::get_object_info([{ref=>$input}]);
			$obj_type =  $info->[0]->[2];
		}

		my $currentparams = Bio::KBase::utilities::args({},[],{
			output_genome => $input.".RAST",
			input_genome => $genomes->[$i],
		    input_contigset => undef,
		    genetic_code => 11,
		    domain => "Bacteria",
		    scientific_name => "unknown taxon"
		});
		my $list = [qw(
			workspace
			scientific_name
			genetic_name
			domain
			call_features_rRNA_SEED
		    call_features_tRNA_trnascan
		    call_selenoproteins
		    call_pyrrolysoproteins
		    call_features_repeat_region_SEED
		    call_features_strep_suis_repeat
		    call_features_strep_pneumo_repeat
		    call_features_crispr
		    call_features_CDS_glimmer3
		    call_features_CDS_prodigal
		    annotate_proteins_kmer_v2
		    kmer_v1_parameters
		    annotate_proteins_similarity
		    resolve_overlapping_features
		    call_features_prophage_phispy
		    retain_old_anno_for_hypotheticals
		)];

		for (my $j=0; $j < @{$list}; $j++) {
			$currentparams->{$list->[$j]} = $params->{$list->[$j]} if (exists $params->{$list->[$j]});
		}

		if ($obj_type =~ /KBaseGenomeAnnotations\.Assembly/) {
			$currentparams->{'input_contigset'} = delete $currentparams->{'input_genome'};
			delete $currentparams->{'retain_old_anno_for_hypotheticals'};
		}
			
		eval {
			my ($output,$message) = $self->annotate_process($currentparams);
			push(@$output_genomes,$output->{ref});
			$htmlmessage .= $message;
		};
		if ($@) {
			$htmlmessage .= $input." failed!<br>\n\n";
		} else {
			$htmlmessage .= $input." succeeded!<br>\n\n";
		}
	}
		
		my $output_genomeset;
		if (defined $params->{output_genome} && $params->{output_genome} gt ' ') { 
			my $output_genomeset = $params->{output_genome};
	        my $genome_set_name = $params->{output_genome};
	        my $genome_set = Bio::KBase::kbaseenv::su_client()->KButil_Build_GenomeSet({
	            workspace_name => $params->{workspace},
	            input_refs => $output_genomes,
	            output_name => $output_genomeset,
	            desc => 'GenomeSet Description'
	        });
		}


	my $path = "/kb/module/work/tmp/microbial_genome_report.$params->{output_genome}";
	open (FH,">$path") || warn("Did not create the output file\n");
	print FH $htmlmessage;
	close FH;
	$htmlmessage .= "<pre>$htmlmessage</pre>\n\n";
    my $reportfile = Bio::KBase::utilities::add_report_file({
    	workspace_name => $params->{workspace},
    	name =>  "microbial_genome_report.$params->{output_genome}",
		path => $path,
		description => 'Microbial Annotation Report'
    });

	Bio::KBase::utilities::print_report_message({
		message => $htmlmessage,html=>1,append => 0
	});
    my $reportout = Bio::KBase::kbaseenv::create_report({
    	workspace_name => $params->{workspace},
    	report_object_name => Bio::KBase::utilities::processid().".report",
    });
	$return = {
    	workspace => $params->{workspace},
    	id => $params->{output_genome},
    	report_ref => $reportout->{"ref"},
    	report_name =>  Bio::KBase::utilities::processid().".report",
    	ws_report_id => Bio::KBase::utilities::processid().".report"
    };
    Bio::KBase::utilities::close_debug();
    #END annotate_genomes
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_genomes:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_genomes');
    }
    return($return);
}




=head2 status 

  $return = $obj->status()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module status. This is a structure including Semantic Versioning number, state and git info.

=back

=cut

sub status {
    my($return);
    #BEGIN_STATUS
    $return = {"state" => "OK", "message" => "", "version" => $VERSION,
               "git_url" => $GIT_URL, "git_commit_hash" => $GIT_COMMIT_HASH};
    #END_STATUS
    return($return);
}

=head1 TYPES



=head2 bool

=over 4



=item Description

A binary boolean


=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 genome_id

=over 4



=item Description

A string representing a genome id.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 contigset_id

=over 4



=item Description

A string representing a ContigSet id.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 workspace_name

=over 4



=item Description

A string representing a workspace name.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 AnnotateGenomeParams

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace has a value which is a string
input_genome has a value which is a RAST_SDK.genome_id
input_contigset has a value which is a RAST_SDK.contigset_id
genetic_code has a value which is an int
domain has a value which is a string
scientific_name has a value which is a string
output_genome has a value which is a string
call_features_rRNA_SEED has a value which is a RAST_SDK.bool
call_features_tRNA_trnascan has a value which is a RAST_SDK.bool
call_selenoproteins has a value which is a RAST_SDK.bool
call_pyrrolysoproteins has a value which is a RAST_SDK.bool
call_features_repeat_region_SEED has a value which is a RAST_SDK.bool
call_features_insertion_sequences has a value which is a RAST_SDK.bool
call_features_strep_suis_repeat has a value which is a RAST_SDK.bool
call_features_strep_pneumo_repeat has a value which is a RAST_SDK.bool
call_features_crispr has a value which is a RAST_SDK.bool
call_features_CDS_glimmer3 has a value which is a RAST_SDK.bool
call_features_CDS_prodigal has a value which is a RAST_SDK.bool
call_features_CDS_genemark has a value which is a RAST_SDK.bool
annotate_proteins_kmer_v2 has a value which is a RAST_SDK.bool
kmer_v1_parameters has a value which is a RAST_SDK.bool
annotate_proteins_similarity has a value which is a RAST_SDK.bool
resolve_overlapping_features has a value which is a RAST_SDK.bool
call_features_prophage_phispy has a value which is a RAST_SDK.bool
retain_old_anno_for_hypotheticals has a value which is a RAST_SDK.bool

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace has a value which is a string
input_genome has a value which is a RAST_SDK.genome_id
input_contigset has a value which is a RAST_SDK.contigset_id
genetic_code has a value which is an int
domain has a value which is a string
scientific_name has a value which is a string
output_genome has a value which is a string
call_features_rRNA_SEED has a value which is a RAST_SDK.bool
call_features_tRNA_trnascan has a value which is a RAST_SDK.bool
call_selenoproteins has a value which is a RAST_SDK.bool
call_pyrrolysoproteins has a value which is a RAST_SDK.bool
call_features_repeat_region_SEED has a value which is a RAST_SDK.bool
call_features_insertion_sequences has a value which is a RAST_SDK.bool
call_features_strep_suis_repeat has a value which is a RAST_SDK.bool
call_features_strep_pneumo_repeat has a value which is a RAST_SDK.bool
call_features_crispr has a value which is a RAST_SDK.bool
call_features_CDS_glimmer3 has a value which is a RAST_SDK.bool
call_features_CDS_prodigal has a value which is a RAST_SDK.bool
call_features_CDS_genemark has a value which is a RAST_SDK.bool
annotate_proteins_kmer_v2 has a value which is a RAST_SDK.bool
kmer_v1_parameters has a value which is a RAST_SDK.bool
annotate_proteins_similarity has a value which is a RAST_SDK.bool
resolve_overlapping_features has a value which is a RAST_SDK.bool
call_features_prophage_phispy has a value which is a RAST_SDK.bool
retain_old_anno_for_hypotheticals has a value which is a RAST_SDK.bool


=end text

=back



=head2 AnnotateGenomeResults

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 workspace_name

=over 4



=item Description

A string representing a workspace name.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 AnnotateGenomeParams

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace has a value which is a string
input_genome has a value which is a RAST_SDK.genome_id
input_contigset has a value which is a RAST_SDK.contigset_id
genetic_code has a value which is an int
domain has a value which is a string
scientific_name has a value which is a string
output_genome has a value which is a string
call_features_rRNA_SEED has a value which is a RAST_SDK.bool
call_features_tRNA_trnascan has a value which is a RAST_SDK.bool
call_selenoproteins has a value which is a RAST_SDK.bool
call_pyrrolysoproteins has a value which is a RAST_SDK.bool
call_features_repeat_region_SEED has a value which is a RAST_SDK.bool
call_features_insertion_sequences has a value which is a RAST_SDK.bool
call_features_strep_suis_repeat has a value which is a RAST_SDK.bool
call_features_strep_pneumo_repeat has a value which is a RAST_SDK.bool
call_features_crispr has a value which is a RAST_SDK.bool
call_features_CDS_glimmer3 has a value which is a RAST_SDK.bool
call_features_CDS_prodigal has a value which is a RAST_SDK.bool
call_features_CDS_genemark has a value which is a RAST_SDK.bool
annotate_proteins_kmer_v2 has a value which is a RAST_SDK.bool
kmer_v1_parameters has a value which is a RAST_SDK.bool
annotate_proteins_similarity has a value which is a RAST_SDK.bool
resolve_overlapping_features has a value which is a RAST_SDK.bool
call_features_prophage_phispy has a value which is a RAST_SDK.bool
retain_old_anno_for_hypotheticals has a value which is a RAST_SDK.bool

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace has a value which is a string
input_genome has a value which is a RAST_SDK.genome_id
input_contigset has a value which is a RAST_SDK.contigset_id
genetic_code has a value which is an int
domain has a value which is a string
scientific_name has a value which is a string
output_genome has a value which is a string
call_features_rRNA_SEED has a value which is a RAST_SDK.bool
call_features_tRNA_trnascan has a value which is a RAST_SDK.bool
call_selenoproteins has a value which is a RAST_SDK.bool
call_pyrrolysoproteins has a value which is a RAST_SDK.bool
call_features_repeat_region_SEED has a value which is a RAST_SDK.bool
call_features_insertion_sequences has a value which is a RAST_SDK.bool
call_features_strep_suis_repeat has a value which is a RAST_SDK.bool
call_features_strep_pneumo_repeat has a value which is a RAST_SDK.bool
call_features_crispr has a value which is a RAST_SDK.bool
call_features_CDS_glimmer3 has a value which is a RAST_SDK.bool
call_features_CDS_prodigal has a value which is a RAST_SDK.bool
call_features_CDS_genemark has a value which is a RAST_SDK.bool
annotate_proteins_kmer_v2 has a value which is a RAST_SDK.bool
kmer_v1_parameters has a value which is a RAST_SDK.bool
annotate_proteins_similarity has a value which is a RAST_SDK.bool
resolve_overlapping_features has a value which is a RAST_SDK.bool
call_features_prophage_phispy has a value which is a RAST_SDK.bool
retain_old_anno_for_hypotheticals has a value which is a RAST_SDK.bool


=end text

=back

call_features_strep_pneumo_repeat has a value which is a RAST_SDK.bool
call_features_crispr has a value which is a RAST_SDK.bool
call_features_CDS_glimmer3 has a value which is a RAST_SDK.bool
call_features_CDS_prodigal has a value which is a RAST_SDK.bool
call_features_CDS_genemark has a value which is a RAST_SDK.bool
annotate_proteins_kmer_v2 has a value which is a RAST_SDK.bool
kmer_v1_parameters has a value which is a RAST_SDK.bool
annotate_proteins_similarity has a value which is a RAST_SDK.bool
resolve_overlapping_features has a value which is a RAST_SDK.bool
call_features_prophage_phispy has a value which is a RAST_SDK.bool
retain_old_anno_for_hypotheticals has a value which is a RAST_SDK.bool


=end text

=back



=head2 AnnotateGenomesResults

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace has a value which is a RAST_SDK.workspace_name
report_name has a value which is a string
report_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace has a value which is a RAST_SDK.workspace_name
report_name has a value which is a string
report_ref has a value which is a string


=end text

=back



=cut

1;
