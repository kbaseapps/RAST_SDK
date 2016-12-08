package RAST_SDK::RAST_SDKImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = "0.1.0";

=head1 NAME

RAST_SDK

=head1 DESCRIPTION

The SDK version of the KBaase Genome Annotation Service.
This wraps genome_annotation which is based off of the SEED annotations.

=cut

#BEGIN_HEADER
use Bio::KBase::AuthToken;
use Bio::KBase::workspace::Client;
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
	return $params;
}

sub util_version {
	my ($self) = @_;
	return "1";
}

sub util_log {
	my($self,$message) = @_;
	print $message."\n";
}

sub util_get_genome {
	my ($self,$workspace,$genomeid) = @_;
	my $output = Bio::KBase::kbaseenv::ga_client()->get_genome_v1({
		genomes => [{
			"ref" => $workspace."/".$genomeid
		}],
		ignore_errors => 1,
		no_data => 0,
		no_metadata => 1
	});
	return $output->{genomes}->[0]->{data};
}

sub util_get_contigs {
	my ($self,$workspace,$objid) = @_;
	my $info = Bio::KBase::kbaseenv::get_object_info([
		Bio::KBase::kbaseenv::configure_ws_id($workspace,$objid)
	],0);
	my $obj;
	if ($info->[0]->[2] =~ /Assembly/) {
		my $output = Bio::KBase::kbaseenv::ac_client()->get_assembly_as_fasta({
			"ref" => $workspace."/".$objid
		});
		my $fasta = "";
		open(my $fh, "<", $output->{path}) || die "Could not find file:".$output->{path};
		while (my $line = <$fh>) {
			$fasta .= $line;
		}
		close($fh);
		#Bio::KBase::utilities::debug(Data::Dumper->Dump([$fasta]));
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
		for (my $i=0; $i < @{$array}; $i++) {
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
		$obj->{md5} = Digest::MD5::md5_hex($str);
		$obj->{_kbasetype} = "Assembly";
	} else {
		$obj = Bio::KBase::kbaseenv::get_objects([
			Bio::KBase::kbaseenv::configure_ws_id($workspace,$objid)
		]);
		$obj = $obj->[0]->{data};
		$obj->{_kbasetype} = "ContigSet";
		$obj->{_reference} = $info->[0]->[6]."/".$info->[0]->[0]."/".$info->[0]->[4];
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

sub annotate {
	my ($self,$parameters) = @_;	
	my $oldfunchash = {};
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
	if (defined($parameters->{input_genome})) {
		$inputgenome = $self->util_get_genome($parameters->{workspace},$parameters->{input_genome});
		for (my $i=0; $i < @{$inputgenome->{features}}; $i++) {
			if (lc($inputgenome->{features}->[$i]->{type}) eq "cds" || lc($inputgenome->{features}->[$i]->{type}) eq "peg") {
				$oldfunchash->{$inputgenome->{features}->[$i]->{id}} = $inputgenome->{features}->[$i]->{function};
				$inputgenome->{features}->[$i]->{function} = "hypothetical protein";
			}
		}
		my $contigref;
		if (defined($inputgenome->{contigset_ref})) {
			$contigref = $inputgenome->{contigset_ref};
		} elsif (defined($inputgenome->{assembly_ref})) {
			$contigref = $inputgenome->{assembly_ref};
		}
		if ($contigref =~ m/^([^\/]+)\/([^\/]+)/) {
			$contigobj = $self->util_get_contigs($1,$2);
		}
		$parameters->{genetic_code} = $inputgenome->{genetic_code};
		$parameters->{domain} = $inputgenome->{domain};
		$parameters->{scientific_name} = $inputgenome->{scientific_name};
	} elsif (defined($parameters->{input_contigset})) {
		$contigobj = $self->util_get_contigs($parameters->{workspace},$parameters->{input_contigset});	
	} else {
		Bio::KBase::utilities::error("Neither contigs nor genome specified!");
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
			$message = "The RAST algorithm was applied to annotating a genome sequence comprised of ".$count." contigs containing ".$size." nucleotides. No initial gene calls were provided.";
		} else {
			$message = "The RAST algorithm was applied to annotating an existing genome: ".$parameters->{scientific_name}.". The sequence for this genome is comprised of ".$count." contigs containing ".$size." nucleotides. Input genome has ".@{$inputgenome->{features}}." existing features.";
		}		
	} else {
		$message = "The RAST algorithm was applied to annotating an existing genome: ".$parameters->{scientific_name}.". No DNA sequence was provided for this genome, therefore new genes cannot be called. We can only functionally annotate the ".@{$inputgenome->{features}}." existing features.";
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
	}
	if (defined($parameters->{call_pyrrolysoproteins}) && $parameters->{call_pyrrolysoproteins} == 1)	{
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
	if (defined($parameters->{call_features_insertion_sequences}) && $parameters->{call_features_insertion_sequences} == 1)	{
		if (length($extragenecalls) == 0) {
			$extragenecalls = "A scan was conducted for the following additional feature types: ";
		} else {
			$extragenecalls .= "; ";
		}
		$extragenecalls .= "insertion sequences";
		push(@{$workflow->{stages}},{name => "call_features_insertion_sequences"});
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
	}
	my $genecalls = "";
	if (defined($parameters->{call_features_CDS_glimmer3}) && $parameters->{call_features_CDS_glimmer3} == 1)	{
		if (@{$inputgenome->{features}} > 0) {
			$inputgenome->{features} = [];
			$message .= " The existing gene features were cleared due to selection of gene calling with Glimmer3, Prodigal, or Genmark.";
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
			Bio::KBase::utilities::error("Cannot call genes on genome with no contigs!");
		}
	}
	if (defined($parameters->{call_features_CDS_prodigal}) && $parameters->{call_features_CDS_prodigal} == 1)	{
		if (@{$inputgenome->{features}} > 0) {
			$inputgenome->{features} = [];
			$message .= " The existing gene features were cleared due to selection of gene calling with Glimmer3, Prodigal, or Genmark.";
		}
		if (length($genecalls) == 0) {
			$genecalls = "Standard gene features were called using: ";
		} else {
			$genecalls .= "; ";
		}
		$genecalls .= "prodigal";
		push(@{$workflow->{stages}},{name => "call_features_CDS_prodigal"});
		if (!defined($contigobj)) {
			Bio::KBase::utilities::error("Cannot call genes on genome with no contigs!");
		}
	}
#	if (defined($parameters->{call_features_CDS_genemark}) && $parameters->{call_features_CDS_genemark} == 1)	{
#		if (@{$inputgenome->{features}} > 0) {
#			$inputgenome->{features} = [];
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
	if (defined($parameters->{find_close_neighbors}) && $parameters->{find_close_neighbors} == 1)	{
		push(@{$workflow->{stages}},{name => "find_close_neighbors"});
	}
	if (defined($parameters->{call_features_prophage_phispy}) && $parameters->{call_features_prophage_phispy} == 1)	{
		push(@{$workflow->{stages}},{name => "call_features_prophage_phispy"});
	}
	if (length($genecalls) > 0) {
		$message .= " ".$genecalls.".";
	}
	if (length($genecalls) > 0) {
		$message .= " ".$extragenecalls.".";
	}
	if (length($annomessage) > 0) {
		$message .= " ".$annomessage.".";
	}
	my $genome = $inputgenome;
	#Bio::KBase::utilities::debug(Data::Dumper->Dump([$inputgenome]));
	my $genehash;
	if (defined($genome->{features})) {
		for (my $i=0; $i < @{$genome->{features}}; $i++) {
			$genehash->{$genome->{features}->[$i]->{id}}->{$genome->{features}->[$i]->{function}} = 1;
		}
	}
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
	if ( defined($contigobj->{contigs}) && scalar(@{$contigobj->{contigs}})>0 ) {
		$genome->{num_contigs} = @{$contigobj->{contigs}};
		$genome->{md5} = $contigobj->{md5};
	}
	#Getting the seed ontology dictionary
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
	my $functionchanges = 0;
	my $proteins = 0;
	my $others = 0;
	my $seedfunctions;
	my $genomefunchash;
	my $seedfunchash;
	if (defined($genome->{features})) {
		for (my $i=0; $i < @{$genome->{features}}; $i++) {
			my $ftr = $genome->{features}->[$i];
			if (defined($genehash) && !defined($genehash->{$ftr->{id}})) {
				$newftrs++;
				if (!defined($genehash->{$ftr->{id}}->{$ftr->{function}})) {
					$functionchanges++;
				}
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
						if (!defined($ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}})) {
							$ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}} = {
								 evidence => [],
								 id => $funchash->{$rolename}->{id},
								 term_name => $funchash->{$rolename}->{name},
								 ontology_ref => $output->[0]->{info}->[6]."/".$output->[0]->{info}->[0]."/".$output->[0]->{info}->[4],
								 term_lineage => [],
							};
						}
						my $found = 0;
						for (my $k=0; $k < @{$ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}}; $k++) {
							if ($ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}->[$k]->{method} eq Bio::KBase::utilities::method()) {
								$ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}->[$k]->{timestamp} = Bio::KBase::utilities::timestamp();
								$ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}->[$k]->{method_version} = $self->util_version();
								$found = 1;
								last;
							}
						}
						if ($found == 0) {
							push(@{$ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}},{
								method => Bio::KBase::utilities::method(),
								method_version => $self->util_version(),
								timestamp => Bio::KBase::utilities::timestamp()
							});
						}
					}
				}
			}
		}
	}
	if (defined($genehash)) {
		$message .= ". In addition to the original ".keys(%{$genehash})." features, ".$newftrs." new features were called";
		$message .= ". Of the original features, ".$functionchanges." were re-annotated by RAST with new functions.";
	}
	$message .= ". Overall, a total of ".$seedfunctions." genes are now annotated with ".keys(%{$genomefunchash})." distinct functions. Of these functions, ".keys(%{$seedfunchash})." are a match for the SEED annotation ontology.";
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
	Bio::KBase::utilities::debug(Bio::KBase::utilities::to_json($contigobj,1));
	Bio::KBase::utilities::debug(Bio::KBase::utilities::to_json($genome,1));
	my $gaout = Bio::KBase::kbaseenv::ga_client()->save_one_genome_v1({
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
		message => $message,
		append => 0,
		html => 0
	});
	Bio::KBase::utilities::print_report_message({
		message => "<p>".$message."</p>",
		append => 0,
		html => 1
	});
	return {"ref" => $gaout->{info}->[6]."/".$gaout->{info}->[0]."/".$gaout->{info}->[4]};
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
$params is an UnspecifiedObject, which can hold any non-null object
$return is a RAST_SDK.AnnotateGenomeResults
AnnotateGenomeResults is a reference to a hash where the following keys are defined:
	workspace has a value which is a RAST_SDK.workspace_name
	id has a value which is a string
workspace_name is a string

</pre>

=end html

=begin text

$params is an UnspecifiedObject, which can hold any non-null object
$return is a RAST_SDK.AnnotateGenomeResults
AnnotateGenomeResults is a reference to a hash where the following keys are defined:
	workspace has a value which is a RAST_SDK.workspace_name
	id has a value which is a string
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
    (defined $params) or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
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
	    call_features_insertion_sequences => 1,
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
	    find_close_neighbors => 0,
	    call_features_prophage_phispy => 1,
	    retain_old_anno_for_hypotheticals => 1
	});
    my $output = $self->annotate($params);
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




=head2 version 

  $return = $obj->version()

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

Return the module version. This is a Semantic Versioning number.

=back

=cut

sub version {
    return $VERSION;
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
find_close_neighbors has a value which is a RAST_SDK.bool
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
find_close_neighbors has a value which is a RAST_SDK.bool
call_features_prophage_phispy has a value which is a RAST_SDK.bool
retain_old_anno_for_hypotheticals has a value which is a RAST_SDK.bool


=end text

=back



=head2 AnnotateGenomeResults

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace has a value which is a RAST_SDK.workspace_name
id has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace has a value which is a RAST_SDK.workspace_name
id has a value which is a string


=end text

=back



=cut

1;
