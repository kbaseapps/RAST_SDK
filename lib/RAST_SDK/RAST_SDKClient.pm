package RAST_SDK::RAST_SDKClient;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

RAST_SDK::RAST_SDKClient

=head1 DESCRIPTION


The SDK version of the KBaase Genome Annotation Service.
This wraps genome_annotation which is based off of the SEED annotations.


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => RAST_SDK::RAST_SDKClient::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my $token = Bio::KBase::AuthToken->new(@args);
	
	if (!$token->error_message)
	{
	    $self->{token} = $token->token;
	    $self->{client}->{token} = $token->token;
	}
        else
        {
	    #
	    # All methods in this module require authentication. In this case, if we
	    # don't have a token, we can't continue.
	    #
	    die "Authentication failed: " . $token->error_message;
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




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
	find_close_neighbors has a value which is a RAST_SDK.bool
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
	find_close_neighbors has a value which is a RAST_SDK.bool
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
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function annotate_genome (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'annotate_genome');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "RAST_SDK.annotate_genome",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'annotate_genome',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method annotate_genome",
					    status_line => $self->{client}->status_line,
					    method_name => 'annotate_genome',
				       );
    }
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
	find_close_neighbors has a value which is a RAST_SDK.bool
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
	find_close_neighbors has a value which is a RAST_SDK.bool
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
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function annotate_genomes (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_genomes:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'annotate_genomes');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "RAST_SDK.annotate_genomes",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'annotate_genomes',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method annotate_genomes",
					    status_line => $self->{client}->status_line,
					    method_name => 'annotate_genomes',
				       );
    }
}
 
  
sub status
{
    my($self, @args) = @_;
    if ((my $n = @args) != 0) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function status (received $n, expecting 0)");
    }
    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
        method => "RAST_SDK.status",
        params => \@args,
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => 'status',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
                          );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method status",
                        status_line => $self->{client}->status_line,
                        method_name => 'status',
                       );
    }
}
   

sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "RAST_SDK.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'annotate_genomes',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method annotate_genomes",
            status_line => $self->{client}->status_line,
            method_name => 'annotate_genomes',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for RAST_SDK::RAST_SDKClient\n";
    }
    if ($sMajor == 0) {
        warn "RAST_SDK::RAST_SDKClient version is $svr_version. API subject to change.\n";
    }
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
report_name has a value which is a string
report_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace has a value which is a RAST_SDK.workspace_name
id has a value which is a string
report_name has a value which is a string
report_ref has a value which is a string


=end text

=back



=head2 GenomeParams

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
input_contigset has a value which is a RAST_SDK.contigset_id
input_genome has a value which is a RAST_SDK.genome_id
output_genome has a value which is a RAST_SDK.genome_id
genetic_code has a value which is an int
domain has a value which is a string
scientific_name has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
input_contigset has a value which is a RAST_SDK.contigset_id
input_genome has a value which is a RAST_SDK.genome_id
output_genome has a value which is a RAST_SDK.genome_id
genetic_code has a value which is an int
domain has a value which is a string
scientific_name has a value which is a string


=end text

=back



=head2 AnnotateGenomesParams

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
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
find_close_neighbors has a value which is a RAST_SDK.bool
call_features_prophage_phispy has a value which is a RAST_SDK.bool
retain_old_anno_for_hypotheticals has a value which is a RAST_SDK.bool

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
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
find_close_neighbors has a value which is a RAST_SDK.bool
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

package RAST_SDK::RAST_SDKClient::RpcClient;
use base 'JSON::RPC::Client';
use POSIX;
use strict;

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
