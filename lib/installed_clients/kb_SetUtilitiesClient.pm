package installed_clients::kb_SetUtilitiesClient;

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

installed_clients::kb_SetUtilitiesClient

=head1 DESCRIPTION


** A KBase module: kb_SetUtilities
**
** This module contains basic utilities for set manipulation, originally extracted
** from kb_util_dylan
**


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => installed_clients::kb_SetUtilitiesClient::RpcClient->new,
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
	my %arg_hash2 = @args;
	if (exists $arg_hash2{"token"}) {
	    $self->{token} = $arg_hash2{"token"};
	} elsif (exists $arg_hash2{"user_id"}) {
	    my $token = Bio::KBase::AuthToken->new(@args);
	    if (!$token->error_message) {
	        $self->{token} = $token->token;
	    }
	}
	
	if (exists $self->{token})
	{
	    $self->{client}->{token} = $self->{token};
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 KButil_Localize_GenomeSet

  $return = $obj->KButil_Localize_GenomeSet($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Localize_GenomeSet_Params
$return is a kb_SetUtilities.KButil_Localize_GenomeSet_Output
KButil_Localize_GenomeSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_ref has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Localize_GenomeSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Localize_GenomeSet_Params
$return is a kb_SetUtilities.KButil_Localize_GenomeSet_Output
KButil_Localize_GenomeSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_ref has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Localize_GenomeSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Localize_GenomeSet
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Localize_GenomeSet (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Localize_GenomeSet:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Localize_GenomeSet');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Localize_GenomeSet",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Localize_GenomeSet',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Localize_GenomeSet",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Localize_GenomeSet',
				       );
    }
}
 


=head2 KButil_Localize_FeatureSet

  $return = $obj->KButil_Localize_FeatureSet($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Localize_FeatureSet_Params
$return is a kb_SetUtilities.KButil_Localize_FeatureSet_Output
KButil_Localize_FeatureSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_ref has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Localize_FeatureSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Localize_FeatureSet_Params
$return is a kb_SetUtilities.KButil_Localize_FeatureSet_Output
KButil_Localize_FeatureSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_ref has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Localize_FeatureSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Localize_FeatureSet
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Localize_FeatureSet (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Localize_FeatureSet:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Localize_FeatureSet');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Localize_FeatureSet",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Localize_FeatureSet',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Localize_FeatureSet",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Localize_FeatureSet',
				       );
    }
}
 


=head2 KButil_Merge_FeatureSet_Collection

  $return = $obj->KButil_Merge_FeatureSet_Collection($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Merge_FeatureSet_Collection_Params
$return is a kb_SetUtilities.KButil_Merge_FeatureSet_Collection_Output
KButil_Merge_FeatureSet_Collection_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_refs has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Merge_FeatureSet_Collection_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Merge_FeatureSet_Collection_Params
$return is a kb_SetUtilities.KButil_Merge_FeatureSet_Collection_Output
KButil_Merge_FeatureSet_Collection_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_refs has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Merge_FeatureSet_Collection_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Merge_FeatureSet_Collection
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Merge_FeatureSet_Collection (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Merge_FeatureSet_Collection:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Merge_FeatureSet_Collection');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Merge_FeatureSet_Collection",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Merge_FeatureSet_Collection',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Merge_FeatureSet_Collection",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Merge_FeatureSet_Collection',
				       );
    }
}
 


=head2 KButil_Slice_FeatureSets_by_Genomes

  $return = $obj->KButil_Slice_FeatureSets_by_Genomes($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Slice_FeatureSets_by_Genomes_Params
$return is a kb_SetUtilities.KButil_Slice_FeatureSets_by_Genomes_Output
KButil_Slice_FeatureSets_by_Genomes_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_featureSet_refs has a value which is a kb_SetUtilities.data_obj_ref
	input_genome_refs has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Slice_FeatureSets_by_Genomes_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Slice_FeatureSets_by_Genomes_Params
$return is a kb_SetUtilities.KButil_Slice_FeatureSets_by_Genomes_Output
KButil_Slice_FeatureSets_by_Genomes_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_featureSet_refs has a value which is a kb_SetUtilities.data_obj_ref
	input_genome_refs has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Slice_FeatureSets_by_Genomes_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Slice_FeatureSets_by_Genomes
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Slice_FeatureSets_by_Genomes (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Slice_FeatureSets_by_Genomes:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Slice_FeatureSets_by_Genomes');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Slice_FeatureSets_by_Genomes",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Slice_FeatureSets_by_Genomes',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Slice_FeatureSets_by_Genomes",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Slice_FeatureSets_by_Genomes',
				       );
    }
}
 


=head2 KButil_Logical_Slice_Two_FeatureSets

  $return = $obj->KButil_Logical_Slice_Two_FeatureSets($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Logical_Slice_Two_FeatureSets_Params
$return is a kb_SetUtilities.KButil_Logical_Slice_Two_FeatureSets_Output
KButil_Logical_Slice_Two_FeatureSets_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_featureSet_ref_A has a value which is a kb_SetUtilities.data_obj_ref
	input_featureSet_ref_B has a value which is a kb_SetUtilities.data_obj_ref
	operator has a value which is a string
	desc has a value which is a string
	output_name has a value which is a kb_SetUtilities.data_obj_name
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Logical_Slice_Two_FeatureSets_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Logical_Slice_Two_FeatureSets_Params
$return is a kb_SetUtilities.KButil_Logical_Slice_Two_FeatureSets_Output
KButil_Logical_Slice_Two_FeatureSets_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_featureSet_ref_A has a value which is a kb_SetUtilities.data_obj_ref
	input_featureSet_ref_B has a value which is a kb_SetUtilities.data_obj_ref
	operator has a value which is a string
	desc has a value which is a string
	output_name has a value which is a kb_SetUtilities.data_obj_name
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Logical_Slice_Two_FeatureSets_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Logical_Slice_Two_FeatureSets
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Logical_Slice_Two_FeatureSets (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Logical_Slice_Two_FeatureSets:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Logical_Slice_Two_FeatureSets');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Logical_Slice_Two_FeatureSets",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Logical_Slice_Two_FeatureSets',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Logical_Slice_Two_FeatureSets",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Logical_Slice_Two_FeatureSets',
				       );
    }
}
 


=head2 KButil_Merge_GenomeSets

  $return = $obj->KButil_Merge_GenomeSets($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Merge_GenomeSets_Params
$return is a kb_SetUtilities.KButil_Merge_GenomeSets_Output
KButil_Merge_GenomeSets_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_refs has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Merge_GenomeSets_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Merge_GenomeSets_Params
$return is a kb_SetUtilities.KButil_Merge_GenomeSets_Output
KButil_Merge_GenomeSets_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_refs has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Merge_GenomeSets_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Merge_GenomeSets
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Merge_GenomeSets (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Merge_GenomeSets:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Merge_GenomeSets');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Merge_GenomeSets",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Merge_GenomeSets',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Merge_GenomeSets",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Merge_GenomeSets',
				       );
    }
}
 


=head2 KButil_Build_GenomeSet

  $return = $obj->KButil_Build_GenomeSet($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Build_GenomeSet_Params
$return is a kb_SetUtilities.KButil_Build_GenomeSet_Output
KButil_Build_GenomeSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_refs has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Build_GenomeSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Build_GenomeSet_Params
$return is a kb_SetUtilities.KButil_Build_GenomeSet_Output
KButil_Build_GenomeSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_refs has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Build_GenomeSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Build_GenomeSet
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Build_GenomeSet (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Build_GenomeSet:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Build_GenomeSet');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Build_GenomeSet",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Build_GenomeSet',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Build_GenomeSet",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Build_GenomeSet',
				       );
    }
}
 


=head2 KButil_Build_GenomeSet_from_FeatureSet

  $return = $obj->KButil_Build_GenomeSet_from_FeatureSet($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Build_GenomeSet_from_FeatureSet_Params
$return is a kb_SetUtilities.KButil_Build_GenomeSet_from_FeatureSet_Output
KButil_Build_GenomeSet_from_FeatureSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_ref has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Build_GenomeSet_from_FeatureSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Build_GenomeSet_from_FeatureSet_Params
$return is a kb_SetUtilities.KButil_Build_GenomeSet_from_FeatureSet_Output
KButil_Build_GenomeSet_from_FeatureSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_ref has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Build_GenomeSet_from_FeatureSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Build_GenomeSet_from_FeatureSet
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Build_GenomeSet_from_FeatureSet (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Build_GenomeSet_from_FeatureSet:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Build_GenomeSet_from_FeatureSet');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Build_GenomeSet_from_FeatureSet",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Build_GenomeSet_from_FeatureSet',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Build_GenomeSet_from_FeatureSet",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Build_GenomeSet_from_FeatureSet',
				       );
    }
}
 


=head2 KButil_Add_Genomes_to_GenomeSet

  $return = $obj->KButil_Add_Genomes_to_GenomeSet($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Add_Genomes_to_GenomeSet_Params
$return is a kb_SetUtilities.KButil_Add_Genomes_to_GenomeSet_Output
KButil_Add_Genomes_to_GenomeSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_genome_refs has a value which is a reference to a list where each element is a kb_SetUtilities.data_obj_ref
	input_genomeset_ref has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Add_Genomes_to_GenomeSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Add_Genomes_to_GenomeSet_Params
$return is a kb_SetUtilities.KButil_Add_Genomes_to_GenomeSet_Output
KButil_Add_Genomes_to_GenomeSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_genome_refs has a value which is a reference to a list where each element is a kb_SetUtilities.data_obj_ref
	input_genomeset_ref has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Add_Genomes_to_GenomeSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Add_Genomes_to_GenomeSet
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Add_Genomes_to_GenomeSet (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Add_Genomes_to_GenomeSet:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Add_Genomes_to_GenomeSet');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Add_Genomes_to_GenomeSet",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Add_Genomes_to_GenomeSet',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Add_Genomes_to_GenomeSet",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Add_Genomes_to_GenomeSet',
				       );
    }
}
 


=head2 KButil_Remove_Genomes_from_GenomeSet

  $return = $obj->KButil_Remove_Genomes_from_GenomeSet($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Remove_Genomes_from_GenomeSet_Params
$return is a kb_SetUtilities.KButil_Remove_Genomes_from_GenomeSet_Output
KButil_Remove_Genomes_from_GenomeSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_genome_refs has a value which is a reference to a list where each element is a kb_SetUtilities.data_obj_ref
	nonlocal_genome_names has a value which is a reference to a list where each element is a kb_SetUtilities.data_obj_name
	input_genomeset_ref has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Remove_Genomes_from_GenomeSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Remove_Genomes_from_GenomeSet_Params
$return is a kb_SetUtilities.KButil_Remove_Genomes_from_GenomeSet_Output
KButil_Remove_Genomes_from_GenomeSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_genome_refs has a value which is a reference to a list where each element is a kb_SetUtilities.data_obj_ref
	nonlocal_genome_names has a value which is a reference to a list where each element is a kb_SetUtilities.data_obj_name
	input_genomeset_ref has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Remove_Genomes_from_GenomeSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Remove_Genomes_from_GenomeSet
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Remove_Genomes_from_GenomeSet (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Remove_Genomes_from_GenomeSet:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Remove_Genomes_from_GenomeSet');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Remove_Genomes_from_GenomeSet",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Remove_Genomes_from_GenomeSet',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Remove_Genomes_from_GenomeSet",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Remove_Genomes_from_GenomeSet',
				       );
    }
}
 


=head2 KButil_Build_ReadsSet

  $return = $obj->KButil_Build_ReadsSet($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Build_ReadsSet_Params
$return is a kb_SetUtilities.KButil_Build_ReadsSet_Output
KButil_Build_ReadsSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_refs has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Build_ReadsSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Build_ReadsSet_Params
$return is a kb_SetUtilities.KButil_Build_ReadsSet_Output
KButil_Build_ReadsSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_refs has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Build_ReadsSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Build_ReadsSet
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Build_ReadsSet (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Build_ReadsSet:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Build_ReadsSet');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Build_ReadsSet",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Build_ReadsSet',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Build_ReadsSet",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Build_ReadsSet',
				       );
    }
}
 


=head2 KButil_Merge_MultipleReadsSets_to_OneReadsSet

  $return = $obj->KButil_Merge_MultipleReadsSets_to_OneReadsSet($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Merge_MultipleReadsSets_to_OneReadsSet_Params
$return is a kb_SetUtilities.KButil_Merge_MultipleReadsSets_to_OneReadsSet_Output
KButil_Merge_MultipleReadsSets_to_OneReadsSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_refs has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Merge_MultipleReadsSets_to_OneReadsSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Merge_MultipleReadsSets_to_OneReadsSet_Params
$return is a kb_SetUtilities.KButil_Merge_MultipleReadsSets_to_OneReadsSet_Output
KButil_Merge_MultipleReadsSets_to_OneReadsSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_refs has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Merge_MultipleReadsSets_to_OneReadsSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Merge_MultipleReadsSets_to_OneReadsSet
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Merge_MultipleReadsSets_to_OneReadsSet (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Merge_MultipleReadsSets_to_OneReadsSet:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Merge_MultipleReadsSets_to_OneReadsSet');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Merge_MultipleReadsSets_to_OneReadsSet",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Merge_MultipleReadsSets_to_OneReadsSet',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Merge_MultipleReadsSets_to_OneReadsSet",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Merge_MultipleReadsSets_to_OneReadsSet',
				       );
    }
}
 


=head2 KButil_Build_AssemblySet

  $return = $obj->KButil_Build_AssemblySet($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Build_AssemblySet_Params
$return is a kb_SetUtilities.KButil_Build_AssemblySet_Output
KButil_Build_AssemblySet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_refs has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Build_AssemblySet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Build_AssemblySet_Params
$return is a kb_SetUtilities.KButil_Build_AssemblySet_Output
KButil_Build_AssemblySet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	input_refs has a value which is a kb_SetUtilities.data_obj_ref
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_ref is a string
data_obj_name is a string
KButil_Build_AssemblySet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=item Description



=back

=cut

 sub KButil_Build_AssemblySet
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Build_AssemblySet (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Build_AssemblySet:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Build_AssemblySet');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Build_AssemblySet",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Build_AssemblySet',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Build_AssemblySet",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Build_AssemblySet',
				       );
    }
}
 


=head2 KButil_Batch_Create_ReadsSet

  $return = $obj->KButil_Batch_Create_ReadsSet($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Batch_Create_ReadsSet_Params
$return is a kb_SetUtilities.KButil_Batch_Create_ReadsSet_Output
KButil_Batch_Create_ReadsSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	name_pattern has a value which is a string
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_name is a string
KButil_Batch_Create_ReadsSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref
data_obj_ref is a string

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Batch_Create_ReadsSet_Params
$return is a kb_SetUtilities.KButil_Batch_Create_ReadsSet_Output
KButil_Batch_Create_ReadsSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	name_pattern has a value which is a string
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_name is a string
KButil_Batch_Create_ReadsSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref
data_obj_ref is a string


=end text

=item Description



=back

=cut

 sub KButil_Batch_Create_ReadsSet
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Batch_Create_ReadsSet (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Batch_Create_ReadsSet:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Batch_Create_ReadsSet');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Batch_Create_ReadsSet",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Batch_Create_ReadsSet',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Batch_Create_ReadsSet",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Batch_Create_ReadsSet',
				       );
    }
}
 


=head2 KButil_Batch_Create_AssemblySet

  $return = $obj->KButil_Batch_Create_AssemblySet($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Batch_Create_AssemblySet_Params
$return is a kb_SetUtilities.KButil_Batch_Create_AssemblySet_Output
KButil_Batch_Create_AssemblySet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	name_pattern has a value which is a string
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_name is a string
KButil_Batch_Create_AssemblySet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref
data_obj_ref is a string

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Batch_Create_AssemblySet_Params
$return is a kb_SetUtilities.KButil_Batch_Create_AssemblySet_Output
KButil_Batch_Create_AssemblySet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	name_pattern has a value which is a string
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_name is a string
KButil_Batch_Create_AssemblySet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref
data_obj_ref is a string


=end text

=item Description



=back

=cut

 sub KButil_Batch_Create_AssemblySet
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Batch_Create_AssemblySet (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Batch_Create_AssemblySet:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Batch_Create_AssemblySet');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Batch_Create_AssemblySet",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Batch_Create_AssemblySet',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Batch_Create_AssemblySet",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Batch_Create_AssemblySet',
				       );
    }
}
 


=head2 KButil_Batch_Create_GenomeSet

  $return = $obj->KButil_Batch_Create_GenomeSet($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_SetUtilities.KButil_Batch_Create_GenomeSet_Params
$return is a kb_SetUtilities.KButil_Batch_Create_GenomeSet_Output
KButil_Batch_Create_GenomeSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	name_pattern has a value which is a string
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_name is a string
KButil_Batch_Create_GenomeSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref
data_obj_ref is a string

</pre>

=end html

=begin text

$params is a kb_SetUtilities.KButil_Batch_Create_GenomeSet_Params
$return is a kb_SetUtilities.KButil_Batch_Create_GenomeSet_Output
KButil_Batch_Create_GenomeSet_Params is a reference to a hash where the following keys are defined:
	workspace_name has a value which is a kb_SetUtilities.workspace_name
	name_pattern has a value which is a string
	output_name has a value which is a kb_SetUtilities.data_obj_name
	desc has a value which is a string
workspace_name is a string
data_obj_name is a string
KButil_Batch_Create_GenomeSet_Output is a reference to a hash where the following keys are defined:
	report_name has a value which is a kb_SetUtilities.data_obj_name
	report_ref has a value which is a kb_SetUtilities.data_obj_ref
data_obj_ref is a string


=end text

=item Description



=back

=cut

 sub KButil_Batch_Create_GenomeSet
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function KButil_Batch_Create_GenomeSet (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to KButil_Batch_Create_GenomeSet:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'KButil_Batch_Create_GenomeSet');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_SetUtilities.KButil_Batch_Create_GenomeSet",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'KButil_Batch_Create_GenomeSet',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method KButil_Batch_Create_GenomeSet",
					    status_line => $self->{client}->status_line,
					    method_name => 'KButil_Batch_Create_GenomeSet',
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
        method => "kb_SetUtilities.status",
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
        method => "kb_SetUtilities.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'KButil_Batch_Create_GenomeSet',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method KButil_Batch_Create_GenomeSet",
            status_line => $self->{client}->status_line,
            method_name => 'KButil_Batch_Create_GenomeSet',
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
        warn "New client version available for installed_clients::kb_SetUtilitiesClient\n";
    }
    if ($sMajor == 0) {
        warn "installed_clients::kb_SetUtilitiesClient version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 workspace_name

=over 4



=item Description

** The workspace object refs are of form:
**
**    objects = ws.get_objects([{'ref': params['workspace_id']+'/'+params['obj_name']}])
**
** "ref" means the entire name combining the workspace id and the object name
** "id" is a numerical identifier of the workspace or object, and should just be used for workspace
** "name" is a string identifier of a workspace or object.  This is received from Narrative.


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



=head2 sequence

=over 4



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



=head2 data_obj_name

=over 4



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



=head2 data_obj_ref

=over 4



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



=head2 bool

=over 4



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



=head2 KButil_Localize_GenomeSet_Params

=over 4



=item Description

KButil_Localize_GenomeSet()
**
**  Method for creating Genome Set with all local Genomes


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_ref has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_ref has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name


=end text

=back



=head2 KButil_Localize_GenomeSet_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Localize_FeatureSet_Params

=over 4



=item Description

KButil_Localize_FeatureSet()
**
**  Method for creating Feature Set with all local Genomes


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_ref has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_ref has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name


=end text

=back



=head2 KButil_Localize_FeatureSet_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Merge_FeatureSet_Collection_Params

=over 4



=item Description

KButil_Merge_FeatureSet_Collection()
**
**  Method for merging FeatureSets


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_refs has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_refs has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Merge_FeatureSet_Collection_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Slice_FeatureSets_by_Genomes_Params

=over 4



=item Description

KButil_Slice_FeatureSets_by_Genomes()
**
**  Method for Slicing a FeatureSet or FeatureSets by a Genome, Genomes, or GenomeSet


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_featureSet_refs has a value which is a kb_SetUtilities.data_obj_ref
input_genome_refs has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_featureSet_refs has a value which is a kb_SetUtilities.data_obj_ref
input_genome_refs has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Slice_FeatureSets_by_Genomes_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Logical_Slice_Two_FeatureSets_Params

=over 4



=item Description

KButil_Logical_Slice_Two_FeatureSets()
**
**  Method for Slicing Two FeatureSets by Venn overlap


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_featureSet_ref_A has a value which is a kb_SetUtilities.data_obj_ref
input_featureSet_ref_B has a value which is a kb_SetUtilities.data_obj_ref
operator has a value which is a string
desc has a value which is a string
output_name has a value which is a kb_SetUtilities.data_obj_name

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_featureSet_ref_A has a value which is a kb_SetUtilities.data_obj_ref
input_featureSet_ref_B has a value which is a kb_SetUtilities.data_obj_ref
operator has a value which is a string
desc has a value which is a string
output_name has a value which is a kb_SetUtilities.data_obj_name


=end text

=back



=head2 KButil_Logical_Slice_Two_FeatureSets_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Merge_GenomeSets_Params

=over 4



=item Description

KButil_Merge_GenomeSets()
**
**  Method for merging GenomeSets


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_refs has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_refs has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Merge_GenomeSets_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Build_GenomeSet_Params

=over 4



=item Description

KButil_Build_GenomeSet()
**
**  Method for creating a GenomeSet


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_refs has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_refs has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Build_GenomeSet_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Build_GenomeSet_from_FeatureSet_Params

=over 4



=item Description

KButil_Build_GenomeSet_from_FeatureSet()
**
**  Method for obtaining a GenomeSet from a FeatureSet


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_ref has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_ref has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Build_GenomeSet_from_FeatureSet_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Add_Genomes_to_GenomeSet_Params

=over 4



=item Description

KButil_Add_Genomes_to_GenomeSet()
**
**  Method for adding a Genome to a GenomeSet


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_genome_refs has a value which is a reference to a list where each element is a kb_SetUtilities.data_obj_ref
input_genomeset_ref has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_genome_refs has a value which is a reference to a list where each element is a kb_SetUtilities.data_obj_ref
input_genomeset_ref has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Add_Genomes_to_GenomeSet_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Remove_Genomes_from_GenomeSet_Params

=over 4



=item Description

KButil_Remove_Genomes_from_GenomeSet()
**
**  Method for removing Genomes from a GenomeSet


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_genome_refs has a value which is a reference to a list where each element is a kb_SetUtilities.data_obj_ref
nonlocal_genome_names has a value which is a reference to a list where each element is a kb_SetUtilities.data_obj_name
input_genomeset_ref has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_genome_refs has a value which is a reference to a list where each element is a kb_SetUtilities.data_obj_ref
nonlocal_genome_names has a value which is a reference to a list where each element is a kb_SetUtilities.data_obj_name
input_genomeset_ref has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Remove_Genomes_from_GenomeSet_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Build_ReadsSet_Params

=over 4



=item Description

KButil_Build_ReadsSet()
**
**  Method for creating a ReadsSet


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_refs has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_refs has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Build_ReadsSet_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Merge_MultipleReadsSets_to_OneReadsSet_Params

=over 4



=item Description

KButil_Merge_MultipleReadsSets_to_OneReadsSet()
**
**  Method for merging multiple ReadsSets into one ReadsSet


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_refs has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_refs has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Merge_MultipleReadsSets_to_OneReadsSet_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Build_AssemblySet_Params

=over 4



=item Description

KButil_Build_AssemblySet()
**
**  Method for creating an AssemblySet


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_refs has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_refs has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Build_AssemblySet_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Batch_Create_ReadsSet_Params

=over 4



=item Description

KButil_Batch_Create_ReadsSet()
**
**  Method for creating a ReadsSet without specifying individual objects


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
name_pattern has a value which is a string
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
name_pattern has a value which is a string
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Batch_Create_ReadsSet_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Batch_Create_AssemblySet_Params

=over 4



=item Description

KButil_Batch_Create_AssemblySet()
**
**  Method for creating an AssemblySet without specifying individual objects


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
name_pattern has a value which is a string
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
name_pattern has a value which is a string
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Batch_Create_AssemblySet_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=head2 KButil_Batch_Create_GenomeSet_Params

=over 4



=item Description

KButil_Batch_Create_GenomeSet()
**
**  Method for creating a GenomeSet without specifying individual objects


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
name_pattern has a value which is a string
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
name_pattern has a value which is a string
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string


=end text

=back



=head2 KButil_Batch_Create_GenomeSet_Output

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a kb_SetUtilities.data_obj_name
report_ref has a value which is a kb_SetUtilities.data_obj_ref


=end text

=back



=cut

package installed_clients::kb_SetUtilitiesClient::RpcClient;
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
