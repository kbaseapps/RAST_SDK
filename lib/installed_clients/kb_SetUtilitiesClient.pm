package installed_clients::kb_SetUtilitiesClient;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
use Time::HiRes;
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
    my %arg_hash = @args;
    $self->{async_job_check_time} = 0.1;
    if (exists $arg_hash{"async_job_check_time_ms"}) {
        $self->{async_job_check_time} = $arg_hash{"async_job_check_time_ms"} / 1000.0;
    }
    $self->{async_job_check_time_scale_percent} = 150;
    if (exists $arg_hash{"async_job_check_time_scale_percent"}) {
        $self->{async_job_check_time_scale_percent} = $arg_hash{"async_job_check_time_scale_percent"};
    }
    $self->{async_job_check_max_time} = 300;  # 5 minutes
    if (exists $arg_hash{"async_job_check_max_time_ms"}) {
        $self->{async_job_check_max_time} = $arg_hash{"async_job_check_max_time_ms"} / 1000.0;
    }
    my $service_version = 'release';
    if (exists $arg_hash{"service_version"}) {
        $service_version = $arg_hash{"service_version"};
    }
    $self->{service_version} = $service_version;

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

sub _check_job {
    my($self, @args) = @_;
# Authentication: ${method.authentication}
    if ((my $n = @args) != 1) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function _check_job (received $n, expecting 1)");
    }
    {
        my($job_id) = @args;
        my @_bad_arguments;
        (!ref($job_id)) or push(@_bad_arguments, "Invalid type for argument 0 \"job_id\" (it should be a string)");
        if (@_bad_arguments) {
            my $msg = "Invalid arguments passed to _check_job:\n" . join("", map { "\t$_\n" } @_bad_arguments);
            Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
                                   method_name => '_check_job');
        }
    }
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_SetUtilities._check_job",
        params => \@args});
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => '_check_job',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
                          );
        } else {
            return $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method _check_job",
                        status_line => $self->{client}->status_line,
                        method_name => '_check_job');
    }
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
    my $job_id = $self->_KButil_Localize_GenomeSet_submit(@args);
    my $async_job_check_time = $self->{async_job_check_time};
    while (1) {
        Time::HiRes::sleep($async_job_check_time);
        $async_job_check_time *= $self->{async_job_check_time_scale_percent} / 100.0;
        if ($async_job_check_time > $self->{async_job_check_max_time}) {
            $async_job_check_time = $self->{async_job_check_max_time};
        }
        my $job_state_ref = $self->_check_job($job_id);
        if ($job_state_ref->{"finished"} != 0) {
            if (!exists $job_state_ref->{"result"}) {
                $job_state_ref->{"result"} = [];
            }
            return wantarray ? @{$job_state_ref->{"result"}} : $job_state_ref->{"result"}->[0];
        }
    }
}

sub _KButil_Localize_GenomeSet_submit {
    my($self, @args) = @_;
# Authentication: required
    if ((my $n = @args) != 1) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function _KButil_Localize_GenomeSet_submit (received $n, expecting 1)");
    }
    {
        my($params) = @args;
        my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
            my $msg = "Invalid arguments passed to _KButil_Localize_GenomeSet_submit:\n" . join("", map { "\t$_\n" } @_bad_arguments);
            Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
                                   method_name => '_KButil_Localize_GenomeSet_submit');
        }
    }
    my $context = undef;
    if ($self->{service_version}) {
        $context = {'service_ver' => $self->{service_version}};
    }
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_SetUtilities._KButil_Localize_GenomeSet_submit",
        params => \@args, context => $context});
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => '_KButil_Localize_GenomeSet_submit',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
            );
        } else {
            return $result->result->[0];  # job_id
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method _KButil_Localize_GenomeSet_submit",
                        status_line => $self->{client}->status_line,
                        method_name => '_KButil_Localize_GenomeSet_submit');
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
    my $job_id = $self->_KButil_Localize_FeatureSet_submit(@args);
    my $async_job_check_time = $self->{async_job_check_time};
    while (1) {
        Time::HiRes::sleep($async_job_check_time);
        $async_job_check_time *= $self->{async_job_check_time_scale_percent} / 100.0;
        if ($async_job_check_time > $self->{async_job_check_max_time}) {
            $async_job_check_time = $self->{async_job_check_max_time};
        }
        my $job_state_ref = $self->_check_job($job_id);
        if ($job_state_ref->{"finished"} != 0) {
            if (!exists $job_state_ref->{"result"}) {
                $job_state_ref->{"result"} = [];
            }
            return wantarray ? @{$job_state_ref->{"result"}} : $job_state_ref->{"result"}->[0];
        }
    }
}

sub _KButil_Localize_FeatureSet_submit {
    my($self, @args) = @_;
# Authentication: required
    if ((my $n = @args) != 1) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function _KButil_Localize_FeatureSet_submit (received $n, expecting 1)");
    }
    {
        my($params) = @args;
        my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
            my $msg = "Invalid arguments passed to _KButil_Localize_FeatureSet_submit:\n" . join("", map { "\t$_\n" } @_bad_arguments);
            Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
                                   method_name => '_KButil_Localize_FeatureSet_submit');
        }
    }
    my $context = undef;
    if ($self->{service_version}) {
        $context = {'service_ver' => $self->{service_version}};
    }
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_SetUtilities._KButil_Localize_FeatureSet_submit",
        params => \@args, context => $context});
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => '_KButil_Localize_FeatureSet_submit',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
            );
        } else {
            return $result->result->[0];  # job_id
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method _KButil_Localize_FeatureSet_submit",
                        status_line => $self->{client}->status_line,
                        method_name => '_KButil_Localize_FeatureSet_submit');
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
    my $job_id = $self->_KButil_Merge_FeatureSet_Collection_submit(@args);
    my $async_job_check_time = $self->{async_job_check_time};
    while (1) {
        Time::HiRes::sleep($async_job_check_time);
        $async_job_check_time *= $self->{async_job_check_time_scale_percent} / 100.0;
        if ($async_job_check_time > $self->{async_job_check_max_time}) {
            $async_job_check_time = $self->{async_job_check_max_time};
        }
        my $job_state_ref = $self->_check_job($job_id);
        if ($job_state_ref->{"finished"} != 0) {
            if (!exists $job_state_ref->{"result"}) {
                $job_state_ref->{"result"} = [];
            }
            return wantarray ? @{$job_state_ref->{"result"}} : $job_state_ref->{"result"}->[0];
        }
    }
}

sub _KButil_Merge_FeatureSet_Collection_submit {
    my($self, @args) = @_;
# Authentication: required
    if ((my $n = @args) != 1) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function _KButil_Merge_FeatureSet_Collection_submit (received $n, expecting 1)");
    }
    {
        my($params) = @args;
        my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
            my $msg = "Invalid arguments passed to _KButil_Merge_FeatureSet_Collection_submit:\n" . join("", map { "\t$_\n" } @_bad_arguments);
            Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
                                   method_name => '_KButil_Merge_FeatureSet_Collection_submit');
        }
    }
    my $context = undef;
    if ($self->{service_version}) {
        $context = {'service_ver' => $self->{service_version}};
    }
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_SetUtilities._KButil_Merge_FeatureSet_Collection_submit",
        params => \@args, context => $context});
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => '_KButil_Merge_FeatureSet_Collection_submit',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
            );
        } else {
            return $result->result->[0];  # job_id
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method _KButil_Merge_FeatureSet_Collection_submit",
                        status_line => $self->{client}->status_line,
                        method_name => '_KButil_Merge_FeatureSet_Collection_submit');
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
    my $job_id = $self->_KButil_Slice_FeatureSets_by_Genomes_submit(@args);
    my $async_job_check_time = $self->{async_job_check_time};
    while (1) {
        Time::HiRes::sleep($async_job_check_time);
        $async_job_check_time *= $self->{async_job_check_time_scale_percent} / 100.0;
        if ($async_job_check_time > $self->{async_job_check_max_time}) {
            $async_job_check_time = $self->{async_job_check_max_time};
        }
        my $job_state_ref = $self->_check_job($job_id);
        if ($job_state_ref->{"finished"} != 0) {
            if (!exists $job_state_ref->{"result"}) {
                $job_state_ref->{"result"} = [];
            }
            return wantarray ? @{$job_state_ref->{"result"}} : $job_state_ref->{"result"}->[0];
        }
    }
}

sub _KButil_Slice_FeatureSets_by_Genomes_submit {
    my($self, @args) = @_;
# Authentication: required
    if ((my $n = @args) != 1) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function _KButil_Slice_FeatureSets_by_Genomes_submit (received $n, expecting 1)");
    }
    {
        my($params) = @args;
        my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
            my $msg = "Invalid arguments passed to _KButil_Slice_FeatureSets_by_Genomes_submit:\n" . join("", map { "\t$_\n" } @_bad_arguments);
            Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
                                   method_name => '_KButil_Slice_FeatureSets_by_Genomes_submit');
        }
    }
    my $context = undef;
    if ($self->{service_version}) {
        $context = {'service_ver' => $self->{service_version}};
    }
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_SetUtilities._KButil_Slice_FeatureSets_by_Genomes_submit",
        params => \@args, context => $context});
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => '_KButil_Slice_FeatureSets_by_Genomes_submit',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
            );
        } else {
            return $result->result->[0];  # job_id
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method _KButil_Slice_FeatureSets_by_Genomes_submit",
                        status_line => $self->{client}->status_line,
                        method_name => '_KButil_Slice_FeatureSets_by_Genomes_submit');
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
    my $job_id = $self->_KButil_Logical_Slice_Two_FeatureSets_submit(@args);
    my $async_job_check_time = $self->{async_job_check_time};
    while (1) {
        Time::HiRes::sleep($async_job_check_time);
        $async_job_check_time *= $self->{async_job_check_time_scale_percent} / 100.0;
        if ($async_job_check_time > $self->{async_job_check_max_time}) {
            $async_job_check_time = $self->{async_job_check_max_time};
        }
        my $job_state_ref = $self->_check_job($job_id);
        if ($job_state_ref->{"finished"} != 0) {
            if (!exists $job_state_ref->{"result"}) {
                $job_state_ref->{"result"} = [];
            }
            return wantarray ? @{$job_state_ref->{"result"}} : $job_state_ref->{"result"}->[0];
        }
    }
}

sub _KButil_Logical_Slice_Two_FeatureSets_submit {
    my($self, @args) = @_;
# Authentication: required
    if ((my $n = @args) != 1) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function _KButil_Logical_Slice_Two_FeatureSets_submit (received $n, expecting 1)");
    }
    {
        my($params) = @args;
        my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
            my $msg = "Invalid arguments passed to _KButil_Logical_Slice_Two_FeatureSets_submit:\n" . join("", map { "\t$_\n" } @_bad_arguments);
            Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
                                   method_name => '_KButil_Logical_Slice_Two_FeatureSets_submit');
        }
    }
    my $context = undef;
    if ($self->{service_version}) {
        $context = {'service_ver' => $self->{service_version}};
    }
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_SetUtilities._KButil_Logical_Slice_Two_FeatureSets_submit",
        params => \@args, context => $context});
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => '_KButil_Logical_Slice_Two_FeatureSets_submit',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
            );
        } else {
            return $result->result->[0];  # job_id
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method _KButil_Logical_Slice_Two_FeatureSets_submit",
                        status_line => $self->{client}->status_line,
                        method_name => '_KButil_Logical_Slice_Two_FeatureSets_submit');
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
    my $job_id = $self->_KButil_Merge_GenomeSets_submit(@args);
    my $async_job_check_time = $self->{async_job_check_time};
    while (1) {
        Time::HiRes::sleep($async_job_check_time);
        $async_job_check_time *= $self->{async_job_check_time_scale_percent} / 100.0;
        if ($async_job_check_time > $self->{async_job_check_max_time}) {
            $async_job_check_time = $self->{async_job_check_max_time};
        }
        my $job_state_ref = $self->_check_job($job_id);
        if ($job_state_ref->{"finished"} != 0) {
            if (!exists $job_state_ref->{"result"}) {
                $job_state_ref->{"result"} = [];
            }
            return wantarray ? @{$job_state_ref->{"result"}} : $job_state_ref->{"result"}->[0];
        }
    }
}

sub _KButil_Merge_GenomeSets_submit {
    my($self, @args) = @_;
# Authentication: required
    if ((my $n = @args) != 1) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function _KButil_Merge_GenomeSets_submit (received $n, expecting 1)");
    }
    {
        my($params) = @args;
        my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
            my $msg = "Invalid arguments passed to _KButil_Merge_GenomeSets_submit:\n" . join("", map { "\t$_\n" } @_bad_arguments);
            Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
                                   method_name => '_KButil_Merge_GenomeSets_submit');
        }
    }
    my $context = undef;
    if ($self->{service_version}) {
        $context = {'service_ver' => $self->{service_version}};
    }
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_SetUtilities._KButil_Merge_GenomeSets_submit",
        params => \@args, context => $context});
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => '_KButil_Merge_GenomeSets_submit',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
            );
        } else {
            return $result->result->[0];  # job_id
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method _KButil_Merge_GenomeSets_submit",
                        status_line => $self->{client}->status_line,
                        method_name => '_KButil_Merge_GenomeSets_submit');
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
    my $job_id = $self->_KButil_Build_GenomeSet_submit(@args);
    my $async_job_check_time = $self->{async_job_check_time};
    while (1) {
        Time::HiRes::sleep($async_job_check_time);
        $async_job_check_time *= $self->{async_job_check_time_scale_percent} / 100.0;
        if ($async_job_check_time > $self->{async_job_check_max_time}) {
            $async_job_check_time = $self->{async_job_check_max_time};
        }
        my $job_state_ref = $self->_check_job($job_id);
        if ($job_state_ref->{"finished"} != 0) {
            if (!exists $job_state_ref->{"result"}) {
                $job_state_ref->{"result"} = [];
            }
            return wantarray ? @{$job_state_ref->{"result"}} : $job_state_ref->{"result"}->[0];
        }
    }
}

sub _KButil_Build_GenomeSet_submit {
    my($self, @args) = @_;
# Authentication: required
    if ((my $n = @args) != 1) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function _KButil_Build_GenomeSet_submit (received $n, expecting 1)");
    }
    {
        my($params) = @args;
        my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
            my $msg = "Invalid arguments passed to _KButil_Build_GenomeSet_submit:\n" . join("", map { "\t$_\n" } @_bad_arguments);
            Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
                                   method_name => '_KButil_Build_GenomeSet_submit');
        }
    }
    my $context = undef;
    if ($self->{service_version}) {
        $context = {'service_ver' => $self->{service_version}};
    }
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_SetUtilities._KButil_Build_GenomeSet_submit",
        params => \@args, context => $context});
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => '_KButil_Build_GenomeSet_submit',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
            );
        } else {
            return $result->result->[0];  # job_id
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method _KButil_Build_GenomeSet_submit",
                        status_line => $self->{client}->status_line,
                        method_name => '_KButil_Build_GenomeSet_submit');
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
    my $job_id = $self->_KButil_Build_GenomeSet_from_FeatureSet_submit(@args);
    my $async_job_check_time = $self->{async_job_check_time};
    while (1) {
        Time::HiRes::sleep($async_job_check_time);
        $async_job_check_time *= $self->{async_job_check_time_scale_percent} / 100.0;
        if ($async_job_check_time > $self->{async_job_check_max_time}) {
            $async_job_check_time = $self->{async_job_check_max_time};
        }
        my $job_state_ref = $self->_check_job($job_id);
        if ($job_state_ref->{"finished"} != 0) {
            if (!exists $job_state_ref->{"result"}) {
                $job_state_ref->{"result"} = [];
            }
            return wantarray ? @{$job_state_ref->{"result"}} : $job_state_ref->{"result"}->[0];
        }
    }
}

sub _KButil_Build_GenomeSet_from_FeatureSet_submit {
    my($self, @args) = @_;
# Authentication: required
    if ((my $n = @args) != 1) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function _KButil_Build_GenomeSet_from_FeatureSet_submit (received $n, expecting 1)");
    }
    {
        my($params) = @args;
        my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
            my $msg = "Invalid arguments passed to _KButil_Build_GenomeSet_from_FeatureSet_submit:\n" . join("", map { "\t$_\n" } @_bad_arguments);
            Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
                                   method_name => '_KButil_Build_GenomeSet_from_FeatureSet_submit');
        }
    }
    my $context = undef;
    if ($self->{service_version}) {
        $context = {'service_ver' => $self->{service_version}};
    }
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_SetUtilities._KButil_Build_GenomeSet_from_FeatureSet_submit",
        params => \@args, context => $context});
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => '_KButil_Build_GenomeSet_from_FeatureSet_submit',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
            );
        } else {
            return $result->result->[0];  # job_id
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method _KButil_Build_GenomeSet_from_FeatureSet_submit",
                        status_line => $self->{client}->status_line,
                        method_name => '_KButil_Build_GenomeSet_from_FeatureSet_submit');
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
	input_genome_refs has a value which is a kb_SetUtilities.data_obj_ref
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
	input_genome_refs has a value which is a kb_SetUtilities.data_obj_ref
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
    my $job_id = $self->_KButil_Add_Genomes_to_GenomeSet_submit(@args);
    my $async_job_check_time = $self->{async_job_check_time};
    while (1) {
        Time::HiRes::sleep($async_job_check_time);
        $async_job_check_time *= $self->{async_job_check_time_scale_percent} / 100.0;
        if ($async_job_check_time > $self->{async_job_check_max_time}) {
            $async_job_check_time = $self->{async_job_check_max_time};
        }
        my $job_state_ref = $self->_check_job($job_id);
        if ($job_state_ref->{"finished"} != 0) {
            if (!exists $job_state_ref->{"result"}) {
                $job_state_ref->{"result"} = [];
            }
            return wantarray ? @{$job_state_ref->{"result"}} : $job_state_ref->{"result"}->[0];
        }
    }
}

sub _KButil_Add_Genomes_to_GenomeSet_submit {
    my($self, @args) = @_;
# Authentication: required
    if ((my $n = @args) != 1) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function _KButil_Add_Genomes_to_GenomeSet_submit (received $n, expecting 1)");
    }
    {
        my($params) = @args;
        my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
            my $msg = "Invalid arguments passed to _KButil_Add_Genomes_to_GenomeSet_submit:\n" . join("", map { "\t$_\n" } @_bad_arguments);
            Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
                                   method_name => '_KButil_Add_Genomes_to_GenomeSet_submit');
        }
    }
    my $context = undef;
    if ($self->{service_version}) {
        $context = {'service_ver' => $self->{service_version}};
    }
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_SetUtilities._KButil_Add_Genomes_to_GenomeSet_submit",
        params => \@args, context => $context});
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => '_KButil_Add_Genomes_to_GenomeSet_submit',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
            );
        } else {
            return $result->result->[0];  # job_id
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method _KButil_Add_Genomes_to_GenomeSet_submit",
                        status_line => $self->{client}->status_line,
                        method_name => '_KButil_Add_Genomes_to_GenomeSet_submit');
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
    my $job_id = $self->_KButil_Build_ReadsSet_submit(@args);
    my $async_job_check_time = $self->{async_job_check_time};
    while (1) {
        Time::HiRes::sleep($async_job_check_time);
        $async_job_check_time *= $self->{async_job_check_time_scale_percent} / 100.0;
        if ($async_job_check_time > $self->{async_job_check_max_time}) {
            $async_job_check_time = $self->{async_job_check_max_time};
        }
        my $job_state_ref = $self->_check_job($job_id);
        if ($job_state_ref->{"finished"} != 0) {
            if (!exists $job_state_ref->{"result"}) {
                $job_state_ref->{"result"} = [];
            }
            return wantarray ? @{$job_state_ref->{"result"}} : $job_state_ref->{"result"}->[0];
        }
    }
}

sub _KButil_Build_ReadsSet_submit {
    my($self, @args) = @_;
# Authentication: required
    if ((my $n = @args) != 1) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function _KButil_Build_ReadsSet_submit (received $n, expecting 1)");
    }
    {
        my($params) = @args;
        my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
            my $msg = "Invalid arguments passed to _KButil_Build_ReadsSet_submit:\n" . join("", map { "\t$_\n" } @_bad_arguments);
            Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
                                   method_name => '_KButil_Build_ReadsSet_submit');
        }
    }
    my $context = undef;
    if ($self->{service_version}) {
        $context = {'service_ver' => $self->{service_version}};
    }
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_SetUtilities._KButil_Build_ReadsSet_submit",
        params => \@args, context => $context});
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => '_KButil_Build_ReadsSet_submit',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
            );
        } else {
            return $result->result->[0];  # job_id
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method _KButil_Build_ReadsSet_submit",
                        status_line => $self->{client}->status_line,
                        method_name => '_KButil_Build_ReadsSet_submit');
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
    my $job_id = $self->_KButil_Merge_MultipleReadsSets_to_OneReadsSet_submit(@args);
    my $async_job_check_time = $self->{async_job_check_time};
    while (1) {
        Time::HiRes::sleep($async_job_check_time);
        $async_job_check_time *= $self->{async_job_check_time_scale_percent} / 100.0;
        if ($async_job_check_time > $self->{async_job_check_max_time}) {
            $async_job_check_time = $self->{async_job_check_max_time};
        }
        my $job_state_ref = $self->_check_job($job_id);
        if ($job_state_ref->{"finished"} != 0) {
            if (!exists $job_state_ref->{"result"}) {
                $job_state_ref->{"result"} = [];
            }
            return wantarray ? @{$job_state_ref->{"result"}} : $job_state_ref->{"result"}->[0];
        }
    }
}

sub _KButil_Merge_MultipleReadsSets_to_OneReadsSet_submit {
    my($self, @args) = @_;
# Authentication: required
    if ((my $n = @args) != 1) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function _KButil_Merge_MultipleReadsSets_to_OneReadsSet_submit (received $n, expecting 1)");
    }
    {
        my($params) = @args;
        my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
            my $msg = "Invalid arguments passed to _KButil_Merge_MultipleReadsSets_to_OneReadsSet_submit:\n" . join("", map { "\t$_\n" } @_bad_arguments);
            Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
                                   method_name => '_KButil_Merge_MultipleReadsSets_to_OneReadsSet_submit');
        }
    }
    my $context = undef;
    if ($self->{service_version}) {
        $context = {'service_ver' => $self->{service_version}};
    }
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_SetUtilities._KButil_Merge_MultipleReadsSets_to_OneReadsSet_submit",
        params => \@args, context => $context});
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => '_KButil_Merge_MultipleReadsSets_to_OneReadsSet_submit',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
            );
        } else {
            return $result->result->[0];  # job_id
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method _KButil_Merge_MultipleReadsSets_to_OneReadsSet_submit",
                        status_line => $self->{client}->status_line,
                        method_name => '_KButil_Merge_MultipleReadsSets_to_OneReadsSet_submit');
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
    my $job_id = $self->_KButil_Build_AssemblySet_submit(@args);
    my $async_job_check_time = $self->{async_job_check_time};
    while (1) {
        Time::HiRes::sleep($async_job_check_time);
        $async_job_check_time *= $self->{async_job_check_time_scale_percent} / 100.0;
        if ($async_job_check_time > $self->{async_job_check_max_time}) {
            $async_job_check_time = $self->{async_job_check_max_time};
        }
        my $job_state_ref = $self->_check_job($job_id);
        if ($job_state_ref->{"finished"} != 0) {
            if (!exists $job_state_ref->{"result"}) {
                $job_state_ref->{"result"} = [];
            }
            return wantarray ? @{$job_state_ref->{"result"}} : $job_state_ref->{"result"}->[0];
        }
    }
}

sub _KButil_Build_AssemblySet_submit {
    my($self, @args) = @_;
# Authentication: required
    if ((my $n = @args) != 1) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function _KButil_Build_AssemblySet_submit (received $n, expecting 1)");
    }
    {
        my($params) = @args;
        my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
            my $msg = "Invalid arguments passed to _KButil_Build_AssemblySet_submit:\n" . join("", map { "\t$_\n" } @_bad_arguments);
            Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
                                   method_name => '_KButil_Build_AssemblySet_submit');
        }
    }
    my $context = undef;
    if ($self->{service_version}) {
        $context = {'service_ver' => $self->{service_version}};
    }
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_SetUtilities._KButil_Build_AssemblySet_submit",
        params => \@args, context => $context});
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => '_KButil_Build_AssemblySet_submit',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
            );
        } else {
            return $result->result->[0];  # job_id
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method _KButil_Build_AssemblySet_submit",
                        status_line => $self->{client}->status_line,
                        method_name => '_KButil_Build_AssemblySet_submit');
    }
}

 
 
sub status
{
    my($self, @args) = @_;
    my $job_id = undef;
    if ((my $n = @args) != 0) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function status (received $n, expecting 0)");
    }
    my $context = undef;
    if ($self->{service_version}) {
        $context = {'service_ver' => $self->{service_version}};
    }
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_SetUtilities._status_submit",
        params => \@args, context => $context});
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => '_status_submit',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
            );
        } else {
            $job_id = $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method _status_submit",
                        status_line => $self->{client}->status_line,
                        method_name => '_status_submit');
    }
    my $async_job_check_time = $self->{async_job_check_time};
    while (1) {
        Time::HiRes::sleep($async_job_check_time);
        $async_job_check_time *= $self->{async_job_check_time_scale_percent} / 100.0;
        if ($async_job_check_time > $self->{async_job_check_max_time}) {
            $async_job_check_time = $self->{async_job_check_max_time};
        }
        my $job_state_ref = $self->_check_job($job_id);
        if ($job_state_ref->{"finished"} != 0) {
            if (!exists $job_state_ref->{"result"}) {
                $job_state_ref->{"result"} = [];
            }
            return wantarray ? @{$job_state_ref->{"result"}} : $job_state_ref->{"result"}->[0];
        }
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
                method_name => 'KButil_Build_AssemblySet',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method KButil_Build_AssemblySet",
            status_line => $self->{client}->status_line,
            method_name => 'KButil_Build_AssemblySet',
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
input_genome_refs has a value which is a kb_SetUtilities.data_obj_ref
input_genomeset_ref has a value which is a kb_SetUtilities.data_obj_ref
output_name has a value which is a kb_SetUtilities.data_obj_name
desc has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace_name has a value which is a kb_SetUtilities.workspace_name
input_genome_refs has a value which is a kb_SetUtilities.data_obj_ref
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
