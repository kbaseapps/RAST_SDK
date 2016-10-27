package KBaseReport::KBaseReportImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org
our $VERSION = "0.0.1";
our $GIT_URL = "https://github.com/kbaseapps/KBaseReport.git";
our $GIT_COMMIT_HASH = "86fbc21b0b98d066579277ea1277396b09d00644";

=head1 NAME

KBaseReport

=head1 DESCRIPTION

Module for a simple WS data object report type.

=cut

#BEGIN_HEADER
use Bio::KBase::AuthToken;
use Bio::KBase::workspace::Client;
use Bio::KBase::HandleService;
use Config::IniFiles;
use JSON;
use MIME::Base64;
use URI;
use HTML::Parser;
use WWW::Mechanize;
use HTML::LinkExtor;
use File::Slurp;
use HTML::TreeBuilder;
use HTML::SimpleLinkExtor;
use Archive::Zip qw( :ERROR_CODES :CONSTANTS );
use List::Util qw(first);
use IO::Compress::Zip qw(zip $ZipError) ;
use File::Copy qw(copy);
use File::Path;
use File::Basename;
use UUID::Random;
use Data::Dumper;

sub format_html_string_base64{
    my ($html_string) = @_;
    my @img_links = $html_string =~ /<img\b(?=\s)(?=(?:[^>=]|='[^']*'|="[^"]*"|=[^'"][^\s>]*)*?\ssrc=['"]([^"]*)['"])(?:[^>=]|='[^']*'|="[^"]*"|=[^'"\s]*)*\s?>/g;
    my %imgLinksHash;
    #print &Dumper (\@img_links);
    for (my $i=0; $i<@img_links; $i++){
        if ($img_links[$i] =~ /base64/){
            next;
        }
        else{
            my @image_type = split /\./, $img_links[$i];
            # prefix get added to the image url data:image/png;base64,
            my $leaderString = "data:image/$image_type[-1];base64,";
            my $img_link = "/kb/module/data/".$img_links[$i];
            my $encoded = MIME::Base64::encode_base64(read_file($img_link));
            chomp $encoded;
            my $dest = $leaderString.$encoded;
            #print "$img_links[$i]\t$dest\n";
            $html_string =~ s/$img_links[$i]/$dest/;
        }
    }

   return $html_string;
}

sub format_images_base64{

    my ($html_file_link) = @_;

    my @f_split = split /\//, $html_file_link;
    my $tmpDir = "/kb/module/work/tmp/tmphtml";
    mkpath([$tmpDir], 1);
    my $out_path = $tmpDir."/embedded_".$f_split[-1];

    open INFILE, "$html_file_link" or die "Couldn't open html file $!\n";
    open OUTFILE, ">$out_path" or die "Couldn't open html out file $!\n";
    my $extor = HTML::SimpleLinkExtor->new();
    $extor->parse_file($html_file_link);
    my @img_links = $extor->img;
    my %imgLinksHash;
    for (my $i=0; $i<@img_links; $i++){
        my @image_type = split /\./, $img_links[$i];
        # prefix get added to the image url data:image/png;base64,
        my $leaderString = "\"data:image/$image_type[-1];base64,";
        if ($img_links[$i] =~ /base64/){
            next;
        }
        else{
            my $encoded = MIME::Base64::encode_base64($img_links[$i]);
            chomp $encoded;
            $imgLinksHash{$img_links[$i]} = $leaderString.$encoded."\"";
        }
    }
    #print &Dumper (\%imgLinksHash);
    while (defined(my $input = <INFILE>)){
        chomp $input;
        if ($input =~ /\<img(.+)/ && $input !~ /base64/ ){
            if($input=~ /src=\"((\w|_|\\|-|\/|\.|:)+)\"/){
                my $old_src = $1;
                if (exists $imgLinksHash{$old_src} ){
                    #print "$imgLinksHash{$old_src}\n";
                    $input =~ s/(<\s*img\s+.*src\s*=\s*)(")?.*?(?(2)")([\s>])/$1$imgLinksHash{$old_src}$3/sig;
                    print OUTFILE "$input\n";
                }
            }
        }
        else{
            print OUTFILE "$input\n"
        }
    }
    close INFILE;
    close OUTFILE;
    return $out_path;
}

sub zip_archive {
    my ($path) = @_;
    my $zip = Archive::Zip->new();
    my $outpath;
    if (-d $path){
        print "processing html folder at $path\n";
        $zip->addTree( $path);
        my $tmpDir = "/kb/module/work/tmp/zippedHtml";
        mkpath([$tmpDir], 1);
        my @folder_name = split /\//, $path;
        $outpath = $tmpDir."/".$folder_name[-1].".zip";
        $zip->writeToFileNamed($outpath);
    }
    elsif (-f $path){
        print "processing html file at $path\n";
        my $tmpDir = "/kb/module/work/tmp/ZippedHtml";
        mkpath([$tmpDir], 1);
        copy $path, $tmpDir;
        $zip->addTree($tmpDir);
        my @folder_name = split /\//, $path;
        $outpath = $tmpDir."/".$folder_name[-1].".zip";
        $zip->writeToFileNamed($outpath);
    }
    else{
        print "Html input file/dir type unknown\n";
        die;
    }
    return $outpath;
}

sub curl_upload_shock {
    my ($file, $shock) = @_;
    my $token = $shock->{token};
    my $url   = $shock->{url};
    my $attr  = q('{"file":"reporter"}');
    my $cmd   = 'curl --connect-timeout 100 -s -X POST -F attributes=@- -F upload=@'.$file." $url/node ";
    $cmd     .= " -H 'Authorization: OAuth $token'";
    my $out   = `echo $attr | $cmd` or die "Connection timeout uploading file to Shock: $file\n";
    my $json  = decode_json($out);
    $json->{status} == 200 or die "Error uploading file: $file\n".$json->{status}." ".$json->{error}->[0]."\n";
    return $json;
}

sub create_handle{
    my ($shock, $shock_out, $handle_service) = @_;
    my $handle;
    $handle->{type} = 'shock';
    $handle->{url}  = $shock->{url};
    $handle->{id}   = $shock_out->{data}->{id};
    my $hid = $handle_service->persist_handle($handle);
    $handle->{hid} = $hid;
    return $handle;
}
sub generate_shock_url{
    my ($handle) = @_;
    return my $url = $handle->{url}."/node/".$handle->{id};
}
#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR

    my $config_file = $ENV{ KB_DEPLOYMENT_CONFIG };
    my $cfg = Config::IniFiles->new(-file=>$config_file);
    my $wsInstance = $cfg->val('KBaseReport','workspace-url');
    die "no workspace-url defined" unless $wsInstance;
    my $ShockInstance = $cfg->val('KBaseReport','shock-url');
    die "no shock-url defined" unless $ShockInstance;
    my $HandleInstance = $cfg->val('KBaseReport','handle-service-url');
    die "no handle-service-url defined" unless $HandleInstance;

    $self->{'workspace-url'} = $wsInstance;
    $self->{'shock-url'} = $ShockInstance;
    $self->{'handle-service-url'} = $HandleInstance;

    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}

=head1 METHODS



=head2 create

  $info = $obj->create($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a KBaseReport.CreateParams
$info is a KBaseReport.ReportInfo
CreateParams is a reference to a hash where the following keys are defined:
	report has a value which is a KBaseReport.Report
	workspace_name has a value which is a string
Report is a reference to a hash where the following keys are defined:
	text_message has a value which is a string
	warnings has a value which is a reference to a list where each element is a string
	objects_created has a value which is a reference to a list where each element is a KBaseReport.WorkspaceObject
	file_links has a value which is a reference to a list where each element is a KBaseReport.LinkedFile
	html_links has a value which is a reference to a list where each element is a KBaseReport.LinkedFile
	direct_html has a value which is a string
	direct_html_link_index has a value which is an int
WorkspaceObject is a reference to a hash where the following keys are defined:
	ref has a value which is a KBaseReport.ws_id
	description has a value which is a string
ws_id is a string
LinkedFile is a reference to a hash where the following keys are defined:
	handle has a value which is a KBaseReport.handle_ref
	description has a value which is a string
	name has a value which is a string
	URL has a value which is a string
handle_ref is a string
ReportInfo is a reference to a hash where the following keys are defined:
	ref has a value which is a KBaseReport.ws_id
	name has a value which is a string

</pre>

=end html

=begin text

$params is a KBaseReport.CreateParams
$info is a KBaseReport.ReportInfo
CreateParams is a reference to a hash where the following keys are defined:
	report has a value which is a KBaseReport.Report
	workspace_name has a value which is a string
Report is a reference to a hash where the following keys are defined:
	text_message has a value which is a string
	warnings has a value which is a reference to a list where each element is a string
	objects_created has a value which is a reference to a list where each element is a KBaseReport.WorkspaceObject
	file_links has a value which is a reference to a list where each element is a KBaseReport.LinkedFile
	html_links has a value which is a reference to a list where each element is a KBaseReport.LinkedFile
	direct_html has a value which is a string
	direct_html_link_index has a value which is an int
WorkspaceObject is a reference to a hash where the following keys are defined:
	ref has a value which is a KBaseReport.ws_id
	description has a value which is a string
ws_id is a string
LinkedFile is a reference to a hash where the following keys are defined:
	handle has a value which is a KBaseReport.handle_ref
	description has a value which is a string
	name has a value which is a string
	URL has a value which is a string
handle_ref is a string
ReportInfo is a reference to a hash where the following keys are defined:
	ref has a value which is a KBaseReport.ws_id
	name has a value which is a string


=end text



=item Description

Create a KBaseReport with a brief summary of an App run.

=back

=cut

sub create
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to create:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'create');
    }

    my $ctx = $KBaseReport::KBaseReportServer::CallContext;
    my($info);
    #BEGIN create
    my $token=$ctx->token;
    my $provenance=$ctx->provenance;
    my $wsClient=Bio::KBase::workspace::Client->new($self->{'workspace-url'},token=>$token);

    print "Creating KBase Report\n";

    if (!defined $params->{report}){

        print "Field report must be defined to save a report";
    }
    if (!defined $params->{workspace_name}){

        print "Field workspace_name must be defined to save a report";
    }

    my $uid = UUID::Random::generate;
    my $reportName = "report_".$uid;

    if (defined $params->{prefix}){
        $reportName=$params->{prefix}.$reportName;
    }
    print "Report Name $reportName\n";

    my $obj_info_list = undef;
    eval {
        $obj_info_list = $wsClient->save_objects({
            'workspace'=>$params->{workspace_name},
            'objects'=>[{
                'type'=>'KBaseReport.Report',
                'data'=>$params->{report},
                'name'=>$reportName,
                'meta'=> {},
                'hidden' => 1,
                'provenance'=>$provenance
            }]
        });
    };
    if ($@) {
        die "Error saving modified genome object to workspace:\n".$@;
    }
    my $wsRef = $obj_info_list->[0]->[6]."/".$obj_info_list->[0]->[0]."/".$obj_info_list->[0]->[4];
    $info = {
        ref => $wsRef,
        name => $obj_info_list->[0]->[1]
    };

    return $info;
    #END create
    my @_bad_returns;
    (ref($info) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"info\" (value was \"$info\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to create:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'create');
    }
    return($info);
}




=head2 create_extended_report

  $info = $obj->create_extended_report($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a KBaseReport.CreateExtendedReportParams
$info is a KBaseReport.ReportInfo
CreateExtendedReportParams is a reference to a hash where the following keys are defined:
	message has a value which is a string
	objects_created has a value which is a reference to a list where each element is a KBaseReport.WorkspaceObject
	warnings has a value which is a reference to a list where each element is a string
	html_links has a value which is a reference to a list where each element is a KBaseReport.File
	direct_html has a value which is a string
	direct_html_link_index has a value which is an int
	file_links has a value which is a reference to a list where each element is a KBaseReport.File
	report_object_name has a value which is a string
	workspace_name has a value which is a string
WorkspaceObject is a reference to a hash where the following keys are defined:
	ref has a value which is a KBaseReport.ws_id
	description has a value which is a string
ws_id is a string
File is a reference to a hash where the following keys are defined:
	path has a value which is a string
	shock_id has a value which is a string
	name has a value which is a string
	description has a value which is a string
ReportInfo is a reference to a hash where the following keys are defined:
	ref has a value which is a KBaseReport.ws_id
	name has a value which is a string

</pre>

=end html

=begin text

$params is a KBaseReport.CreateExtendedReportParams
$info is a KBaseReport.ReportInfo
CreateExtendedReportParams is a reference to a hash where the following keys are defined:
	message has a value which is a string
	objects_created has a value which is a reference to a list where each element is a KBaseReport.WorkspaceObject
	warnings has a value which is a reference to a list where each element is a string
	html_links has a value which is a reference to a list where each element is a KBaseReport.File
	direct_html has a value which is a string
	direct_html_link_index has a value which is an int
	file_links has a value which is a reference to a list where each element is a KBaseReport.File
	report_object_name has a value which is a string
	workspace_name has a value which is a string
WorkspaceObject is a reference to a hash where the following keys are defined:
	ref has a value which is a KBaseReport.ws_id
	description has a value which is a string
ws_id is a string
File is a reference to a hash where the following keys are defined:
	path has a value which is a string
	shock_id has a value which is a string
	name has a value which is a string
	description has a value which is a string
ReportInfo is a reference to a hash where the following keys are defined:
	ref has a value which is a KBaseReport.ws_id
	name has a value which is a string


=end text



=item Description

A more complex function to create a report that enables the user to specify files and html view that the report should link to

=back

=cut

sub create_extended_report
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to create_extended_report:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'create_extended_report');
    }

    my $ctx = $KBaseReport::KBaseReportServer::CallContext;
    my($info);
    #BEGIN create_extended_report
    my $workspace_name = $params->{workspace_name};
    my $token=$ctx->token;
    my $provenance=$ctx->provenance;
    my $wsClient=Bio::KBase::workspace::Client->new($self->{'workspace-url'},token=>$token);
    my $file_link_arr = $params->{file_links};
    my $html_link_arr = $params->{html_links};

	my $shock = { url => $self->{'shock-url'}, token => $token };
    my $handle_service = Bio::KBase::HandleService->new($self->{'handle-service-url'});
    my @file_arr;
    my @html_arr;
    my $html_string = format_html_string_base64($params->{direct_html});
    #print "$html_string\n";
    #die;

    for (my $i=0; $i< @$file_link_arr; $i++){
        my $shock_out;
        my $handle_return;
        my $LinkedFile = {};
        if (defined $file_link_arr->[$i]->{shock_id}){
            if ($handle_service) {
                $shock_out->{data}->{id} = $file_link_arr->[$i]->{shock_id};
                $handle_return = create_handle($shock, $shock_out, $handle_service);
            }
            my $url = generate_shock_url($handle_return);
            $LinkedFile = {
                handle => $handle_return->{hid},
                description => $file_link_arr->[$i]->{description},
                name => $file_link_arr->[$i]->{name},
                URL => $url
            };
        }
        else{
            if ((-f $file_link_arr->[$i]->{path}) && (defined $file_link_arr->[$i]->{path})){

                $shock_out = curl_upload_shock ($file_link_arr->[$i]->{path}, $shock);
                if ($handle_service) {
                    $handle_return = create_handle($shock, $shock_out, $handle_service);
                }
                my $url = generate_shock_url($handle_return);
                $LinkedFile = {
                    handle => $handle_return->{hid},
                    description => $file_link_arr->[$i]->{description},
                    name => $file_link_arr->[$i]->{name},
                    URL => $url
                };
            }
            else{
                print "Error!! cannot access file at $file_link_arr->[$i]->{path}\n";
                next;
            }
        }
        push (@file_arr, $LinkedFile);
    }

    for (my $i=0; $i< @$html_link_arr; $i++){
        my $shock_out;
        my $handle_return;
        my $LinkedFile = {};
        if (defined $html_link_arr->[$i]->{shock_id}){
            if ($handle_service) {
                $shock_out->{data}->{id} = $html_link_arr->[$i]->{shock_id};
                $handle_return = create_handle($shock, $shock_out, $handle_service);
            }
            my $url = generate_shock_url($handle_return);
            $LinkedFile = {
                handle => $handle_return->{hid},
                description => $html_link_arr->[$i]->{description},
                name => $html_link_arr->[$i]->{name},
                URL => $url
            };
        }

        else{

            if (defined $html_link_arr->[$i]->{path}){
                my $out_link = zip_archive ($html_link_arr->[$i]->{path});
                if ( (-f $out_link) && ($handle_service) )  {
                    $shock_out = curl_upload_shock ($out_link, $shock);
                    $handle_return = create_handle($shock, $shock_out, $handle_service);
                }
                my $url = generate_shock_url($handle_return);
                $LinkedFile = {
                    handle => $handle_return->{hid},
                    description => $html_link_arr->[$i]->{description},
                    name =>  $html_link_arr->[$i]->{name},
                    URL => $url
                };
                print "$url\n";
            }
            else{
                print "Error!! cannot access html file at $html_link_arr->[$i]->{path}\n";
                next;
            }
        }
        push (@html_arr, $LinkedFile);
    }

    my $report = {
        text_message => $params->{message},
        file_links => \@file_arr,
        html_links => \@html_arr,
        direct_html => $html_string,
        objects_created => $params->{objects_created}
    };

	my $obj_info_list = undef;
    eval {
        $obj_info_list = $wsClient->save_objects({
            'workspace'=>$workspace_name,
            'objects'=>[{
                'type'=>'KBaseReport.Report',
                'data'=>$report,
                'name'=>$params->{report_object_name},
                'hidden' => 0,
                'provenance'=>$provenance
            }]
        });
    };
    if ($@) {
        die "Error saving modified genome object to workspace:\n".$@;
    }
    print &Dumper ($obj_info_list);
    my $wsRef = $obj_info_list->[0]->[6]."/".$obj_info_list->[0]->[0]."/".$obj_info_list->[0]->[4];
    $info = {
    	ref => $wsRef,
    	name => $params->{report_object_name}

    };

    return $info;
    #END create_extended_report
    my @_bad_returns;
    (ref($info) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"info\" (value was \"$info\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to create_extended_report:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'create_extended_report');
    }
    return($info);
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



=head2 ws_id

=over 4



=item Description

@id ws


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



=head2 handle_ref

=over 4



=item Description

Reference to a handle
@id handle


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



=head2 WorkspaceObject

=over 4



=item Description

Represents a Workspace object with some brief description text
that can be associated with the object.
@optional description


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
ref has a value which is a KBaseReport.ws_id
description has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
ref has a value which is a KBaseReport.ws_id
description has a value which is a string


=end text

=back



=head2 LinkedFile

=over 4



=item Description

Represents a file or html archive that the report should like to
@optional description


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
handle has a value which is a KBaseReport.handle_ref
description has a value which is a string
name has a value which is a string
URL has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
handle has a value which is a KBaseReport.handle_ref
description has a value which is a string
name has a value which is a string
URL has a value which is a string


=end text

=back



=head2 Report

=over 4



=item Description

A simple Report of a method run in KBase.
It only provides for now a way to display a fixed width text output summary message, a
list of warnings, and a list of objects created (each with descriptions).
@optional warnings file_links html_links direct_html direct_html_link_index
@metadata ws length(warnings) as Warnings
@metadata ws length(text_message) as Size(characters)
@metadata ws length(objects_created) as Objects Created


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
text_message has a value which is a string
warnings has a value which is a reference to a list where each element is a string
objects_created has a value which is a reference to a list where each element is a KBaseReport.WorkspaceObject
file_links has a value which is a reference to a list where each element is a KBaseReport.LinkedFile
html_links has a value which is a reference to a list where each element is a KBaseReport.LinkedFile
direct_html has a value which is a string
direct_html_link_index has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
text_message has a value which is a string
warnings has a value which is a reference to a list where each element is a string
objects_created has a value which is a reference to a list where each element is a KBaseReport.WorkspaceObject
file_links has a value which is a reference to a list where each element is a KBaseReport.LinkedFile
html_links has a value which is a reference to a list where each element is a KBaseReport.LinkedFile
direct_html has a value which is a string
direct_html_link_index has a value which is an int


=end text

=back



=head2 CreateParams

=over 4



=item Description

Provide the report information.  The structure is:
    params = {
        report: {
            text_message: '',
            warnings: ['w1'],
            objects_created: [ {
                ref: 'ws/objid',
                description: ''
            }]
        },
        workspace_name: 'ws'
    }


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report has a value which is a KBaseReport.Report
workspace_name has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report has a value which is a KBaseReport.Report
workspace_name has a value which is a string


=end text

=back



=head2 ReportInfo

=over 4



=item Description

The reference to the saved KBaseReport.  The structure is:
    reportInfo = {
        ref: 'ws/objid/ver',
        name: 'myreport.2262323452'
    }


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
ref has a value which is a KBaseReport.ws_id
name has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
ref has a value which is a KBaseReport.ws_id
name has a value which is a string


=end text

=back



=head2 File

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
path has a value which is a string
shock_id has a value which is a string
name has a value which is a string
description has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
path has a value which is a string
shock_id has a value which is a string
name has a value which is a string
description has a value which is a string


=end text

=back



=head2 CreateExtendedReportParams

=over 4



=item Description

Parameters used to create a more complex report with file and html links
The following arguments allow the user to specify the classical data fields in the report object:
string message - simple text message to store in report object
list <WorkspaceObject> objects_created;
list <string> warnings - a list of warning messages in simple text
The following argument allows the user to specify the location of html files/directories that the report widget will render <or> link to:
list <fileRef> html_links - a list of paths or shock node IDs pointing to a single flat html file or to the top level directory of a website
The report widget can render one html view directly. Set one of the following fields to decide which view to render:
string direct_html - simple html text that will be rendered within the report widget
int  direct_html_link_index - use this to specify the index of the page in html_links to view directly in the report widget (ignored if html_string is set)
The following argument allows the user to specify the location of files that the report widget should link for download:
list <fileRef> file_links - a list of paths or shock node IDs pointing to a single flat file
The following parameters indicate where the report object should be saved in the workspace:
string report_object_name - name to use for the report object (job ID is used if left unspecified)
string workspace_name - name of workspace where object should be saved


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
message has a value which is a string
objects_created has a value which is a reference to a list where each element is a KBaseReport.WorkspaceObject
warnings has a value which is a reference to a list where each element is a string
html_links has a value which is a reference to a list where each element is a KBaseReport.File
direct_html has a value which is a string
direct_html_link_index has a value which is an int
file_links has a value which is a reference to a list where each element is a KBaseReport.File
report_object_name has a value which is a string
workspace_name has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
message has a value which is a string
objects_created has a value which is a reference to a list where each element is a KBaseReport.WorkspaceObject
warnings has a value which is a reference to a list where each element is a string
html_links has a value which is a reference to a list where each element is a KBaseReport.File
direct_html has a value which is a string
direct_html_link_index has a value which is an int
file_links has a value which is a reference to a list where each element is a KBaseReport.File
report_object_name has a value which is a string
workspace_name has a value which is a string


=end text

=back



=cut

1;
