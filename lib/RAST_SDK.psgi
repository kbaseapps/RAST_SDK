use RAST_SDK::RAST_SDKImpl;

use RAST_SDK::RAST_SDKServer;
use Plack::Middleware::CrossOrigin;



my @dispatch;

{
    my $obj = RAST_SDK::RAST_SDKImpl->new;
    push(@dispatch, 'RAST_SDK' => $obj);
}


my $server = RAST_SDK::RAST_SDKServer->new(instance_dispatch => { @dispatch },
				allow_get => 0,
			       );

my $handler = sub { $server->handle_input(@_) };

$handler = Plack::Middleware::CrossOrigin->wrap( $handler, origins => "*", headers => "*");
