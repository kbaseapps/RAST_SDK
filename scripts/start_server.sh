#!/bin/bash
if [ -L $0 ] ; then
script_dir=$(cd "$(dirname "$(readlink $0)")"; pwd -P)
else
script_dir=$(cd "$(dirname "$0")"; pwd -P)
fi
export KB_DEPLOYMENT_CONFIG=$script_dir/../deploy.cfg
export PERL5LIB=$script_dir/../lib:$PERL5LIB
plackup $script_dir/../lib/RAST_SDK.psgi
