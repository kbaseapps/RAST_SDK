#!/bin/bash
if [ -L $0 ] ; then
script_dir=$(cd "$(dirname "$(readlink $0)")"; pwd -P)
else
script_dir=$(cd "$(dirname "$0")"; pwd -P)
fi
export PERL5LIB=$script_dir/../lib:$PERL5LIB
perl $script_dir/../lib/RAST_SDK/RAST_SDKServer.pm $1 $2 $3
