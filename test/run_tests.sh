#!/bin/bash
echo "Running $0 with args $@"
if [ -L $0 ] ; then
test_dir=$(cd "$(dirname "$(readlink $0)")"; pwd -P) # for symbolic link
else
test_dir=$(cd "$(dirname "$0")"; pwd -P) # for normal file
fi
base_dir=$(cd $test_dir && cd .. && pwd);
export KB_DEPLOYMENT_CONFIG=$base_dir/deploy.cfg
export KB_AUTH_TOKEN=`cat /kb/module/work/token`
export PERL5LIB=$base_dir/lib:$PERL5LIB
cd $base_dir
prove -I $test_dir -lvrm $test_dir
