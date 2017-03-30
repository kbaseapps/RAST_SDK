#!/bin/bash
script_dir=$(dirname "$(readlink -f "$0")")
export KB_DEPLOYMENT_CONFIG=$script_dir/../deploy.cfg
export KB_AUTH_TOKEN=`cat /kb/module/work/token`
export PERL5LIB=$script_dir/../lib:$PATH:$PERL5LIB
#Used to mock out certain binaries (e.g. kmer_guts)
#export PATH=$script_dir/../test/:$PATH
cd $script_dir/../test
perl -e 'opendir my $dh, "."; my @l = grep { /\\.pl$/ } readdir $dh; foreach my $s (@l) { print("Running ".$s."\n"); system "perl", $s; }'
