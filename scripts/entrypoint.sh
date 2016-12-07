#!/bin/bash

. /kb/deployment/user-env.sh

python ./scripts/prepare_deploy_cfg.py ./deploy.cfg ./work/config.properties

export PATH=$PATH:$KB_TOP/services/genome_annotation/bin:$KB_TOP/services/cdmi_api/bin
export KB_DEPLOYMENT_CONFIG=/kb/module/deploy.cfg

( cd /kb/deployment/services/kmer_annotation_figfam/;./start_service & )
if [ $# -eq 0 ] ; then
  sh ./scripts/start_server.sh
elif [ "${1}" = "test" ] ; then
  echo "Run Tests"
  make test
elif [ "${1}" = "async" ] ; then
  sh ./scripts/run_async.sh
elif [ "${1}" = "init" ] ; then
  echo "Initialize module"
  cd /data
  curl http://bioseed.mcs.anl.gov/~chenry/kmer.tgz|tar xzf -
  ln -s /data/kmer/Release70 /data/kmer/ACTIVE/Release70
  ln -s /data/kmer/Release70 /data/kmer/DEFAULT
  if [ -d kmer ] ; then
  	cd ..
  	touch __READY__
  fi
elif [ "${1}" = "bash" ] ; then
  bash
elif [ "${1}" = "report" ] ; then
  export KB_SDK_COMPILE_REPORT_FILE=./work/compile_report.json
  make compile
else
  echo Unknown
fi
