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
  mkdir kmer
  cd kmer
  curl -X GET https://ci.kbase.us/services/shock-api/node/a4d6c083-af46-4f1d-aea8-e3c4d4a1598b?download --user kbasetest:@Suite525|tar xzf -
  if [ -d Data.2 ] ; then
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
