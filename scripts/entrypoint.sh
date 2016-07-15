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
  /kb/deployment/services/kmer_annotation_figfam/bin/kmer-figfam-update-data Release70
  cd ..
  curl ftp://ftp.theseed.org/KmerClassification/Data.may1.tgz|tar xzf -
  touch __READY__
elif [ "${1}" = "bash" ] ; then
  bash
elif [ "${1}" = "report" ] ; then
  export KB_SDK_COMPILE_REPORT_FILE=./work/compile_report.json
  make compile
else
  echo Unknown
fi
