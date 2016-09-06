FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.

# RUN apt-get update
RUN cpanm -i Config::IniFiles

ADD ./bootstrap bootstrap

RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/dev_container/modules && \
  rm -rf jars && \
  git clone https://github.com/kbase/jars && \
  rm -rf kb_sdk && \
  git clone https://github.com/kbase/kb_sdk -b develop && \
  cd /kb/dev_container/modules/jars && \
  make deploy && \
  cd /kb/dev_container/modules/kb_sdk && \
  make && make deploy
RUN \
  . /kb/dev_container/user-env.sh && \
  cd /kb/dev_container/modules && \
  rm -rf data_api && \
  git clone https://github.com/kbase/data_api -b develop && \
  pip install --upgrade /kb/dev_container/modules/data_api

RUN \
    cd bootstrap && \
    ./build.search_for_rnas /kb/runtime/ && \
    ./build.glimmer /kb/runtime/ && \
    ./build.elph /kb/runtime/ && \
    ./build.prodigal /kb/runtime/ && \
    cd .. && rm -rf bootstrap

# Build kb_seed
RUN cd /kb/dev_container/modules && \
    rm -rf kb_seed strep_repeats kmer_annotation_figfam genome_annotation && \
    git clone https://github.com/kbase/kb_seed && \
    git clone https://github.com/kbase/strep_repeats && \
    git clone https://github.com/kbase/kmer_annotation_figfam && \
    git clone https://github.com/kbase/genome_annotation && \
    . /kb/dev_container/user-env.sh && \
    cd kb_seed && make && make TARGET=/kb/deployment deploy && cd .. && \
    cd strep_repeats && make && make TARGET=/kb/deployment deploy && cd ..&& \
    cd kmer_annotation_figfam && make && make TARGET=/kb/deployment deploy && cd ..&& \
    cd genome_annotation && make && make TARGET=/kb/deployment deploy && cd .. && \
    sed -i 's/print .*keeping.*/#ignore/'  /kb/deployment/lib/GenomeTypeObject.pm


#RUN sed -i 's/capture_stderr/tee_stderr/' /kb/deployment/lib/Bio/KBase/GenomeAnnotation/GenomeAnnotationImpl.pm

RUN \
    cd /kb/deployment/services/kmer_annotation_figfam/ && \
    sed 's|$KB_TOP/deployment.cfg|/kb/module/deploy.cfg|' -i ./start_service  && \
    sed 's/8/1/' -i ./start_service 

RUN mkdir /data && \
    mkdir /data/Data.may1 && \
    mkdir /data/kmer

#RUN sed -i 's/->port/->port, Passive=>1/' /kb/deployment/plbin/kmer-figfam-update-data.pl
                                                           
# -----------------------------------------

COPY ./ /kb/module

RUN mkdir -p /kb/module/work

WORKDIR /kb/module

RUN make

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
