FROM kbase/kbase:sdkbase2.latest

ADD ./bootstrap bootstrap

RUN \
    apt-get -y update && \
    apt-get -y install nano bioperl uuid-runtime && \
    cd bootstrap && \
    ./build.search_for_rnas /kb/runtime/ && \
    ./build.glimmer /kb/runtime/ && \
    ./build.elph /kb/runtime/ && \
    ./build.prodigal /kb/runtime/ && \
    ./build.phispy /kb/runtime/ && \
    cd .. && rm -rf bootstrap && \
    # Add random forest for phispy \
    wget https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-12.tar.gz && \
    R CMD INSTALL ./randomForest_4.6-12.tar.gz && \
    rm randomForest_4.6-12.tar.gz && \
    # Build kb_seed \
    cd /kb/dev_container/modules && \
    rm -rf idserver kb_seed strep_repeats kmer_annotation_figfam genome_annotation && \
    git clone https://github.com/kbase/kb_seed && \
    git clone https://github.com/kbase/strep_repeats && \
    git clone https://github.com/kbase/kmer_annotation_figfam && \
    git clone https://github.com/kbase/genome_annotation && \
    git clone https://github.com/kbase/idserver && \
    . /kb/dev_container/user-env.sh && \
    cd kb_seed && git checkout 20190314 && make && make TARGET=/kb/deployment deploy && cd .. && \
    cd strep_repeats && make && make TARGET=/kb/deployment deploy && cd .. && \
    cd kmer_annotation_figfam && make && make TARGET=/kb/deployment deploy && cd .. && \
    cd genome_annotation && make && make TARGET=/kb/deployment deploy && cd .. && \
    cd idserver && make && make TARGET=/kb/deployment deploy && cd .. && \
    # local file edits \
    sed -i 's/print .*keeping.*/#ignore/'  /kb/deployment/lib/GenomeTypeObject.pm && \
    cd /kb/deployment/services/kmer_annotation_figfam/ && \
    chmod -R 777 /kb/deployment/services/kmer_annotation_figfam && \
    sed 's|$KB_TOP/deployment.cfg|/kb/module/deploy.cfg|' -i ./start_service  && \
    sed 's|$KB_TOP/services/kmer_annotation_figfam|/tmp/|' -i ./start_service  && \
    sed 's/8/1/' -i ./start_service

# install perl deps
COPY ./cpanfile /kb/module/cpanfile
WORKDIR /kb/module
RUN cpanm --installdeps .

COPY ./ /kb/module

    # try to force tmp files to be written to right spot
RUN rm -rf /tmp && ln -s /kb/module/work/tmp /tmp && \
    mkdir /data                 && \
    mkdir /data/Data.may1       && \
    mkdir /data/kmer            && \
    mkdir -p /kb/module/work    && \
    chmod -R 777 /kb/module     && \
    cd /kb/module               && \
    make

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
