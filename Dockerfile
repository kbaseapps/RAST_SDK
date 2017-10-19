FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.

# RUN apt-get update
RUN cpanm -i Config::IniFiles && \
    cpanm -i UUID::Random && \
    cpanm -i HTML::SimpleLinkExtor && \
    cpanm -i WWW::Mechanize --force && \
    cpanm -i MIME::Base64 && \
    apt-get -y install nano

ADD ./bootstrap bootstrap

RUN \
    cd bootstrap && \
    ./build.search_for_rnas /kb/runtime/ && \
    ./build.glimmer /kb/runtime/ && \
    ./build.elph /kb/runtime/ && \
    ./build.prodigal /kb/runtime/ && \
    ./build.phispy /kb/runtime/ && \
    cd .. && rm -rf bootstrap

# Add random forest for phispy
RUN \
    wget https://cran.r-project.org/src/contrib/randomForest_4.6-12.tar.gz && \
    R CMD INSTALL ./randomForest_4.6-12.tar.gz && \
    rm randomForest_4.6-12.tar.gz


# Build kb_seed
RUN cd /kb/dev_container/modules && \
    rm -rf idserver kb_seed strep_repeats kmer_annotation_figfam genome_annotation && \
    git clone https://github.com/kbase/kb_seed && \
    git clone https://github.com/kbase/strep_repeats && \
    git clone https://github.com/kbase/kmer_annotation_figfam && \
    git clone https://github.com/kbase/genome_annotation && \
    git clone https://github.com/kbase/idserver && \
    . /kb/dev_container/user-env.sh && \
    cd kb_seed && make && make TARGET=/kb/deployment deploy && cd .. && \
    cd strep_repeats && make && make TARGET=/kb/deployment deploy && cd ..&& \
    cd kmer_annotation_figfam && make && make TARGET=/kb/deployment deploy && cd ..&& \
    cd genome_annotation && make && make TARGET=/kb/deployment deploy && cd .. && \
    cd idserver && make && make TARGET=/kb/deployment deploy && cd .. && \
    sed -i 's/print .*keeping.*/#ignore/'  /kb/deployment/lib/GenomeTypeObject.pm


#RUN sed -i 's/capture_stderr/tee_stderr/' /kb/deployment/lib/Bio/KBase/GenomeAnnotation/GenomeAnnotationImpl.pm

RUN \
    cpanm install Set::IntervalTree && \
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
RUN chmod -R 777 /kb/module
RUN chmod -R 777 /kb/deployment/services/kmer_annotation_figfam
WORKDIR /kb/module

RUN make

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
