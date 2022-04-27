FROM kbase/rast_base:1.9.1

WORKDIR /kb/module

COPY ./ /kb/module

    # try to force tmp files to be written to right spot
RUN rm -rf /tmp && ln -s /kb/module/work/tmp /tmp && \
    mkdir -p /data                 && \
    mkdir -p /data/Data.may1       && \
    mkdir -p /data/kmer            && \
    mkdir -p /kb/module/work    && \
    chmod -R 777 /kb/module     && \
    cd /kb/module               && \
    make

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]

