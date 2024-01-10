FROM condaforge/mambaforge
ENV VERSION='2.1.4'
ENV CORES='16'

RUN apt update && \
    ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y tzdata && \
    dpkg-reconfigure --frontend noninteractive tzdata && \
    apt-get clean
RUN mamba install -c bioconda -y snakemake==7.32.4 fastp==0.23.4 && \
    cd /opt && \
    wget https://github.com/Kennedy-Lab-UW/Duplex-Seq-Pipeline/archive/refs/tags/v$VERSION.tar.gz && \
    tar zxvf v$VERSION.tar.gz && rm v$VERSION.tar.gz && \
    ln -s Duplex-Seq-Pipeline-$VERSION Duplex-Seq-Pipeline && \
    cd Duplex-Seq-Pipeline && \
    bash setupDS.sh $CORES && \
    rm -rf test
RUN mamba install -y -c bioconda snpeff==5.2 snpsift==5.2
ENV PATH=/opt/Duplex-Seq-Pipeline:$PATH