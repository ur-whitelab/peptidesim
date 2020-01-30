FROM whitelab/plumed-gromacs

MAINTAINER Andrew White

USER root

#Install packmol
RUN apt-get install -y gfortran make && apt-get clean
ADD packmol.tar.gz /usr/share
RUN cd /usr/share/packmol && ./configure && make && ln -s /usr/share/packmol/packmol  /usr/bin/packmol

#Get tools for testing
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 git mercurial subversion && \
    apt-get clean

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy

RUN pip install --no-cache-dir git+https://github.com/whitead/GromacsWrapper tox;

USER whitelab

RUN echo "import gromacs; gromacs.config.setup()" | python

#set default
ADD test.sh /home/whitelab
ADD interact.sh /home/whitelab
CMD ["bash", "/home/whitelab/test.sh"]