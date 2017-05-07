FROM whitelab/plumed-gromacs

MAINTAINER Andrew White

USER root

#Install packmol
RUN apt-get install -y gfortran make && apt-get clean
ADD packmol.tar.gz /usr/share
RUN cd /usr/share/packmol && ./configure && make && ln -s /usr/share/packmol/packmol  /usr/bin/packmol

#Get tools for testing
RUN apt-get install -y python-tox python-setuptools python-pip python-numpy python-scipy && apt-get clean
RUN apt-get install -y emacs

ENV USER=whitelab
USER whitelab
WORKDIR /home/whitelab

#setup gromacswrapper
RUN pip install --no-cache-dir --user GromacsWrapper && echo "import gromacs; gromacs.config.setup()" | python
RUN sed -i 's/tools = /tools = gmx gmx_mpi/' /home/whitelab/.gromacswrapper.cfg


    

#set default
ADD test.sh /home/whitelab
CMD ["bash", "/home/whitelab/test.sh"]