FROM whitelab/plumed-gromacs

MAINTAINER Andrew White

USER root

#Install packmol
RUN apt-get install -y gfortran make && apt-get clean
ADD packmol.tar.gz /usr/share
RUN cd /usr/share/packmol && ./configure && make && ln -s /usr/share/packmol/packmol  /usr/bin/packmol

#Get tools for testing
RUN apt-get install -y python-tox python-setuptools python-pip python-numpy python-scipy && apt-get clean

#set default
ADD test.sh /home/whitelab
WORKDIR /home/$USER
USER whitelab
CMD ["bash", "/home/whitelab/test.sh"]