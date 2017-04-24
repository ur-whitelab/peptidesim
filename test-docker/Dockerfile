FROM ubuntu:xenial

MAINTAINER Andrew White

USER root

#set default shell to bash
RUN echo "dash dash/sh boolean false" | debconf-set-selections
RUN dpkg-reconfigure -f noninteractive dash

#install gromacs

#enable contributed packages
RUN apt-get update
RUN apt-get install -y gromacs-openmpi gromacs-data ssh fftw3
#RUN apt-get install -y -f -t wheezy gromacs=4.5.5-2 gromacs-data=4.5.5-2

#Install packmol
RUN apt-get install -y gfortran make
ADD packmol.tar.gz /usr/share
RUN cd /usr/share/packmol && ./configure && make && ln -s /usr/share/packmol/packmol  /sbin/packmol

#Get tools for testing
RUN apt-get install -y python-tox python-setuptools python-pip python-numpy python-scipy emacs

#create user
RUN useradd -ms /bin/bash tester

#set default
ADD test.sh /home/tester
CMD ["bash", "/home/tester/test.sh"]