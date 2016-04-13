FROM jupyter/scipy-notebook

MAINTAINER Andrew White

USER root

#make it easy to execuate as root
RUN echo '%sudo ALL=(ALL) NOPASSWD: ALL' > /etc/sudoers.d/00nopass && chmod 0440 /etc/sudoers.d/00nopass
RUN adduser jovyan sudo
RUN passwd jovyan -d

#set default shell to bash
RUN echo "dash dash/sh boolean false" | debconf-set-selections
RUN dpkg-reconfigure -f noninteractive dash

#install gromacs

#enable contributed packages
RUN sed -i 's/main/main contrib/g' /etc/apt/sources.list
RUN apt-get update
RUN apt-get install -y gromacs gromacs-data

#configure the GromacsWrapper
RUN source activate python2 && pip install git+git://github.com/whitead/GromacsWrapper.git@develop && echo "import gromacs;  gromacs.config.setup(); " | python

ADD gromacs.cfg /root/.gromacswrapper.cfg

#Install packmol
RUN apt-get install -y gfortran
ADD packmol.tar.gz /usr/share
RUN cd /usr/share/packmol && ./configure && make && ln -s /usr/share/packmol/packmol  /sbin/packmol

USER jovyan

