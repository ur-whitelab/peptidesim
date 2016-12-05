#!/bin/bash

cp -R /usr/share/peptidesim /root
cd /root/peptidesim/package
pip install -r requirements.txt
tox
#install for interactive use
source activate python2
pip install . --user
tail -f /dev/null
