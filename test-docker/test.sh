#!/bin/bash

cp -R /usr/share/peptidesim /root
cd /root/peptidesim/package
pip install -r requirements.txt
tox
#tail -f /dev/null
