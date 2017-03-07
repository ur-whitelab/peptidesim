#!/bin/bash

cp -R /usr/share/peptidesim /home/tester
cd /home/tester/peptidesim/package
tox
#install for interactive use
#pip install . --user
#tail -f /dev/null
