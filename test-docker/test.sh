#!/bin/bash

cp -R ~/peptidesim ~/scratch
cd ~/scratch/peptidesim/package
tox
#install for interactive use
#pip install . --user
#tail -f /dev/null
