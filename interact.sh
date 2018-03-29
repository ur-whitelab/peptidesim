#!/bin/bash


echo "In container peptidesim directroy run: pip install -e ."
echo "Then change to /home/whitelab/scratch and run"
echo "~/.local/bin/pytest ../peptidesim/package"
echo "No need to rebuild container or re-install. Just edit code"
docker run --rm -it -v `pwd`:/home/whitelab/peptidesim peptidesim/test bash
