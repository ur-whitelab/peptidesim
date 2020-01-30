#!/bin/bash
docker run --rm -it -v `pwd`:/home/whitelab/peptidesim peptidesim/test bash ../interact.sh
