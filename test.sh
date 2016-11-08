#!/bin/bash

docker build -t peptidesim/test test-docker
docker run -v `pwd`:/usr/share/peptidesim peptidesim/test
