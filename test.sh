#!/bin/bash

RESULT=$(docker run -v "`pwd`:/usr/share/peptidesim" peptidesim/test)
docker stop $RESULT
