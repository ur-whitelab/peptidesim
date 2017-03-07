#!/bin/bash

docker run -it -u tester -v "`pwd`:/usr/share/peptidesim" peptidesim/test
