#!/bin/bash

docker run -it -u tester -v "`pwd`:/home/tester/peptidesim" peptidesim/test
