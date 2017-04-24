#!/bin/bash

RESULT=$(docker run -d -v "`pwd`:/home/tester/peptidesim" peptidesim/test)
echo "run the following: "
echo "[wait until docker container is finished (check logs)]"
echo "sudo docker logs $RESULT"
echo "[enter the docker container]"
echo "sudo docker exec -i -t $RESULT /bin/bash"
echo "[switch to python 2]"
echo "source activate python2"
echo "[make it possible to use text-editor]"
echo "export TERM=linux"
echo "after exiting docker, be sure to terminate it!"
