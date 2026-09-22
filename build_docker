#!/bin/bash


if [ "$#" -lt 1 ]; then
    PATH="./python-ismrmrd-server"
else
    PATH="$1"
fi

if [ "$#" -lt 2 ]; then
    TAG="bart"
else
    TAG="$2"
fi

cd $PATH
cd docker
cd $TAG

/usr/bin/docker build -t $TAG .
