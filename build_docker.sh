#!/bin/bash


if [ "$#" -lt 1 ]; then
    PATH="./python-ismrmrd-server"
else
    PATH="$1"
fi

if [ "$#" -lt 2 ]; then
    TAG="bart_server"
else
    TAG="$2"
fi

cd $PATH
cd docker


/usr/bin/docker build --no-cache -t $TAG .
