#!/bin/bash

# input & output
if [ "$#" -lt 1 ]; then
    echo "Input file missing"
    exit 1
else
    IN_FILE="$1"
fi

if [ "$#" -lt 2 ]; then
    OUT_FILE="recon/out.h5"
    rm $OUT_FILE
else
    OUT_FILE="$2"
fi

python python-ismrmrd-server/client.py -c powergrid_pulseq -o $OUT_FILE -G images $IN_FILE
