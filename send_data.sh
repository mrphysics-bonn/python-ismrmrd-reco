#!/bin/bash

# input & output
if [ "$#" -lt 1 ]; then
    echo "Input file missing"
    exit 1
else
    IN_FILE="$1"
fi

if [ "$#" -lt 2 ]; then
    OUT_FILE="debug/out.h5"
    rm $OUT_FILE
else
    OUT_FILE="$2"
fi

PARAM_PATH=${BASH_SOURCE%/*}/python-ismrmrd-server/parameter_maps

if [[ $IN_FILE == *.dat ]]; then
    rm /tmp/tmp.h5
    siemens_to_ismrmrd -f $IN_FILE --user-map $PARAM_PATH/IsmrmrdParameterMap_Siemens_pulseq.xml --user-stylesheet $PARAM_PATH/IsmrmrdParameterMap_Siemens_pulseq.xsl -o /tmp/tmp.h5
    IN_FILE="/tmp/tmp.h5"
fi

client.py -a 127.0.0.1 -c bart_jemris -o $OUT_FILE -G images $IN_FILE
