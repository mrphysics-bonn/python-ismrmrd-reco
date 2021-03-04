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
else
    OUT_FILE="$2"
fi

PARAM_PATH=${BASH_SOURCE%/*}/python-ismrmrd-server/parameter_maps

if [[ $IN_FILE == *.dat ]]; then
    rm /tmp/tmp.h5
    siemens_to_ismrmrd -f $IN_FILE --user-map $PARAM_PATH/IsmrmrdParameterMap_Siemens_pulseq.xml --user-stylesheet $PARAM_PATH/IsmrmrdParameterMap_Siemens_pulseq.xsl -o /tmp/tmp.h5
    IN_FILE="/tmp/tmp.h5"
fi

gadgetron_ismrmrd_client -a 127.0.0.1 -c bart_pulseq -f $IN_FILE -o $OUT_FILE
