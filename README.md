# Reconstruction pipeline for Pulseq acquistion data using the BART toolbox
## Example sequence and data

The example sequence is located in "pulseq_sequence". Data acquired with this example sequence is located in "example_data".
For file size reasons only one slice & noise data is stored in the example data. Reconstructed data is stored in "debug/out.h5"

## Set up docker image and start server

- install docker and add user to docker group (`sudo usermod -aG docker username`)
- run `./build_docker`
- run `./start_docker`
- run `start_server` within docker container
- results (in npy format) and a log file are stored in the `debug` folder

## Reconstruction via Gadgetron ISMRMRD Client

Reconstruction can be started via the script send_data.sh:
- install Gadgetron (https://github.com/gadgetron/gadgetron)
- install ismrmrd and siemens_to_ismrmrd (https://github.com/ismrmrd/siemens_to_ismrmrd)
- usage in terminal for example data: send_data.sh "ismrmrd_file/siemens_twix_file"

## Reconstruction of Pulseq data

If you want to reconstruct Pulseq data, an additional protocol file has to be provided. This protocol file has to contain all necessary information for reconstruction.
The protocol file has to be located in "dependency/pulseq_protocols". The protocol file name has to be stored in the tFree parameter of the Siemens protocol. The necessary parameters for a reconstruction are listed in the functions "insert_hdr" and "insert_acq" in "bart_pulseq.py".
