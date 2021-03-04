## Linux host machine
- install docker and add user to docker group (`sudo usermod -aG docker username`)
- run `./build_docker`
- run `./start_docker`
- run `start_server` within docker container
- results (in npy format) and a log file are stored in the `debug` folder

## Gadgetron ISMRMRD Client
Reconstruction can also be started via the script send_data.sh without the IDEA virtual machine (more convenient):  
- install Gadgetron from https://github.com/gadgetron/gadgetron  
- Set correct location of ParameterMaps for siemens_to_ismrmrd conversion in send_data.sh
- set default out file path in send_data.sh
- usage in terminal: send_data.sh "siemens_twix_file" "reco_out_file"

## Reconstruction of Pulseq data
If you want to reconstruct Pulseq data, an additional protocol has to be provided. This protocol file has to contain all necessary information for reconstruction.  
The protocol file has to be located in "dependency/pulseq_protocols". The necessary parameters for a reconstruction are listed in the functions "insert_hdr" and "insert_acq" in "bart_pulseq.py".