# Reconstruction pipeline for Pulseq acquistion data using the BART toolbox
## Example sequence and data

This repository contains a reconstruction pipeline for MRI data acquired with Pulseq [1]. The data is reconstructed using the BART MRI Toolbox [2]. The data was acquired with a spiral sequence and the reconstruction uses the GIRF predicted [3] spiral k-space trajectory.
The example Pulseq sequence is located in the folder "pulseq_sequence". Data acquired with this example sequence is stored in the ISMRMRD format [4] in "example_data".
For file size reasons the example data contains only one slice & noise data. Reconstructed data is located in "debug/out.h5"

## Set up docker image and start server
- Clone the repository and do `git submodule update --init`
- Install docker and add user to docker group (`sudo groupadd docker`, `sudo usermod -aG docker $USER` and `newgrp docker`)
- Run `./build_docker.sh` from the project folder
- Run `./start_docker` from the project folder
- Run `start_server` within docker container
- Results (in npy format) and a log file are stored in the `debug` folder

## Sending data via Client

Reconstruction can be started via 2 options:
1. Use the provided `client.py` from the python-ismrmrd-server folder:
2. Use the Gadgetron ISMRMRD Client (Gadgetron has to be installed: https://github.com/gadgetron/gadgetron)

Via `client.py`:
- Run `python client.py -c bart_pulseq ../example_data/pulseq_gre_dataset.h5`. The option "-c" submits the configuration for the current reconstruction, which is evaluated in `server.py` and starts the respective reconstruction script.

Via Gadgetron ISMRMRD client:
- Reconstruction can be started via the script `send_data.sh`, run `./send_data.sh example_data/pulseq_gre_dataset.h5` from the project folder
- Optional for Siemens Twix files: Install ismrmrd and siemens_to_ismrmrd (https://github.com/ismrmrd/siemens_to_ismrmrd) for Siemens datasets. After installation, send_data.sh can also handle Siemens raw data acquired with the Pulseq sequence.
## Reconstruction of Pulseq data

For reconstruction of Pulseq data, an additional protocol file has to be provided. This protocol file has to contain all necessary information for reconstruction such as counters, flags and other metadata.
The protocol file has to be located in "dependency/pulseq_protocols". The protocol file name has to be stored in the free text parameter "tFree" of the protocol, if the sequence is executed on a Siemens scanner. The protocol name is transferred to the first user defined string parameter of the ISMRMRD file by the siemens_to_ismrmrd convertert. The necessary parameters for a reconstruction are listed in the functions "insert_hdr" and "insert_acq" in "bart_pulseq.py".

## References

[1] Layton, K. J. et. al. Pulseq: A rapid and hardware-independent pulse sequence prototyping framework, MRM, 2017;77(4):1544-1552, http://pulseq.github.io/

[2] BART Toolbox for Computational Magnetic Resonance Imaging, DOI: 10.5281/zenodo.592960, https://mrirecon.github.io/bart

[3] Vannesjo, S. J. et al. Gradient System Characterization by Impulse Response Measurements with a Dynamic Field Camera. MRM
2013;69:583-593

[4] Inati, J. I. et. al. ISMRM Raw data format: A proposed standard for MRI raw datasets, MRM, 2017;77(1):411-421, https://ismrmrd.github.io
