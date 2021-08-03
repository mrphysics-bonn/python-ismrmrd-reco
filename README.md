# Reconstruction pipeline for Pulseq acquistion data using the BART toolbox
## Example sequence and data

This repository contains a reconstruction pipeline for MRI data acquired with Pulseq [1]. The data is reconstructed using the BART MRI Toolbox [2]. The example data provided in this repository was acquired with a spiral sequence and the reconstruction uses the GIRF predicted [3] spiral k-space trajectory.
The example Pulseq sequence is located in the folder "pulseq_sequence". Data acquired with this example sequence is stored in the ISMRMRD format [4] in "example_data".
For file size reasons the example data contains only one slice & noise data. Reconstructed data is located in "debug/out.h5"

## Set up docker image and start reconstruction server

- Clone the repository and run `git submodule update --init`
- Install docker and add user to docker group (execute `sudo groupadd docker`, `sudo usermod -aG docker $USER` and `newgrp docker` after docker installation)
- Run `./build_docker.sh` from the project folder. This builds the docker image on your system.

The container can be started by executing `./start_docker` or `./start_docker_it` from the project folder:
- `./start_docker` starts the container and runs the reconstruction server in background until it is killed with `docker kill #containerID`, where "#containerID" is the ID of the container (check with `docker ps`)
- `./start_docker_it` starts an interactive docker session in the bash shell. Run `start_server` within the container to start the reconstruction server. If you leave the session, the container is killed.

## Sending data via client

Reconstruction can be started via the provided `client.py` from the "python-ismrmrd-server" folder:

- Run `python client.py -c bart_pulseq ../example_data/pulseq_gre_dataset.h5`. The option "-c" submits the configuration for the current reconstruction, which is evaluated in `server.py` and starts the respective reconstruction script.
- The scripts `send_data_pulseq.sh` and `send_data_jemris.sh` can be used for sending data. For example, the above command reduces to `./send_data_pulseq.sh example_data/pulseq_gre_dataset.h5`.
- Required only for reconstruction of JEMRIS simulation data: Install the client in your conda environment by running `pip install .` from the "python-ismrmrd-server" folder. This lets you execute the client from anywhere.
- Debug files (in npy format) and a log file are stored in the "debug" folder

## Reconstruction of Pulseq data

For reconstruction of Pulseq data, an additional protocol file has to be provided. This protocol file has to contain all necessary information for reconstruction such as counters, flags and other metadata. The protocol file has to be located in "dependency/pulseq_protocols". 

If the sequence is executed on a Siemens scanner, the following steps are necessary:
- The protocol file name has to be stored in the free text parameter "tFree" of the raw data protocol.
- The raw data has to be converted to the ISMRMRD format with the siemens_to_ismrmrd converter (https://github.com/ismrmrd/siemens_to_ismrmrd). After installation of the converter, `send_data_pulseq.sh` can handle Siemens raw data acquired with the Pulseq sequence.
- The protocol name is transferred to the first user defined string parameter of the ISMRMRD file by the converter. The parameter maps used for conversion are stored in "python-ismrmrd-server/parameters_maps".
- The necessary parameters for reconstruction are listed in the functions "insert_hdr" and "insert_acq" in "python-ismrmrd-server/pulseq_prot.py".

## Reconstruction of JEMRIS data

Reconstruction of JEMRIS simulation data can be started within JEMRIS by selecting the "-r" option in the simulation. However, the following prerequisites have to be met:
- The reconstruction server has to be running.
- The `client.py` has to be installed in the conda environment, where the JEMRIS simulation is executed.

Reconstruction of already simulated data can also be started by running `send_data_jemris.sh` as described above.

## References

[1] Layton, K. J. et. al. Pulseq: A rapid and hardware-independent pulse sequence prototyping framework, MRM, 2017;77(4):1544-1552, http://pulseq.github.io/

[2] BART Toolbox for Computational Magnetic Resonance Imaging, DOI: 10.5281/zenodo.592960, https://mrirecon.github.io/bart

[3] Vannesjo, S. J. et al. Gradient System Characterization by Impulse Response Measurements with a Dynamic Field Camera. MRM
2013;69:583-593

[4] Inati, J. I. et. al. ISMRM Raw data format: A proposed standard for MRI raw datasets, MRM, 2017;77(1):411-421, https://ismrmrd.github.io
