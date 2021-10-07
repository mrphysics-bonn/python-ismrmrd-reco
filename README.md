# Reconstruction pipeline for Pulseq acquistion data using the BART toolbox
## Example sequence and data

This repository contains a reconstruction pipeline for MRI data acquired with Pulseq [1]. The data is reconstructed using the BART MRI Toolbox [2].
The example Pulseq sequences are located in the folder "pulseq_sequence". Data acquired with this example sequence is stored in the ISMRMRD format [3] in "example_data". The non-Cartesian example data provided in this repository was acquired with a spiral sequence and the reconstruction uses the GIRF predicted [4] spiral k-space trajectory.
For file size reasons the example data contains only one slice & noise data. Reconstructed data will be located in "debug/out.h5".

## Set up docker image and start reconstruction server

- Clone the repository and run `git submodule update --init`
- Install docker and add user to docker group (execute `sudo groupadd docker`, `sudo usermod -aG docker $USER` and `newgrp docker` after docker installation)
- Run `./build_docker.sh` from the project folder. This builds the docker image on your system.
- Alternatively the docker image can be installed from Dockerhub with `docker pull mavel101/bart-reco-server`.

The default docker container contains only CPU based reconstructions. A Docker container with GPU support can be build with: `./build_docker.sh python-ismrmrd-server/ bart_cuda`
Note that this container are of larger size and that the GPU versions need nvidia-docker installed (https://github.com/NVIDIA/nvidia-docker).

The container can be started by executing `./start_docker` or `./start_docker_it` from the project folder:
- `./start_docker` starts the container and runs the reconstruction server in background until it is killed with `docker kill #containerID`, where "#containerID" is the ID of the container (check with `docker ps`)
- `./start_docker_it` starts an interactive docker session in the bash shell. Run `start_server` within the container to start the reconstruction server. If you leave the session, the container is killed.
- Use `./start_docker_gpu` or `./start_docker_it_gpu` for GPU support (nvidia-docker has to be installed)

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

## Extend & Modify existing reconstruction

If the scripts `./start_docker` or `./start_docker_it` were used to start the docker container, the reconstruction codebase in the `python-ismrmrd-server` subdirectory is mounted in the active container. Changes to the reconstruction scripts will immediately applied in a new reconstruction started. If subscripts are altered, a restart of the reconstruction server might be necessary.
New reconstruction scripts can be added in `server.py`, where a new configuration name should be assigned. The new reconstruction can be started sending data via `client.py` with the new configuration name by using the "-c" option (see above).

More reconstruction scripts and Dockerfiles (e.g. for the PowerGrid reconstruction toolbox [5] for B0 correction support) can be found in the original repository: https://github.com/pehses/python-ismrmrd-server/tree/pulseq

## References

[1] Layton, K. J. et. al. Pulseq: A rapid and hardware-independent pulse sequence prototyping framework, MRM, 2017;77(4):1544-1552, http://pulseq.github.io/

[2] BART Toolbox for Computational Magnetic Resonance Imaging, DOI: 10.5281/zenodo.592960, https://mrirecon.github.io/bart

[3] Inati, J. I. et. al. ISMRM Raw data format: A proposed standard for MRI raw datasets, MRM, 2017;77(1):411-421, https://ismrmrd.github.io

[4] Vannesjo, S. J. et al. Gradient System Characterization by Impulse Response Measurements with a Dynamic Field Camera. MRM
2013;69:583-593

[5] Cerjanic, A. et al. PowerGrid: A open source library for accelerated iterative magnetic resonance image reconstruction, Proc. Intl. Soc. Mag. Res. Med., 2016
