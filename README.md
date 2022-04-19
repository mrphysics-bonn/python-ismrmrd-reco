# Reconstruction pipeline for Pulseq/JEMRIS data with example Sequences and Datasets

This repository contains a reconstruction pipeline for MRI data acquired with Pulseq [1]. The data is reconstructed using the BART MRI Toolbox [2].

## Quick start guide - from sequence to image

In this quick start, Cartesian and spiral sequences are created with PyPulseq. The raw data collected with this sequence, is reconstructed and images are displayed. A Python and a Docker installation are required. Information on how to install Docker can be found below.

Creating a sequence:
1. Install the necessary dependencies by running `conda env create -f seqdev.yml`, which will create a new Python environment. Activate the environment with `conda activate seqdev`.
2. Run the Python scripts "write_spiral.py" or "write_cartesian.py" in the directory "example_sequences/pypulseq". At the top of both scripts, protocol parameters and the sequence filename can be changed. 
3. A Pulseq file (.seq) is created in the same directory and an ISMRMRD metadata file (.h5) is created in the folder "dependency/metadata". This metadata file is important for the reconstruction, as raw data obtained from Pulseq sequences does not contain any information, on how the kspace was acquired.

Running a reconstruction:
1. Pull the reconstruction container from Dockerhub: `docker pull mavel101/bart-reco-server`.
2. Start the container by running `./start_docker mavel101/bart-reco-server`. The reconstruction server is now running in the background.
3. Example ISMRMRD raw data files are located in "example_data". Raw data conversion for Siemens data to ISMRMRD is described below.
4. Run a reconstruction by sending the data to the server.  
Example Pulseq reconstruction: `./send_data_pulseq.sh example_data/scanner/raw_spiralout_gre_fatsat_7T.h5 recon/out.h5`. 
Example JEMRIS reconstruction: `./send_data_jemris.sh example_data/simu/signals_spiralout_clean_slc30.h5 recon/out.h5`.  
5. The metadata is automatically merged to the raw data during the reconstruction process. Logging information and debug files can be found in the "debug" folder.
6. In this example, the reconstructed image is stored in "recon/out.h5". The image can be viewed by running the Python script "plot_img.py". Images are stored in ISMRMRD image format. Image files will not be overwritten, but new images will be appended to existing files.

JEMRIS example sequences can be found in "example_sequences/jemris". Installation instructions and documentation regarding JEMRIS can be found on the projects website: https://github.com/JEMRIS/jemris/

## Example sequences and data

The relevant files for reconstruction are placed in subfolders:  

- "example_data": Contains raw datasets from a real MR scanner and from simulation with JEMRIS [3], that can be reconstructed with the pipeline.
- "example_sequences": Contains the Pulseq sequence files, raw data was acquired with, as well as the source code for Python/PyPulseq sequences (incl. ISMRMRD [4] metadata creation) and XML files for JEMRIS sequences
- "dependency": Contains reconstruction dependencies, mainly the ISMRMRD [4] metadata files
- "recon": Contains reconstructed images in hdf5 file. 

The non-Cartesian example data provided in this repository was acquired with a spiral sequence and the reconstruction uses the GIRF predicted [5] spiral k-space trajectory.

## Set up docker image and start reconstruction server

- Install docker and add your user to the docker group (execute `sudo groupadd docker`, `sudo usermod -aG docker $USER` and `newgrp docker` after docker installation)
- A working docker image can be installed from Dockerhub with `docker pull mavel101/bart-reco-server`.

If you want to build the docker image from the latest Dockerfile in this repository, the following steps are required:
- Clone the repository and run `git submodule update --init`
- Run `./build_docker.sh` from the project folder. This builds the docker image on your system.

The default docker container contains only CPU based reconstructions. A Docker container with GPU support can be build with: `./build_docker.sh python-ismrmrd-server/ bart_cuda`
Note that this container is of larger size and that the GPU version need nvidia-docker installed (https://github.com/NVIDIA/nvidia-docker).

The container can be started by executing `./start_docker` or `./start_docker_it` from the project folder:
- `./start_docker` starts the container and runs the reconstruction server in background until it is killed with `docker kill #containerID`, where "#containerID" is the ID of the container (check with `docker ps`)
- `./start_docker_it` starts an interactive docker session in the bash shell. Run `start_server` within the container to start the reconstruction server. If you leave the session, the container is killed.
- Use `./start_docker_gpu` or `./start_docker_it_gpu` for GPU support (nvidia-docker has to be installed)

## Sending data via client

Reconstruction can be started via the provided `client.py` from the "python-ismrmrd-server" folder. A conda environment with the client can be installed, by using the provided `ismrmrd_client.yaml` using the command `conda env create -f ismrmrd_client.yml`. Afterwards it is activated by running `conda activate ismrmrd_client`.
Alternatively, the client can be installed in the current environment by running `pip install .` from the "python-ismrmrd-server" folder. The dependencies numpy, h5py and ismrmrd-python will automatically be installed. 

To run an example spiral reconstruction:

- The scripts `send_data_pulseq.sh` (data from MR scanner) and `send_data_jemris.sh` (data from JEMRIS simulations) can be used as a shortcut for sending the data. For example, the command `./send_data_pulseq.sh example_data/scanner/raw_spiralout_gre_fatsat_7T.h5 recon/out.h5` reconstructs a spiral dataset.
- By default, reconstructed data will be located in "recon/out.h5".
- Reconstructed images can be plotted with the `plot_img.py` script.
- Debug files (in npy format) and a log file are stored in the "debug" folder
- The above command extends to `python python-ismrmrd-server/client.py -c bart_pulseq -g images -o recon/out.h5 example_data/scanner/raw_spiralout_gre_fatsat_7T.h5`. The option "-c" selects the configuration name for the current reconstruction, which is evaluated in `server.py` and starts the respective reconstruction script. Available options are "bart_pulseq" for Pulseq reconstructions and "bart_jemris" for JEMRIS reconstructions. The option -o defines the image output path and option -G defines the group name of the output images in the hdf5 file.
- The script plot_img.py can be used to view reconstructed images via the Python library matplotlib

## Further information on acquisition and reconstruction
### Acquisition and reconstruction of Pulseq data

For reconstruction of Pulseq data, a ISMRMRD metadata file has to be provided. This metadata file has to contain all necessary information for reconstruction such as counters, flags and other metadata. The metadata file has to be located in "dependency/metadata". 

If the sequence is executed on a Siemens scanner, the following steps are necessary:
- The metadata file name has to be stored in the free text parameter "tFree" of the raw data file. It will be converted to the first user defined string parameter in the ISMRMRD file.
- The raw data has to be converted to the ISMRMRD format with the siemens_to_ismrmrd converter (https://github.com/ismrmrd/siemens_to_ismrmrd). After installation of the converter, `send_data_pulseq.sh` can handle Siemens raw data acquired with the Pulseq sequence.
- The metadata filename is transferred to the first user defined string parameter of the ISMRMRD file by the converter. The parameter maps used for conversion are stored in "python-ismrmrd-server/parameters_maps".
- The necessary parameters for reconstruction are listed in the functions "insert_hdr" and "insert_acq" in "python-ismrmrd-server/pulseq_helper.py".

### Reconstruction of JEMRIS simulation data

Reconstruction of JEMRIS simulation data can be started within JEMRIS by selecting the "-r" option in the simulation or by starting the BART recon in the GUI. However, the following prerequisites have to be met:
- The reconstruction server has to be running.
- The `client.py` (and its dependencies, see above) has to be installed in the conda environment, where the JEMRIS simulation is executed. Install the client by running `pip install .` from the "python-ismrmrd-server" folder. This lets you execute the client from anywhere. It is recommended to use the provided `ismrmrd_client.yml` to create a conda environment for the client.

Reconstruction of already simulated data can also be started by running `send_data_jemris.sh` as described in "Sending data via client".

### Extend & Modify existing reconstruction

If the scripts `./start_docker` or `./start_docker_it` were used to start the docker container, the reconstruction codebase in the `python-ismrmrd-server` subdirectory is mounted in the active container. Changes to the reconstruction scripts will immediately applied in a new reconstruction started. If subscripts are altered, a restart of the reconstruction server might be necessary.
Some important scripts are explained in more detail:  
- `server.py`: New reconstruction scripts can be added here and a new configuration name should be assigned. The new reconstruction can be started sending data via `client.py` with the new configuration name by using the "-c" option (see above).
- `pulseq_helper.py`: Contains the funtions for mirroring metadata information from the metadata file to the streamed raw data.
- `bart_pulseq.py`: Launches a BART reconstruction pipeline for Pulseq data, depending on the trajectory type 

More reconstruction scripts and Dockerfiles including the PowerGrid reconstruction toolbox [6] for static offresonance correction support can be found in the "pulseq" branch of the sub-repository `python-ismrmrd-server` (https://github.com/pehses/python-ismrmrd-server/tree/pulseq).

## References

[1] Layton, K. J. et. al. Pulseq: A rapid and hardware-independent pulse sequence prototyping framework, MRM, 2017;77(4):1544-1552, http://pulseq.github.io/

[2] BART Toolbox for Computational Magnetic Resonance Imaging, DOI: 10.5281/zenodo.592960, https://mrirecon.github.io/bart

[3] Inati, J. I. et. al. ISMRM Raw data format: A proposed standard for MRI raw datasets, MRM, 2017;77(1):411-421, https://ismrmrd.github.io

[4] Stöcker, T. et. al. High-Performance Computing MRI Simulations, MRM, 2010;64:186-193, https://www.jemris.org/

[5] Vannesjo, S. J. et al. Gradient System Characterization by Impulse Response Measurements with a Dynamic Field Camera. MRM
2013;69:583-593

[6] Cerjanic, A. et al. PowerGrid: A open source library for accelerated iterative magnetic resonance image reconstruction, Proc. Intl. Soc. Mag. Res. Med., 2016, http://mrfil.github.io/PowerGrid/
