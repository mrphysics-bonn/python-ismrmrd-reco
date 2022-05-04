# Reconstruction pipeline for MRI raw data acquired with Pulseq or JEMRIS

This repository contains a reconstruction pipeline for MRI data acquired with Pulseq [1] or with the MR simulator JEMRIS [2]. The data is reconstructed using the BART MRI Toolbox [3].

## Installation

A [Python](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and a [Docker](https://docs.docker.com/get-docker/) installation are required to run the reconstruction server. Additionally, the following steps have to be done:
1. After Docker installation, add your user to the docker group (execute `sudo groupadd docker`, `sudo usermod -aG docker $USER` and `newgrp docker`)
2. Pull the Docker image of the reconstruction server from Dockerhub: `docker pull mavel101/bart-reco-server`.  
3. Python dependencies for the reconstruction are defined in "ismrmrd_client.yml". Run `conda env create -f ismrmrd_client.yml` to install the environment.  

Optional dependencies for the creation of MR sequences with PyPulseq [4] are defined in "seqdev.yml". Run `conda env create -f seqdev.yml` to install a new Python environment or run `conda env update -f seqdev.yml -n ismrmrd_client` to add the dependencies to the "ismrmrd_client" environment.

## Quick start guide - from sequence to image

In this quick start, Cartesian and spiral sequences are created with PyPulseq. The raw data collected with this sequence, is reconstructed and images are displayed. 

Creating a sequence:
1. Activate the Python environment with `conda activate seqdev`.
2. Run the Python scripts with `python write_spiral.py` or `python write_cartesian.py` in the directory "example_sequences/pypulseq". At the top of both scripts, protocol parameters and the sequence filename can be changed. 
3. A Pulseq file (.seq) is created in the same directory and an MRD (originally ISMRMRD) metadata file (.h5) is created in the folder "dependency/metadata". This metadata file is important for the reconstruction, as raw data obtained from Pulseq sequences does not contain any information, on how the kspace was acquired.

Running a reconstruction:
1. Start the Docker container by running `./start_docker mavel101/bart-reco-server`. The reconstruction server is now running in the background.
2. Activate the Python environment with `conda activate ismrmrd_client`.
3. Run a reconstruction by sending the data to the server.  
Example Pulseq reconstruction: `./send_data_pulseq.sh example_data/scanner/raw_spiralout_gre_fatsat_7T.h5 recon/out.h5`. 
Example JEMRIS reconstruction: `./send_data_jemris.sh example_data/simu/signals_spiralout_clean_slc30.h5 recon/out.h5`.  
5. Logging information and debug files can be found in the "debug" folder.
6. In this example, the reconstructed image is stored in "recon/out.h5". The image can be viewed by running the Python script "plot_img.py". Images are stored in MRD image format. Image files will not be overwritten, but new images will be appended to existing files.
7. More example raw data files are located in "example_data". Raw data conversion for Siemens data to MRD [5] is described below.

## Example sequences and data

The relevant files for reconstruction are placed in subfolders:  

- "example_data": Contains raw datasets from real MR scanners and from simulation with JEMRIS [2], that can be reconstructed with the pipeline.
- "example_sequences": Contains the Pulseq sequence files, raw data was acquired with, as well as the source code for Python/PyPulseq sequences (incl. MRD metadata creation) and XML files for JEMRIS sequences
- "dependency": Contains reconstruction dependencies, mainly the MRD metadata files
- "recon": Contains reconstructed images in hdf5 file. 

The non-Cartesian example data provided in this repository was acquired with a spiral sequence and the reconstruction uses the GIRF predicted [6] spiral k-space trajectory.  
Additional example sequences, installation instructions and documentation regarding JEMRIS can be found on the projects website: https://github.com/JEMRIS/jemris/.  
Pulseq sequence files can be converted to the GE compatible format TOPPE, using the converter at https://github.com/toppeMRI/PulseGEq.

## Further information on the reconstruction server
### Docker image (latest build, GPU support, startup scripts)

If you want to build the docker image from the latest Dockerfile in this repository, the following steps are required:
- Clone the repository and run `git submodule update --init`
- Run `./build_docker.sh` from the project folder. This builds the docker image on your system.

The default docker image contains only CPU based reconstructions. A Docker image with GPU support can be build with: `./build_docker.sh python-ismrmrd-server/ bart_cuda`
Note that this image is of larger size and that the GPU version needs nvidia-docker to be installed (https://github.com/NVIDIA/nvidia-docker).

The container can be started by executing `./start_docker` or `./start_docker_it` from the project folder:
- `./start_docker` starts the container and runs the reconstruction server in background until it is killed with `docker kill #containerID`, where "#containerID" is the ID of the container (check with `docker ps`)
- `./start_docker_it` starts an interactive docker session in the bash shell. Run `start_server` within the container to start the reconstruction server. If you leave the session, the container is killed.
- Use `./start_docker_gpu` or `./start_docker_it_gpu` for GPU support (nvidia-docker has to be installed)

### Sending data via the client

Reconstruction can be started via the provided `client.py` from the "python-ismrmrd-server" folder. A conda environment with the client can be installed, by using the provided `ismrmrd_client.yml` using the command `conda env create -f ismrmrd_client.yml`. Afterwards it is activated by running `conda activate ismrmrd_client`.
Alternatively, the client can be installed in the current environment by running `pip install .` from the "python-ismrmrd-server" folder. The dependencies numpy, h5py and ismrmrd-python will automatically be installed. 

To run an example spiral reconstruction:

- The scripts `send_data_pulseq.sh` (data from MR scanner) and `send_data_jemris.sh` (data from JEMRIS simulations) can be used as a shortcut for sending the data. For example, the command `./send_data_pulseq.sh example_data/scanner/raw_spiralout_gre_fatsat_7T.h5 recon/out.h5` reconstructs a spiral dataset.
- By default, reconstructed data will be located in "recon/out.h5".
- Reconstructed images can be plotted with the `plot_img.py` script.
- Debug files (in npy format) and a log file are stored in the "debug" folder
- The above command extends to `python python-ismrmrd-server/client.py -c bart_pulseq -g images -o recon/out.h5 example_data/scanner/raw_spiralout_gre_fatsat_7T.h5`. The option "-c" selects the configuration name for the current reconstruction, which is evaluated in `server.py` and starts the respective reconstruction script. Available options are "bart_pulseq" for Pulseq reconstructions and "bart_jemris" for JEMRIS reconstructions. The option -o defines the image output path and option -G defines the group name of the output images in the hdf5 file.
- The script plot_img.py can be used to view reconstructed images via the Python library matplotlib

### Reconstruction of Pulseq data

For reconstruction of Pulseq data, a MRD metadata file has to be provided. This metadata file has to contain all necessary information for reconstruction such as counters, flags and other metadata. The metadata file has to be located in "dependency/metadata".  
IMPORTANT: The metadata filename has to be saved in the first user defined string parameter ("userParameterString") of the raw data ISMRMRD file in order to access the metadata in the reconstruction.  

For Siemens data this is can automatically be done at file conversion:
1. The metadata file should have the same name as the sequence.
2. The sequence/metadata filename has to be stored in the free text parameter "tFree" of the raw data file. This is automatically done in the newest Pulseq interpreter sequence (v1.4.0).
3. The raw data is converted to the MRD format with the siemens_to_ismrmrd converter (https://github.com/ismrmrd/siemens_to_ismrmrd). After installation of the converter, `send_data_pulseq.sh` can handle Siemens raw data acquired with the Pulseq sequence. The shell script uses the parameter maps in "python-ismrmrd-server/parameter_maps" for file conversion.
4. Metadata merging is done by the functions "insert_hdr" and "insert_acq" in "python-ismrmrd-server/pulseq_helper.py".

### Reconstruction of JEMRIS simulation data

Reconstruction of JEMRIS simulation data can be started within JEMRIS by selecting the "-r" option, when running JEMRIS in the command line or by starting the recon in the GUI with the "start reco" button. However, the following prerequisites have to be met:
- The Docker image of the reconstruction server has to be pulled from Dockerhub (`docker pull mavel101/bart-reco-server`). If JEMRIS is running from the command line, the reconstruction server also has to be started.
- The `client.py` (and its dependencies) has to be installed in the conda environment, where the JEMRIS simulation is executed. This lets you execute the client from anywhere. It is recommended to use the provided `ismrmrd_client.yml` to create a conda environment for the client. Manual installation of the client is possible by running `pip install .` from the "python-ismrmrd-server" folder.

Reconstruction of already simulated data can also be started by running `send_data_jemris.sh` as described in "Sending data via client".

### Extend & Modify existing reconstruction

If the scripts `./start_docker` or `./start_docker_it` were used to start the docker container, the reconstruction codebase in the `python-ismrmrd-server` subdirectory is mounted in the active container. Changes to the reconstruction scripts will immediately applied in a new reconstruction started. If subscripts are altered, a restart of the reconstruction server might be necessary.
Some important scripts are explained in more detail:  
- `server.py`: New reconstruction scripts can be added here and a new configuration name should be assigned. The new reconstruction can be started sending data via `client.py` with the new configuration name by using the "-c" option (see above).
- `pulseq_helper.py`: Contains the funtions for mirroring metadata information from the metadata file to the streamed raw data.
- `bart_pulseq.py`: Launches a BART reconstruction pipeline for Pulseq data, depending on the trajectory type 

## Static offresonance correction support

More reconstruction scripts and Dockerfiles including the PowerGrid reconstruction toolbox [7] for static offresonance correction support can be found in the "pulseq" branch of the sub-repository `python-ismrmrd-server` (https://github.com/pehses/python-ismrmrd-server/tree/pulseq).

## References

[1] Layton, K. J. et. al. Pulseq: A rapid and hardware-independent pulse sequence prototyping framework, MRM, 2017;77(4):1544-1552, http://pulseq.github.io/

[2] Stöcker, T. et. al. High-Performance Computing MRI Simulations, MRM, 2010;64:186-193, https://www.jemris.org/

[3] BART Toolbox for Computational Magnetic Resonance Imaging, DOI: 10.5281/zenodo.592960, https://mrirecon.github.io/bart

[3] BART Toolbox for Computational Magnetic Resonance Imaging, DOI: 10.5281/zenodo.592960, https://mrirecon.github.io/bart

[4] Ravi, Keerthi, Sairam Geethanath, and John Vaughan. "PyPulseq: A Python Package for MRI Pulse Sequence Design." Journal of Open Source Software 4.42 (2019): 1725., https://github.com/imr-framework/pypulseq

[5] Inati, J. I. et. al. ISMRM Raw data format: A proposed standard for MRI raw datasets, MRM, 2017;77(1):411-421, https://ismrmrd.github.io

[6] Vannesjo, S. J. et al. Gradient System Characterization by Impulse Response Measurements with a Dynamic Field Camera. MRM
2013;69:583-593

[7] Cerjanic, A. et al. PowerGrid: A open source library for accelerated iterative magnetic resonance image reconstruction, Proc. Intl. Soc. Mag. Res. Med., 2016, http://mrfil.github.io/PowerGrid/
