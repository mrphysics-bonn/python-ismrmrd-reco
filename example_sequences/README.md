## Sequence and MRD metadata creation with Python/PyPulseq and JEMRIS

This folder contains sequence files (.seq) for execution with Pulseq, as well as the JEMRIS (.xml) source files and the Python source code sequence and MRD metadata creation. JEMRIS can be found at: https://github.com/JEMRIS/jemris/.  

A Python environment with all dependencies for running the scripts write_cartesian.py (2-echo GRE) & write_spiral.py (2D spiral GRE) can be installed by running `conda env create -f seqdev.yml` from the top directory.

In this environment, a custom PyPulseq version is installed (https://github.com/mavel101/pypulseq/tree/dev_mv).
