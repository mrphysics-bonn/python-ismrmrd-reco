## Sequence and ISMRMRD protocol creation with Python/PyPulseq and JEMRIS

This folder contains sequence files (.seq) for execution with Pulseq, as well as the JEMRIS (.xml) source files and the Python source code sequence and ISMRMRD protocol creation. JEMRIS can be found at: https://github.com/JEMRIS/jemris/tree/ismrmrd-export.  

For the Python scripts write_cartesian.py (B0 mapping) & write_spiral.py the following dependencies are needed:
- Custom version of PyPulseq:
```console
git clone https://github.com/mavel101/pypulseq
cd pypulseq
git switch dev_mv
git checkout 7b29514
pip install .
```
- ISMRMRD Python API:
```console
git clone https://github.com/ismrmrd/ismrmrd-python
cd ismrmrd-python
git checkout v1.9.3
pip install .
```
- spiraltraj for the Spiral sequence (write_spiral.py):
```console
git clone https://github.com/mrphysics-bonn/spiral-pypulseq-example
cd spiral-pypulseq-example
cd spiraltraj
pip install .
```
A shorter Python example without installing dependencies can be found at: https://github.com/mrphysics-bonn/spiral-pypulseq-example.
