spiraltraj
======================================
spiraltraj is a python wrapper of a C/C++ spiral MRI gradient trajectory optimization code (gen_spiral.cpp contains the python interface).

vdspiral calculates a simple Archimedean spiral (with possible variable density adjustments) as a parameterized curve that is then passed to a gradient optimization algorithm.

The gradient optimization is performed using "Time Optimal Gradient Design" C-Code published by Miki Lustig ("https://people.eecs.berkeley.edu/~mlustig/Software.html", files: mtg_functions.h, spline.h). I lazily and incompletely translated the code to C++ a long time ago (this was a bit unnecessary, really).

The gradient optimization algorithm is published in:
Lustig M, Seung-Jean K, Pauly JM, IEEE Trans Med Imaging. 2008 Jun;27(6):866-73. doi: 10.1109/TMI.2008.922699.

## Author
Questions & feedback go to philipp.ehses@dzne.de

## Dependencies
- python 3.7 (previous version may or may not work)
- python development header files ("Python.h" & dependencies)
- gcc (other compilers may or may not work)

## Installation
Simply clone the repository and run
```
python setup.py install
```
from spiraltraj's base directory.

## Usage

```python
import spiraltraj

help(spiraltraj)
```
should get you started.