# gmsh4mrst

Automatically create Gmsh meshes for use in MRST.

`gmsh4mrst` is created to bridge the gap between Python package Gmsh and MATLAB module MRST, and was developed as part of my Bachelor's thesis during spring 2022.

## Citation
If you end up using `gmsh4mrst` in any research, I would love to hear about it. Feel free to cite the bachelor thesis introducing the software as well:
Berg, Andreas B. (2022). _Combining Gmsh and MRST - Developing software for more efficient grid creation in two dimensions_

The thesis can be found [here](https://github.com/BollaBerg/MA2002-Bachelor-project/blob/main/Berg%2C%20Andreas%20Bjelland%20(2022)%20-%20Combining%20Gmsh%20and%20MRST.pdf).

## Installation
The Python part of the module is available on [PyPi](https://pypi.org/project/gmsh4mrst/).

The MATLAB part of the module must be manually added to your MATLAB Path.

## Requirements
In order to use the Python part of the module, simply install the package (`pip install gmsh4mrst`).

In order to use the MATLAB part of the module, ensure `MRST 2022a` is on your MATLAB Path, with module `gmsh` active, and install the Python part of the module. Ensure that your default Python (i.e. the Python that comes up when using the command `python`) has `gmsh4mrst` installed, and is one of the Python versions compitable with MATLAB (i.e. `Python 3.8` or `Python 3.9` when using MATLAB 2022a).

## Common issues
### MATLAB
#### SEGFAULT
This sometimes happens when calling Python from MATLAB. One thing that has previously worked to fix this is ensuring Python is run Out-of-process, rather than In-process. This can be ensured by running `pyenv("ExecutionMode","OutOfProcess")` from within MATLAB. Note that this must be done before running any Python code from MATLAB, so may require a restart of MATLAB.