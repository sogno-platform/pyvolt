# State Estimation Algorithm
The algorithm is implemented in Python and derived from Marco Pau's original version.  
This project uses CIMpy to read CIM data into Python objects: https://git.rwth-aachen.de/acs/core/cim/cimpy

## Installation

Install package requirements

```
pip install -r requirements.txt
```

Install state estimation package in development mode with (being in the directory of the repository)

```
python setup.py develop
```

Besides, install cimpy package in development mode with

```
cd dependencies/cimpy
python setup.py develop
```

## Getting started

To get started you find executable examples under [examples/quickstart](examples/quickstart)