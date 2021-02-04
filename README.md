# PyVolt
 
This package provides algorithms for running powerflow simulations, state estimation and FLISR in Python. <br>
It uses [CIMpy](https://github.com/cim-iec/cimpy) to read network data based on the Common Information Model.

## Installation

### User 

Install the pyvolt package with

    $ python setup.py install

### Developer

Install the pyvolt package in development mode with

    $ python setup.py develop

### Docker setup

Docker is a quick way to set up a development environment in case you do not have a python distribution installed.

Start a docker container that mounts the cloned files:

On Windows

    $ docker run -it -v ${pwd}:/pyvolt python bash  

On Linux

    $ docker run -it -v $(pwd):/pyvolt python bash 

Then, you can change into the pyvolt folder and run the setup file as described above

    $ cd /pyvolt
    $ python setup.py develop


## Getting started

To get started you find executable examples under [examples/quickstart](examples/quickstart)