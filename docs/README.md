# Documentation

A good place to start is the examples/quickstart/ folder of the repo.

## Inputs

Basically, pyvolt requires to kinds of data:
- grid topology
- measurements or load flow results
- uncertainty values for measurements

The grid data can be imported from CIM / CGMES XML files (e.g. examples/sample_data/CIGRE-MV-NoTap/*.xml).
cimpy is used for the import of grid data.

Measurement values can be read from a CSV file (e.g. examples/sample_data/CIGRE-MV-NoTap/*.csv).

Uncertainties are specified per measurement type. The measurement types are defined in: pyvolt/measurement.py -> class MeasType

## Outputs

The state estimation returns a Results object (pyvolt/results.py).
While voltages are available upon completion of the state estimation, other variables like currents and power can be calculated on demand if required.