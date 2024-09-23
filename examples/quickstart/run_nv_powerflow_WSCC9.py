import logging
from pathlib import Path
from pyvolt import network
from pyvolt import nv_powerflow
import numpy
import cimpy
import os
from pyvolt.network import BusType


logging.basicConfig(filename='run_nv_powerflow.log', level=logging.INFO, filemode='w')

#python starts as module in subdirectory, 2 folders up to set the new path
this_file_folder =  Path(__file__).parents[2]
p = str(this_file_folder)+"/examples/sample_data/WSCC9"
xml_path = Path(p)


xml_files = [os.path.join(xml_path, "WSCC-09_Dyn_Fourth_EQ.xml"),
             os.path.join(xml_path, "WSCC-09_Dyn_Fourth_SV.xml"),
             os.path.join(xml_path, "WSCC-09_Dyn_Fourth_TP.xml")]

# Read cim files and create new network.System object
res = cimpy.cim_import(xml_files, "cgmes_v2_4_15")
system = network.System()
base_apparent_power = 100  # MW
system.load_cim_data(res['topology'], base_apparent_power)

# Set node 1 as slack
system.nodes[4].type=BusType.SLACK

# Execute power flow analysis
results_pf, num_iter = nv_powerflow.solve(system)

# Print node voltages
print("Powerflow converged in " + str(num_iter) + " iterations.\n")
print("Results: \n")
voltages = []
for node in results_pf.nodes:
    print('Magnitude {}={}, Phase {}={}'.format(node.topology_node.uuid, numpy.abs(node.voltage_pu), node.topology_node.uuid, numpy.angle(node.voltage_pu, deg=True)))
    voltages.append(node.voltage_pu)

"""
# validation against matpower
matpower_voltages = {
    "BUS1": [1.040, 0.000],
    "BUS2": [1.025, 9.693],
    "BUS3": [1.025, 4.881],
    "BUS4": [0.996, -2.306],
    "BUS5": [0.966, -3.737],
    "BUS6": [1.004, 2.107],
    "BUS7": [0.979, 0.836],
    "BUS8": [0.997, 3.974],
    "BUS9": [0.951, 4.138]
}

# assertion
for idx, node in enumerate(results_pf.nodes):
    print(matpower_voltages[node.topology_node.uuid][0] - numpy.round(numpy.abs(node.voltage_pu),3))
    print(matpower_voltages[node.topology_node.uuid][1] - numpy.round(numpy.angle(node.voltage_pu, deg=True),3))
"""    
