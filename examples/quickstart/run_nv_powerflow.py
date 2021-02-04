import logging
from pyvolt import network
from pyvolt import nv_powerflow
from pyvolt import results
import cimpy
import os
from pathlib import Path

logging.basicConfig(filename='CIGRE.log', level=logging.INFO, filemode='w')

cim_dir = Path('.') / 'sample_data' / 'CIGRE-MV-NoTap'
cim_files = cim_dir.glob('*.xml')
cim_list = []
for file in cim_files:
    cim_list.append(str(file.absolute()))
print(cim_list)

# read cim files and create new network.System object
res, _, _ = cimpy.cim_import(cim_list, "cgmes_v2_4_15")
system = network.System()
base_apparent_power = 25  # MW
system.load_cim_data(res, base_apparent_power)

# Execute power flow analysis
results_pf, num_iter = nv_powerflow.solve(system)

# print node voltages
print("Powerflow converged in " + str(num_iter) + " iterations.\n")
print("Results: \n")
for node in results_pf.nodes:
    print('{}={}'.format(node.topology_node.uuid, node.voltage_pu))
    #print('{}={}'.format(node.topology_node.uuid, node.voltage))

