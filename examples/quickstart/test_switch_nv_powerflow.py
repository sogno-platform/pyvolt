import logging
import numpy as np
from pyvolt import network
from pyvolt import nv_powerflow
import cimpy
import os


logging.basicConfig(filename='test_switch_nv_powerflow.log', level=logging.INFO, filemode='w')

this_file_folder = os.path.dirname(os.path.realpath(__file__))
xml_path = os.path.realpath(os.path.join(this_file_folder, "..", "sample_data", "CIGRE-MV-NoTap-WithBreaker"))
xml_files = [os.path.join(xml_path, "20191126T1535Z_YYY_EQ_.xml"),
             os.path.join(xml_path, "20191126T1535Z_XX_YYY_SV_.xml"),
             os.path.join(xml_path, "20191126T1535Z_XX_YYY_TP_.xml")]

# Read cim files and create new network.System object
res = cimpy.cim_import(xml_files, "cgmes_v2_4_15")
system = network.System()
base_apparent_power = 25  # MW
system.load_cim_data(res['topology'], base_apparent_power)

# Open breaker
system.breakers[-1].open_breaker()
system.Ymatrix_calc()

# Execute power flow analysis
results_pf, num_iter = nv_powerflow.solve(system)

# Print node voltages
print("Powerflow converged in " + str(num_iter) + " iterations.\n")
print("Results:")
for node in results_pf.nodes:
    print('{}={}'.format(node.topology_node.name, np.absolute(node.voltage)))
    #print('{}={}'.format(node.topology_node.name, np.absolute(node.voltage_pu)))
print("\n")

print("-------------------------------------------------------------")
print("\n")

# Close breaker
system.breakers[-1].close_breaker()
system.Ymatrix_calc()

# Execute power flow analysis
results_pf, num_iter = nv_powerflow.solve(system)

# Print node voltages
print("Powerflow converged in " + str(num_iter) + " iterations.\n")
print("Results:")
for node in results_pf.nodes:
    print('{}={}'.format(node.topology_node.name, np.absolute(node.voltage)))
    #print('{}={}'.format(node.topology_node.name, np.absolute(node.voltage_pu)))
