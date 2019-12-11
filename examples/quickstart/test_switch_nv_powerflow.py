import logging
import numpy as np
from acs.state_estimation import network
from acs.state_estimation import nv_powerflow
from acs.state_estimation import nv_state_estimator
from acs.state_estimation import measurement
from acs.state_estimation import results
import cimpy
import os

logging.basicConfig(filename='CIGRE.log', level=logging.INFO, filemode='w')

xml_path = r".\sample_data\CIGRE-MV-NoTap-WithBreaker"
xml_files = [xml_path + r"\20191126T1535Z_YYY_EQ_.xml",
             xml_path + r"\20191126T1535Z_XX_YYY_SV_.xml",
             xml_path + r"\20191126T1535Z_XX_YYY_TP_.xml"]

xml_files_abs = []
for file in xml_files:
    xml_files_abs.append(os.path.abspath(file))

# read cim files and create new network.Systen object
res, _ = cimpy.cim_import(xml_files_abs, "cimgen_v2_4_15")
system = network.System()
base_apparent_power = 25  # MW
system.load_cim_data(res, base_apparent_power)

#open breaker
system.breakers[-1].open_breaker()
system.Ymatrix_calc()

# Execute power flow analysis
results_pf, num_iter = nv_powerflow.solve(system)

# print node voltages
print("Powerflow converged in " + str(num_iter) + " iterations.\n")
print("Results:")
for node in results_pf.nodes:
    print('{}={}'.format(node.topology_node.name, np.absolute(node.voltage)))
    #print('{}={}'.format(node.topology_node.name, np.absolute(node.voltage_pu)))
print("\n")

print("-------------------------------------------------------------")
print("\n")

#close breaker
system.breakers[-1].close_breaker()
system.Ymatrix_calc()

# Execute power flow analysis
results_pf, num_iter = nv_powerflow.solve(system)

# print node voltages
print("Powerflow converged in " + str(num_iter) + " iterations.\n")
print("Results:")
for node in results_pf.nodes:
    print('{}={}'.format(node.topology_node.name, np.absolute(node.voltage)))
    #print('{}={}'.format(node.topology_node.name, np.absolute(node.voltage_pu)))