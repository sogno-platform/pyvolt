import logging
from acs.state_estimation import network
from acs.state_estimation import nv_powerflow
from acs.state_estimation import results
import cimpy
import os
from villas.dataprocessing.readtools import read_timeseries_dpsim


#print(os.path.dirname(__file__))
#os.chdir(os.path.dirname(__file__))

logging.basicConfig(filename='CIGRE.log', level=logging.INFO, filemode='w')

xml_path = r".\sample_data"
xml_files = [xml_path + r"\Rootnet_FULL_NE_06J16h_DI.xml",
                 xml_path + r"\Rootnet_FULL_NE_06J16h_EQ.xml",
                 xml_path + r"\Rootnet_FULL_NE_06J16h_SV.xml",
                 xml_path + r"\Rootnet_FULL_NE_06J16h_TP.xml"]

xml_files_abs = []
for file in xml_files:
    xml_files_abs.append(os.path.abspath(file))

# read cim files and create new network.Systen object
res, _ = cimpy.cim_import(xml_files_abs, "cimgen_v2_4_15")
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
print("\n\n\n")

#load results of csv file:
csv_file = xml_path + r".\CIGRE-MV-NoTap.csv"
results_CIGRE = results.Results(system)
results_CIGRE.read_data_dpsim(csv_file)

# Show numerical comparison
print(results_CIGRE.get_voltages(pu=False)/1000-results_pf.get_voltages(pu=False))
