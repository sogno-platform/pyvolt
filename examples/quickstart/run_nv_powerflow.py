import logging
from pyvolt import network
from pyvolt import nv_powerflow
from pyvolt import results
import cimpy
import os
from villas.dataprocessing.readtools import read_timeseries_dpsim
from pathlib import Path

logging.basicConfig(filename='CIGRE.log', level=logging.INFO, filemode='w')

this_file_folder = Path(__file__).resolve().parent
xml_path = this_file_folder / "sample_data" / "CIGRE-MV-NoTap"
xml_files = [xml_path / "Rootnet_FULL_NE_06J16h_DI.xml",
             xml_path / "Rootnet_FULL_NE_06J16h_EQ.xml",
             xml_path / "Rootnet_FULL_NE_06J16h_SV.xml",
             xml_path / "Rootnet_FULL_NE_06J16h_TP.xml"]

xml_files_abs = []
for file in xml_files:
    xml_files_abs.append(os.path.abspath(file))

# read cim files and create new network.Systen object
res = cimpy.cim_import(xml_files_abs, "cgmes_v2_4_15")
system = network.System()
base_apparent_power = 25  # MW
system.load_cim_data(res['topology'], base_apparent_power)

# Execute power flow analysis
results_pf, num_iter = nv_powerflow.solve(system)

# print node voltages
print("Powerflow converged in " + str(num_iter) + " iterations.\n")
print("Results: \n")
for node in results_pf.nodes:
    print('{}={}'.format(node.topology_node.uuid, node.voltage_pu))
    #print('{}={}'.format(node.topology_node.uuid, node.voltage))

