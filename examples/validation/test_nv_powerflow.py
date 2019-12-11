import os
import logging
from acs.state_estimation import network
from acs.state_estimation import nv_powerflow
from acs.state_estimation import results
import cimpy

logging.basicConfig(filename='CIGRE.log', level=logging.INFO, filemode='w')

xml_path = r"..\quickstart\sample_data\CIGRE-MV-NoTap"
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
results_pf, num_iter_cim = nv_powerflow.solve(system)

# print node voltages
print("results_pf.voltages: ")
for node in results_pf.nodes:
    print('{}={}'.format(node.topology_node.uuid, node.voltage))
print("\n")

# Show numerical comparison 
loadflow_results_path = r"..\quickstart\sample_data\CIGRE-MV-NoTap"

loadflow_results_file = loadflow_results_path + r"\CIGRE-MV-NoTap.csv"
results_dpsim = results.Results(system)
results_dpsim.read_data_dpsim(file_name=loadflow_results_file)

print("numerical comparison : results_pf.voltages - results_dpsim.voltages ")
for pf_node in results_pf.nodes:
    dpsim_node = results_dpsim.get_node(uuid=pf_node.topology_node.uuid)
    diff = pf_node.voltage - dpsim_node.voltage / 1000
    print('pf_node.{}-dpsim_node.{} = {}'.format(pf_node.topology_node.uuid, dpsim_node.topology_node.uuid, diff))
