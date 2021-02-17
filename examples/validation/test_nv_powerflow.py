import os
import logging
from pyvolt import network
from pyvolt import nv_powerflow
from pyvolt import results
import cimpy


logging.basicConfig(filename='test_nv_powerflow.log', level=logging.INFO, filemode='w')

this_file_folder = os.path.dirname(os.path.realpath(__file__))
xml_path = os.path.realpath(os.path.join(this_file_folder, "..", "sample_data", "CIGRE-MV-NoTap"))
xml_files = [os.path.join(xml_path, "Rootnet_FULL_NE_06J16h_DI.xml"),
             os.path.join(xml_path, "Rootnet_FULL_NE_06J16h_EQ.xml"),
             os.path.join(xml_path, "Rootnet_FULL_NE_06J16h_SV.xml"),
             os.path.join(xml_path, "Rootnet_FULL_NE_06J16h_TP.xml")]

# Read cim files and create new network.System object
res = cimpy.cim_import(xml_files, "cgmes_v2_4_15")
system = network.System()
base_apparent_power = 25  # MW
system.load_cim_data(res['topology'], base_apparent_power)

# Execute power flow analysis
results_pf, num_iter_cim = nv_powerflow.solve(system)

# Print node voltages
print("results_pf.voltages: ")
for node in results_pf.nodes:
    print('{}={}'.format(node.topology_node.uuid, node.voltage))
print("\n")

# Perform numerical comparison
loadflow_results_file = os.path.realpath(os.path.join(this_file_folder,
                                                      "..",
                                                      "sample_data",
                                                      "CIGRE-MV-NoTap",
                                                      "CIGRE-MV-NoTap.csv"))

results_dpsim = results.Results(system)
results_dpsim.read_data(loadflow_results_file)

print("numerical comparison : results_pf.voltages - results_dpsim.voltages ")
for pf_node in results_pf.nodes:
    dpsim_node = results_dpsim.get_node(uuid=pf_node.topology_node.uuid)
    diff = pf_node.voltage - dpsim_node.voltage / 1000
    print('pf_node.{}-dpsim_node.{} = {}'.format(pf_node.topology_node.uuid, dpsim_node.topology_node.uuid, diff))
