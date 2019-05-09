import logging
from  acs.state_estimation import network
from  acs.state_estimation import nv_powerflow
from  acs.state_estimation import results
import cimpy

logging.basicConfig(filename='CIGRE.log', level=logging.INFO, filemode='w')

cim_xml_path = r".\examples\quickstart\sample_data"
cim_xml_files=[cim_xml_path + r"\Rootnet_FULL_NE_06J16h_DI.xml", 
		   cim_xml_path + r"\Rootnet_FULL_NE_06J16h_EQ.xml",
		   cim_xml_path + r"\Rootnet_FULL_NE_06J16h_SV.xml",
		   cim_xml_path + r"\Rootnet_FULL_NE_06J16h_TP.xml"]

#read cim files and create new network.Systen object
res=cimpy.cimread(cim_xml_files)
system = network.System()
base_apparent_power = 25    #MW
system.load_cim_data(res, base_apparent_power)

#Execute power flow analysis
results_pf, num_iter_cim = nv_powerflow.solve(system)

#print node voltages
print("results_pf.voltages: ")
for node in results_pf.nodes:
    print('{}={}'.format(node.topology_node.uuid, node.voltage))
print("\n\n\n")