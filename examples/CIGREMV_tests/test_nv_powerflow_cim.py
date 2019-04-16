import sys
import logging

sys.path.append("../../acs/state_estimation")
import network
import nv_powerflow_cim

sys.path.append("../../../cimpy")
import cimpy

logging.basicConfig(filename='CIGRE.log', level=logging.INFO, filemode='w')

cim_xml_path = r"C:\Users\Martin\Desktop\hiwi\git\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results"
cim_xml_files=[cim_xml_path + r"\Rootnet_FULL_NE_06J16h_DI.xml", 
		   cim_xml_path + r"\Rootnet_FULL_NE_06J16h_EQ.xml",
		   cim_xml_path + r"\Rootnet_FULL_NE_06J16h_SV.xml",
		   cim_xml_path + r"\Rootnet_FULL_NE_06J16h_TP.xml"]

#read cim files and create new network.Systen object
res=cimpy.cimread(cim_xml_files)
system = network.System()
system.load_cim_data(res, 20)

#print node voltages
for node in system.nodes:
    print('{}={}'.format(node.uuid, node.voltage))

print()

#print node powers
for node in system.nodes:
    print('{}={}'.format(node.uuid, node.power))

results, num_iter_cim = nv_powerflow_cim.solve(system)

print()
print("voltages:")
#print node voltages
for node in results.nodes:
    print('{}={}'.format(node.topology_node.uuid, node.voltage))

#results.print_voltages_polar()