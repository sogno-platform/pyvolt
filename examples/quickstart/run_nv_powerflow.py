import logging
from pyvolt import network
from pyvolt import nv_powerflow
from pyvolt import results
import cimpy
import os
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
voltages = []
for node in results_pf.nodes:
    print('{}={}'.format(node.topology_node.uuid, node.voltage_pu))
    #print('{}={}'.format(node.topology_node.uuid, node.voltage))
    voltages.append(node.voltage_pu)

voltages_ref = [(1-7.970485900477431e-27j), (0.9521818868802214-0.11692768153747995j), (0.9642955926931457-0.09862127081290231j), (0.8796973782245792-0.15318580971335868j), (0.8758799767843979-0.1554670528853566j), (0.8774339089327876-0.1545390713879984j), (0.8704521134131005-0.1574589214738466j), (0.8719578342107204-0.15661905676450852j), (0.8731024990161087-0.1561497673751958j), (0.8727417602545003-0.15623359013901236j), (0.8740410220986286-0.1565630049874985j), (0.9563016474701451-0.09917826765833906j), (0.9592141716734833-0.09896267637101246j), (0.8702137462858025-0.15760036065945185j), (0.9239489705253996-0.13105032262255972j)]
assert (voltages == voltages_ref), "Results do not match reference results." 