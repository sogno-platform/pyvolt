import logging
from pathlib import Path
from pyvolt import network
from pyvolt import nv_powerflow
import numpy
import cimpy
import os


logging.basicConfig(filename='run_nv_powerflow.log', level=logging.INFO, filemode='w')

#python starts as module in subdirectory, 2 folders up to set the new path
this_file_folder =  Path(__file__).parents[2]
p = str(this_file_folder)+"/examples/sample_data/CIGRE-MV-NoTap"
xml_path = Path(p)


xml_files = [os.path.join(xml_path, "CIGRE-MV-NoTap_DI.xml"),
             os.path.join(xml_path, "CIGRE-MV-NoTap_EQ.xml"),
             os.path.join(xml_path, "CIGRE-MV-NoTap_SV.xml"),
             os.path.join(xml_path, "CIGRE-MV-NoTap_TP.xml")]

# Read cim files and create new network.System object
res = cimpy.cim_import(xml_files, "cgmes_v2_4_15")
system = network.System()
base_apparent_power = 25  # MW
system.load_cim_data(res['topology'], base_apparent_power)

# Execute power flow analysis
results_pf, num_iter = nv_powerflow.solve(system)

# Print node voltages
print("Powerflow converged in " + str(num_iter) + " iterations.\n")
print("Results: \n")
voltages = []
for node in results_pf.nodes:
    print('{}={}'.format(node.topology_node.uuid, node.voltage_pu))
    #print('{}={}'.format(node.topology_node.uuid, node.voltage))
    voltages.append(node.voltage_pu)

voltages_ref = [(1-7.970485900477431e-27j), (0.9521818868802214-0.11692768153747995j),
                (0.9642955926931457-0.09862127081290231j), (0.8796973782245792-0.15318580971335868j),
                (0.8758799767843979-0.1554670528853566j), (0.8774339089327876-0.1545390713879984j),
                (0.8704521134131005-0.1574589214738466j), (0.8719578342107204-0.15661905676450852j),
                (0.8731024990161087-0.1561497673751958j), (0.8727417602545003-0.15623359013901236j),
                (0.8740410220986286-0.1565630049874985j), (0.9563016474701451-0.09917826765833906j),
                (0.9592141716734833-0.09896267637101246j), (0.8702137462858025-0.15760036065945185j),
                (0.9239489705253996-0.13105032262255972j)]
epsilon = 1e-4
numpy.testing.assert_array_almost_equal(voltages, voltages_ref), "Results do not match reference results."
