import logging
import numpy as np
import matplotlib.pyplot as plt

import cimpy
from  acs.state_estimation import network
from  acs.state_estimation import nv_state_estimator
from  acs.state_estimation import measurement
from  acs.state_estimation import results

logging.basicConfig(filename='CIGRE.log', level=logging.INFO, filemode='w')

cim_xml_path = r"..\quickstart\sample_data"
cim_xml_files=[cim_xml_path + r"\Rootnet_FULL_NE_06J16h_DI.xml", 
		   cim_xml_path + r"\Rootnet_FULL_NE_06J16h_EQ.xml",
		   cim_xml_path + r"\Rootnet_FULL_NE_06J16h_SV.xml",
		   cim_xml_path + r"\Rootnet_FULL_NE_06J16h_TP.xml"]

#read cim files and create new network.Systen object
res=cimpy.cimread(cim_xml_files)
system = network.System()
base_apparent_power = 25    #MW
system.load_cim_data(res, base_apparent_power)

#read Input-Ergebnisdatei and store it in a results.Results object
loadflow_results_path = r"..\quickstart\sample_data"
loadflow_results_file = loadflow_results_path + r"\CIGRE-MV-NoTap.csv" 
powerflow_results = results.Results(system)
powerflow_results.read_data_dpsim(loadflow_results_file)

# --- State Estimation with Ideal Measurements ---
""" Write here the percent uncertainties of the measurements""" 
V_unc = 0
I_unc = 0
Sinj_unc = 0
S_unc = 0
Pmu_mag_unc = 0
Pmu_phase_unc = 0

# Create measurements object
"""use all node voltages as measures"""
measurements_set = measurement.Measurents_set()
for node in powerflow_results.nodes:
	measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag, np.absolute(node.voltage_pu), Pmu_mag_unc)
for node in powerflow_results.nodes:
	measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase, np.angle(node.voltage_pu), Pmu_phase_unc)
measurements_set.meas_creation()

# Perform state estimation
state_estimation_results_ideal = nv_state_estimator.DsseCall(system, measurements_set) 

# Show numerical comparison
Vest_ideal = state_estimation_results_ideal.get_voltages(pu=False)
Vtrue = powerflow_results.get_voltages(pu=False)
print(Vest_ideal - Vtrue)

# --- State Estimation with Non-Ideal Measurements ---
""" Write here the percent uncertainties of the measurements""" 
Pmu_mag_unc = 1

# Create measurements data structures
"""use all node voltages as measures"""
measurements_set = measurement.Measurents_set()
for node in powerflow_results.nodes:
	measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag, np.absolute(node.voltage_pu), Pmu_mag_unc)
for node in powerflow_results.nodes:
	measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase, np.angle(node.voltage_pu), Pmu_phase_unc)
measurements_set.meas_creation()

# Perform state estimation
state_estimation_results_real = nv_state_estimator.DsseCall(system, measurements_set) 

# Show numerical comparison
Vest_real = state_estimation_results_real.get_voltages(pu=False)
print(Vest_real - Vtrue)

# Plot comparison
line_width = 6
fontsize = 26

plt.figure()
ax = plt.subplot()

nodes_uuid = [system.nodes[elem].uuid for elem in range(len(system.nodes))]

# Reorder and rescale results
idx_filter = np.argsort([int(uuid[1:]) for uuid in nodes_uuid])[1:]
nodes_uuid_filtered = [nodes_uuid[idx] for idx in idx_filter]
Vtrue_filtered = [abs(Vtrue[idx]/1000) for idx in idx_filter]
Vest_ideal_filtered = [abs(Vest_ideal[idx]/1000) for idx in idx_filter]
Vest_real_filtered = [abs(Vest_real[idx]/1000) for idx in idx_filter]


plt.plot(Vest_ideal_filtered, linewidth=line_width, linestyle='-', label="state estimator (ideal measurements)")
plt.plot(Vtrue_filtered, linewidth=line_width, linestyle=':', label="DPsim load flow results")
plt.plot(Vest_real_filtered, linewidth=line_width, linestyle='-', label="state estimator (non-ideal measurements)")

plt.xticks(range(len(system.nodes)), fontsize=fontsize)
plt.yticks(fontsize=fontsize)
ax.set_xticklabels(nodes_uuid_filtered)
plt.ylabel("Node voltage [kV]", fontsize=fontsize)
plt.xlim([0,len(system.nodes)-2])
plt.legend(fontsize=fontsize)
plt.show()