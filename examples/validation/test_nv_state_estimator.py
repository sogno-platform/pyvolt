import os
import logging
import numpy as np
import matplotlib.pyplot as plt

import cimpy
from pyvolt import network
from pyvolt import nv_state_estimator
from pyvolt import measurement
from pyvolt import results


logging.basicConfig(filename='test_nv_state_estimator.log', level=logging.INFO, filemode='w')

this_file_folder = os.path.dirname(os.path.realpath(__file__))
xml_path = os.path.realpath(os.path.join(this_file_folder, "..", "sample_data", "CIGRE-MV-NoTap"))
xml_files = [os.path.join(xml_path, "CIGRE-MV-NoTap_DI.xml"),
             os.path.join(xml_path, "CIGRE-MV-NoTap_EQ.xml"),
             os.path.join(xml_path, "CIGRE-MV-NoTap_SV.xml"),
             os.path.join(xml_path, "CIGRE-MV-NoTap_TP.xml")]

# Read cim files and create new network.System object
res = cimpy.cim_import(xml_files, "cgmes_v2_4_15")
system = network.System()
base_apparent_power = 25  # MW
system.load_cim_data(res['topology'], base_apparent_power)

# Read input result file and store it in a results.Results object
loadflow_results_file = os.path.realpath(os.path.join(this_file_folder,
                                                      "..",
                                                      "sample_data",
                                                      "CIGRE-MV-NoTap",
                                                      "CIGRE-MV-NoTap.csv"))

powerflow_results = results.Results(system)
powerflow_results.read_data(loadflow_results_file)

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
measurements_set = measurement.MeasurementSet()
for node in powerflow_results.nodes:
    measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag,
                                        np.absolute(node.voltage_pu), Pmu_mag_unc)
for node in powerflow_results.nodes:
    measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase,
                                        np.angle(node.voltage_pu), Pmu_phase_unc)
measurements_set.meas_creation()

# Perform state estimation
state_estimation_results_ideal = nv_state_estimator.DsseCall(system, measurements_set)

# Perform numerical comparison
Vest_ideal = state_estimation_results_ideal.get_voltages(pu=False)
Vtrue = powerflow_results.get_voltages(pu=False)
print(Vest_ideal - Vtrue)

# --- State Estimation with Non-Ideal Measurements ---
""" Write here the percent uncertainties of the measurements"""
Pmu_mag_unc = 1

# Create measurements data structures
"""use all node voltages as measures"""
measurements_set = measurement.MeasurementSet()
for node in powerflow_results.nodes:
    measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag,
                                        np.absolute(node.voltage_pu), Pmu_mag_unc)
for node in powerflow_results.nodes:
    measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase,
                                        np.angle(node.voltage_pu), Pmu_phase_unc)
measurements_set.meas_creation()

# Perform state estimation
state_estimation_results_real = nv_state_estimator.DsseCall(system, measurements_set)

# Perform numerical comparison
Vest_real = state_estimation_results_real.get_voltages(pu=False)
print(Vest_real - Vtrue)

# Plot comparison
line_width = 6
fontsize = 26

plt.figure()
ax = plt.subplot()

nodes_uuid = [system.nodes[elem].uuid for elem in range(len(system.nodes))]

# Reorder and rescale results
idx_filter = np.argsort([int(uuid[1:]) for uuid in nodes_uuid])#[1:]
nodes_uuid_filtered = [nodes_uuid[idx] for idx in idx_filter]
Vtrue_filtered = [abs(Vtrue[idx] / 1000) for idx in idx_filter]
Vest_ideal_filtered = [abs(Vest_ideal[idx] / 1000) for idx in idx_filter]
Vest_real_filtered = [abs(Vest_real[idx] / 1000) for idx in idx_filter]

plt.plot(Vest_ideal_filtered, linewidth=line_width, linestyle='solid', label="state estimator (ideal measurements)")
plt.plot(Vtrue_filtered, linewidth=line_width, linestyle='dashed', label="DPsim load flow results")
plt.plot(Vest_real_filtered, linewidth=line_width, linestyle='dotted', label="state estimator (non-ideal measurements)")

plt.xticks(range(len(system.nodes)), fontsize=fontsize)
plt.yticks(fontsize=fontsize)
ax.set_xticklabels(nodes_uuid_filtered)
plt.ylabel("Node voltage [kV]", fontsize=fontsize)
plt.xlim([0, len(system.nodes) - 2])
plt.legend(fontsize=fontsize)
plt.show()
