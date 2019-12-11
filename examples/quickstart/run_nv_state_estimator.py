# This example will provide a test case running the state estimator.
# The state estimation is performed based on the results using the nv_powerflow implementation

import logging
import numpy as np
from acs.state_estimation import network
from acs.state_estimation import nv_powerflow
from acs.state_estimation import nv_state_estimator
from acs.state_estimation import measurement
from acs.state_estimation import results
import cimpy
import os

#os.chdir(os.path.dirname(__file__))

logging.basicConfig(filename='CIGRE.log', level=logging.INFO, filemode='w')

xml_path = r".\sample_data\CIGRE-MV-NoTap"
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

# --- State Estimation ---
""" Write here the percent uncertainties of the measurements"""
V_unc = 0
I_unc = 0
Sinj_unc = 0
S_unc = 0
Pmu_mag_unc = 1
Pmu_phase_unc = 0

# Create measurements data structures
"""use all node voltages as measures"""
measurements_set = measurement.MeasurementSet()
for node in results_pf.nodes:
    measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_mag,
                                        np.absolute(node.voltage_pu), Pmu_mag_unc)
for node in results_pf.nodes:
    measurements_set.create_measurement(node.topology_node, measurement.ElemType.Node, measurement.MeasType.Vpmu_phase,
                                        np.angle(node.voltage_pu), Pmu_phase_unc)
measurements_set.meas_creation()

# Perform state estimation
state_estimation_results = nv_state_estimator.DsseCall(system, measurements_set)

# print node voltages
print("state_estimation_results.voltages: ")
for node in state_estimation_results.nodes:
    print('{}={}'.format(node.topology_node.uuid, node.voltage))
