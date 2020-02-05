import logging
import numpy as np
from acs.state_estimation import network
from acs.state_estimation import nv_powerflow
from acs.state_estimation import nv_state_estimator
from acs.state_estimation import measurement
from acs.state_estimation import results
import cimpy
import os

logging.basicConfig(filename='CIGRE.log', level=logging.INFO, filemode='w')

xml_path = r".\sample_data\CIGRE-MV-NoTap-WithBreaker"
xml_files = [xml_path + r"\20191126T1535Z_YYY_EQ_.xml",
             xml_path + r"\20191126T1535Z_XX_YYY_SV_.xml",
             xml_path + r"\20191126T1535Z_XX_YYY_TP_.xml"]

xml_files_abs = []
for file in xml_files:
    xml_files_abs.append(os.path.abspath(file))

# read cim files and create new network.Systen object
res, _ = cimpy.cim_import(xml_files_abs, "cgmes_v2_4_15")
system = network.System()
base_apparent_power = 25  # MW
system.load_cim_data(res, base_apparent_power)

#open breaker
system.breakers[-1].open_breaker()
system.Ymatrix_calc()

# Execute power flow analysis
results_pf, num_iter = nv_powerflow.solve(system)

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
    print('{}={}'.format(node.topology_node.name, node.voltage))