# The state estimation is performed based on the results using the nv_powerflow implementation.

import os
import logging
from pathlib import Path
import numpy as np
import pandas as pd
import json

import cimpy
from pyvolt import network
from pyvolt import nv_powerflow
from pyvolt import nv_state_estimator
from pyvolt import measurement
from pyvolt import results

def create_measurements(system):
    """Execute power flow analysis to generate voltage measurement file"""
    results_pf, num_iter_cim = nv_powerflow.solve(system)

    voltages = {}
    # write to csv file
    #for node in results_pf.nodes:
    #    voltages[node.topology_node.uuid + '.V_mag'] = [np.absolute(node.voltage_pu)]
    #
    #voltages_pd = pd.DataFrame.from_dict(voltages)
    #voltages_pd.to_json('voltage_mag_pu.csv')

    # write to json file
    for node in results_pf.nodes:
        voltages[node.topology_node.uuid] = {}
        voltages[node.topology_node.uuid]['Vpmu_mag'] = np.absolute(node.voltage_pu)
        voltages[node.topology_node.uuid]['Vpmu_phase'] = np.angle(node.voltage_pu)

    with open("voltage_mag_pu.json", "w") as outfile:
        json.dump(voltages, outfile)

    return results_pf

logging.basicConfig(filename='run_nv_state_estimator.log', level=logging.INFO, filemode='w')

this_file_folder = os.path.dirname(os.path.realpath(__file__))
xml_path = os.path.realpath(os.path.join(this_file_folder, "..", "sample_data", "CIGRE-MV-NoTap"))
xml_files = [os.path.join(xml_path, "CIGRE-MV-NoTap_DI.xml"),
             os.path.join(xml_path, "CIGRE-MV-NoTap_EQ.xml"),
             os.path.join(xml_path, "CIGRE-MV-NoTap_SV.xml"),
             os.path.join(xml_path, "CIGRE-MV-NoTap_TP.xml")]

# Read cim files and create new network.System object
res = cimpy.cim_import(xml_files, "cgmes_v2_4_15")
system = network.System()
base_apparent_power = 25  # in MW to compute pu values
system.load_cim_data(res['topology'], base_apparent_power)

# Read measurement config file
with open('run_nv_state_estimator_meas_config.json') as meas_config_file:
    meas_config = json.load(meas_config_file)
print(meas_config)

# Create artificial measurements from powerflow (optional if real measurements available)
results_pf = create_measurements(system)

# read measurements from file
with open('voltage_mag_pu.json') as measurements_file:
    meas_data = json.load(measurements_file)
print(meas_data)

# Create measurements data structures
# use all node voltages as measures
meas_set = measurement.MeasurementSet(system)
meas_set.read_measurements_from_dict(meas_data, meas_config)
# Add noise to the artificial measurements
meas_set.meas_creation()

# Perform state estimation
se_results = nv_state_estimator.DsseCall(system, meas_set)

# Print node voltages
print("state_estimation_results.voltages: ")
for node in se_results.nodes:
    print('{}={}'.format(node.topology_node.uuid, node.voltage))

Vest_ideal = se_results.get_voltages(pu=False)
Vtrue = results_pf.get_voltages(pu=False)
print(Vest_ideal - Vtrue)
