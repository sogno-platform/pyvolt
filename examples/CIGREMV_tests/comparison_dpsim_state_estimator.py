import numpy as np
import math
import sys
import copy
import logging
import matplotlib.pyplot as plt

sys.path.append("../../../cimpy")
import cimpy

sys.path.append("../../acs/state_estimation")
import network
import nv_powerflow_cim
from measurement_generator import *
sys.path.append("../../../dataprocessing")
from villas.dataprocessing.readtools import *
import nv_state_estimator_cim

logging.basicConfig(filename='CIGRE.log', level=logging.INFO)

cim_xml_path = r"D:\git\data\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results"
cim_xml_files=[cim_xml_path + r"\Rootnet_FULL_NE_06J16h_DI.xml", 
		   cim_xml_path + r"\Rootnet_FULL_NE_06J16h_EQ.xml",
		   cim_xml_path + r"\Rootnet_FULL_NE_06J16h_SV.xml",
		   cim_xml_path + r"\Rootnet_FULL_NE_06J16h_TP.xml"]

loadflow_results_path = r"D:\git\data\reference-results\DPsim\StaticPhasor"
loadflow_results_file = loadflow_results_path + r"\CIGRE-MV-NoTap.csv" 

res=cimpy.cimread(cim_xml_files)
cimpy.setNodes(res)
cimpy.setPowerTransformerEnd(res)
system = network.System()
system.load_cim_data(res)
Ymatr, Adj = network.Ymatrix_calc(system)

#read Input-Ergebnisdatei
ts_dpsim = read_timeseries_csv(loadflow_results_file)

#order readed data according to system.nodes
Vtrue=np.zeros(len(ts_dpsim), dtype=np.complex_)
for elem in range(len(system.nodes)): 
	Vtrue[elem] = ts_dpsim[system.nodes[elem].uuid].values[0]

""" Write here the indexes of the nodes/branches where the measurements are"""
V_idx = np.array([])
I_idx = np.array([])
Sinj_idx = np.array([])
S1_idx = np.array([])
S2_idx = np.array([])
Ipmu_idx = np.array([])
Vpmu_idx= np.linspace(1, len(Vtrue), len(Vtrue))

# --- State Estimation with Ideal Measurements ---
""" Write here the percent uncertainties of the measurements""" 
V_unc = 0
I_unc = 0
Sinj_unc = 0
S_unc = 0
Pmu_mag_unc = 0
Pmu_phase_unc = 0

# Create measurements data structures
V = Measurements(V_idx,V_unc)
I = Measurements(I_idx,I_unc)
Sinj = Measurements(Sinj_idx,Sinj_unc)
S1 = Measurements(S1_idx,S_unc)
S2 = Measurements(S2_idx,S_unc)
Ipmu_mag = Measurements(Ipmu_idx,Pmu_mag_unc)
Ipmu_phase = Measurements(Ipmu_idx,Pmu_phase_unc)
Vpmu_mag = Measurements(Vpmu_idx,Pmu_mag_unc)
Vpmu_phase = Measurements(Vpmu_idx,Pmu_phase_unc)

meas_ideal = Measurement_set(V, I, Sinj, S1, S2, Vpmu_mag, Vpmu_phase, Ipmu_mag, Ipmu_phase)
zdata_ideal = Zdata_init(meas_ideal)
zdata_ideal.mtype = np.ones(meas_ideal.V.num)
zdatameas_ideal =  Zdata_init(meas_ideal)

zdata_ideal, zdatameas_ideal = Zdatatrue_creation(zdata_ideal, zdatameas_ideal, meas_ideal, system, Vtrue, np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))
zdatameas_ideal = Zdatameas_creation(zdata_ideal, zdatameas_ideal)

# Perform state estimation
Vest_ideal, Iest, Iinjest, S1est, S2est, Sinjest = nv_state_estimator_cim.DsseCall(system, zdatameas_ideal, Ymatr, Adj)

# Show numerical comparison
print(Vest_ideal.complex - Vtrue)

# --- State Estimation with Non-Ideal Measurements ---
""" Write here the percent uncertainties of the measurements""" 
Pmu_mag_unc = 1

# Create measurements data structures
Vpmu_mag = Measurements(Vpmu_idx,Pmu_mag_unc)

meas_real = Measurement_set(V, I, Sinj, S1, S2, Vpmu_mag, Vpmu_phase, Ipmu_mag, Ipmu_phase)
zdata_real = Zdata_init(meas_real)
zdata_real.mtype = np.ones(meas_real.V.num)
zdatameas_real =  Zdata_init(meas_real)

zdata_real, zdatameas_real = Zdatatrue_creation(zdata_real, zdatameas_real, meas_real, system, Vtrue, np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))
zdatameas_real = Zdatameas_creation(zdata_real, zdatameas_real)

# Perform state estimation
Vest_real, Iest, Iinjest, S1est, S2est, Sinjest = nv_state_estimator_cim.DsseCall(system, zdatameas_real, Ymatr, Adj)

# Show numerical comparison
print(Vest_real.complex - Vtrue)

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
Vest_ideal_filtered = [abs(Vest_ideal.complex[idx]/1000) for idx in idx_filter]
Vest_real_filtered = [abs(Vest_real.complex[idx]/1000) for idx in idx_filter]


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