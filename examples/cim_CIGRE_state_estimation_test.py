import numpy as np
import math
import sys
import copy
import logging

sys.path.append("../../cimpy")
import cimpy

sys.path.append("../acs/state_estimation")
import network
import nv_powerflow_cim
import cim_py_95bus_meas_data
import nv_state_estimator_cim

logging.basicConfig(filename='CIGRE.log', level=logging.INFO)

class Measurements:
	def __init__(self, index, unc):
		self.index = index.astype(int)
		self.unc = unc
		self.num = len(index)
		
class Measurement_set:
	def __init__(self, V, I, Sinj, S1, S2, Vpmu_mag, Vpmu_phase, Ipmu_mag, Ipmu_phase):
		self.V = V
		self.Sinj = Sinj
		self.S1 = S1
		self.S2 = S2
		self.I = I
		self.Vpmu_mag = Vpmu_mag
		self.Vpmu_phase = Vpmu_phase 
		self.Ipmu_mag = Ipmu_mag
		self.Ipmu_phase = Ipmu_phase
		
class Zdata_init:
	def __init__(self, meas):
		nmeas = meas.V.num + meas.I.num + 2*meas.Sinj.num + 2*meas.S1.num + 2*meas.S2.num + meas.Ipmu_mag.num + meas.Ipmu_phase.num + meas.Vpmu_mag.num + meas.Vpmu_phase.num
		self.mtype = np.zeros(nmeas)
		self.mval = np.zeros(nmeas)
		self.mbranch = np.zeros(nmeas)
		self.mfrom = np.zeros(nmeas)
		self.mto = np.zeros(nmeas)
		self.mstddev = np.zeros(nmeas)

xml_files=[r"..\..\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_DI.xml", 
		   r"..\..\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_EQ.xml",
		   r"..\..\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_SV.xml",
		   r"..\..\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_TP.xml"]

res=cimpy.cimread(xml_files)
cimpy.setNodes(res)
cimpy.setPowerTransformerEnd(res)
system = network.System()
system.load_cim_data(res)

def read_data(file_name):
	file = []
	with open(file_name, "r") as filestream:
		for line in filestream:
			file.append([x.strip() for x in line.split(',')])
	
	#remove "time"
	file[0].pop(0) 
	file[1].pop(0)
	
	results={}
	for elem in range(len(file[0])):
		node, type = file[0][elem].split('.')
		if not node in results:
			results[node] = 0.0 + 0.0j
		if type=="real":
			results[node] += float(file[1][elem])
		elif type=="imag":
			results[node] += complex(0,float(file[1][elem]))
	
	return results
	
file = r"..\..\reference-results\DPsim\StaticPhasor\CIGRE-MV-NoTap.csv" 
measurements = read_data(file)

#order readed data according to system.nodes
Vtrue=np.zeros(len(measurements), dtype=np.complex_)
for elem in range(len(system.nodes)): 
	Vtrue[elem] = measurements[system.nodes[elem].uuid]

""" Write here the indexes of the nodes/branches where the measurements are"""
V_idx = np.linspace(1, len(measurements), len(measurements))
I_idx = np.array([])
Sinj_idx = np.array([])
S1_idx = np.array([])
S2_idx = np.array([])
Ipmu_idx = np.array([])
#Vpmu_idx = np.array([])
Vpmu_idx=copy.deepcopy(V_idx)

""" Write here the percent uncertainties of the measurements""" 
V_unc = 0
I_unc = 0
Sinj_unc = 0
S_unc = 0
Pmu_mag_unc = 0
Pmu_phase_unc = 0

V = Measurements(V_idx,V_unc)
I = Measurements(I_idx,I_unc)
Sinj = Measurements(Sinj_idx,Sinj_unc)
S1 = Measurements(S1_idx,S_unc)
S2 = Measurements(S2_idx,S_unc)
Ipmu_mag = Measurements(Ipmu_idx,Pmu_mag_unc)
Ipmu_phase = Measurements(Ipmu_idx,Pmu_phase_unc)
Vpmu_mag = Measurements(Vpmu_idx,Pmu_mag_unc)
Vpmu_phase = Measurements(Vpmu_idx,Pmu_phase_unc)
meas = Measurement_set(V, I, Sinj, S1, S2, Vpmu_mag, Vpmu_phase, Ipmu_mag, Ipmu_phase)
zdata = Zdata_init(meas)
zdata.mtype = np.ones(meas.V.num)
zdatameas =  Zdata_init(meas)

zdata, zdatameas = cim_py_95bus_meas_data.Zdatatrue_creation(zdata, zdatameas, meas, system, Vtrue, np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))

MC_trials = 1
iter_counter = np.zeros(MC_trials)
Ymatr_cim, Adj_cim = network.Ymatrix_calc(system)

for mciter in range(0,MC_trials):	  
	zdatameas = cim_py_95bus_meas_data.Zdatameas_creation(zdata, zdatameas)
	Vest_cim, Iest_cim, Iinjest_cim, S1est_cim, S2est_cim, Sinjest_cim = nv_state_estimator_cim.DsseCall(system, zdatameas, Ymatr_cim, Adj_cim)

print(Vest_cim.complex - Vtrue)