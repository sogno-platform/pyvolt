import numpy as np
import math
import sys
import copy 

sys.path.append("../../acs/state_estimation")
import py_95bus_network_data
import py_95bus_meas_data
import nv_state_estimator
import nv_powerflow

import network
import nv_powerflow_cim
import nv_state_estimator_cim

class PerUnit:
    def __init__(self, S, V):
        self.S = S
        self.V = V
        self.I = S/V
        self.Z = S/(V**2)
        
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
		
		
""" Insert here per unit values of the grid for power and voltage """
S = 100*(10**6)
V = (11*(10**3))/math.sqrt(3)
slackV = 1.02

Base = PerUnit(S,V)
branch, node = py_95bus_network_data.Network_95_nodes(Base, slackV)
Vtrue, Itrue, Iinjtrue, S1true, S2true, Sinjtrue, num_iter = nv_powerflow.solve(branch, node)

system = network.load_python_data(node, branch, node.type)
res, num_iter = nv_powerflow_cim.solve(system)

""" Write here the indexes of the nodes/branches where the measurements are"""
V_idx = np.array([1,11,55])
I_idx = np.array([13,36])
Sinj_idx = np.linspace(2,node.num,node.num-1)
S1_idx = np.array([1,11,28,55,59])
S2_idx = np.array([10,54])
Ipmu_idx = np.array([13,36])
Vpmu_idx = np.array([1,11,55])

""" Write here the percent uncertainties of the measurements""" 
V_unc = 1
I_unc = 2
Sinj_unc = 5
S_unc = 2
Pmu_mag_unc = 0.7
Pmu_phase_unc = 0.7

""" create measurements set (for nv_state_estimator_cim.py)"""
measurements_set = nv_state_estimator_cim.Measurents_set()
for i in [0,10,54]:
	PowerflowNode = res.get_node(i)
	measurements_set.create_measurement(PowerflowNode.topology_node, nv_state_estimator_cim.ElemType.Node, nv_state_estimator_cim.MeasType.V, PowerflowNode.voltage, V_unc)
	measurements_set.create_measurement(PowerflowNode.topology_node, nv_state_estimator_cim.ElemType.Node, nv_state_estimator_cim.MeasType.Vpmu, PowerflowNode.voltage, Pmu_mag_unc)
for i in [12,35]:
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch,  nv_state_estimator_cim.MeasType.I, res.branches[i].current, I_unc)
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch,  nv_state_estimator_cim.MeasType.Ipmu, res.branches[i].current, Pmu_mag_unc)
for i in range(1,len(res.nodes)):
	PowerflowNode = res.get_node(i)
	measurements_set.create_measurement(PowerflowNode.topology_node, nv_state_estimator_cim.ElemType.Node, nv_state_estimator_cim.MeasType.Sinj, PowerflowNode.power, Sinj_unc)
for i in [0,10,27,54,58]:
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch,  nv_state_estimator_cim.MeasType.S1, res.branches[i].power, I_unc)
for i in [9,53]:
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch,  nv_state_estimator_cim.MeasType.S2, res.branches[i].power2, I_unc)

print("########################   TEST   ####################")
Meas_mag = np.array([])
Meas_phase = np.array([])
for meas in measurements_set.measurements_set:
	if meas.meas_type == nv_state_estimator_cim.MeasType.Vpmu:
		Meas_mag = np.append(Meas_mag, np.absolute(meas.meas_value))
		Meas_phase = np.append(Meas_phase, np.angle(meas.meas_value))
		#Meas_mag = np.append(Meas_mag, meas.meas_value.real)
		#Meas_phase = np.append(Meas_phase, meas.meas_value.imag)
print(Meas_mag)
print(Meas_phase)

	
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
zdatameas =  Zdata_init(meas)

zdata, zdatameas = py_95bus_meas_data.Zdatatrue_creation(zdata, zdatameas, meas, branch, Vtrue, Itrue, Iinjtrue, S1true, S2true, Sinjtrue)