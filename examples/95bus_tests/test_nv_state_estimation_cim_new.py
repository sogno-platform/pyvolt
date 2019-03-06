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
	std_dev = np.absolute(PowerflowNode.voltage)*(V_unc/300)
	measurements_set.create_measurement(PowerflowNode.topology_node, nv_state_estimator_cim.ElemType.Node, nv_state_estimator_cim.MeasType.V_mag, np.absolute(PowerflowNode.voltage), std_dev)
for i in range(1,len(res.nodes)):
	PowerflowNode = res.get_node(i)
	std_dev = np.absolute(PowerflowNode.power.real*(Sinj_unc/300))
	measurements_set.create_measurement(PowerflowNode.topology_node, nv_state_estimator_cim.ElemType.Node, nv_state_estimator_cim.MeasType.Sinj_real, PowerflowNode.power.real, std_dev)
for i in range(1,len(res.nodes)):
	PowerflowNode = res.get_node(i)
	std_dev = np.absolute(PowerflowNode.power.imag*(Sinj_unc/300))
	measurements_set.create_measurement(PowerflowNode.topology_node, nv_state_estimator_cim.ElemType.Node, nv_state_estimator_cim.MeasType.Sinj_imag, PowerflowNode.power.imag, std_dev)
for i in [0,10,27,54,58]:
	std_dev = np.absolute(res.branches[i].power.real*(S_unc/300))
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch, nv_state_estimator_cim.MeasType.S1_real, res.branches[i].power.real, std_dev)
for i in [9,53]:
	std_dev = np.absolute(res.branches[i].power2.real*(S_unc/300))
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch, nv_state_estimator_cim.MeasType.S2_real, res.branches[i].power2.real, std_dev)
for i in [0,10,27,54,58]:	
	std_dev = np.absolute(res.branches[i].power.imag*(S_unc/300))
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch, nv_state_estimator_cim.MeasType.S1_imag, res.branches[i].power.imag, std_dev)
for i in [9,53]:
	std_dev = np.absolute(res.branches[i].power2.imag*(S_unc/300))
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch, nv_state_estimator_cim.MeasType.S2_imag, res.branches[i].power2.imag, std_dev)
for i in [12,35]:
	std_dev = np.absolute(res.branches[i].current)*(I_unc/300)
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch, nv_state_estimator_cim.MeasType.I_mag, np.absolute(res.branches[i].current), std_dev)
for i in [0,10,54]:
	PowerflowNode = res.get_node(i)
	std_dev = np.absolute(PowerflowNode.voltage)*(Pmu_mag_unc/300)
	measurements_set.create_measurement(PowerflowNode.topology_node, nv_state_estimator_cim.ElemType.Node, nv_state_estimator_cim.MeasType.Vpmu_mag, np.absolute(PowerflowNode.voltage), std_dev)
for i in [0,10,54]:
	PowerflowNode = res.get_node(i)
	std_dev = Pmu_phase_unc/300
	measurements_set.create_measurement(PowerflowNode.topology_node, nv_state_estimator_cim.ElemType.Node, nv_state_estimator_cim.MeasType.Vpmu_phase, np.angle(PowerflowNode.voltage), std_dev)
for i in [12,35]:
	std_dev = np.absolute(res.branches[i].current)*(Pmu_mag_unc/300)
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch, nv_state_estimator_cim.MeasType.Ipmu_mag, np.absolute(res.branches[i].current), std_dev)
for i in [12,35]:
	std_dev = Pmu_phase_unc/300
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch, nv_state_estimator_cim.MeasType.Ipmu_phase, np.angle(res.branches[i].current), std_dev)

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
zdatameas = Zdata_init(meas)

zdata, zdatameas = py_95bus_meas_data.Zdatatrue_creation(zdata, zdatameas, meas, branch, Vtrue, Itrue, Iinjtrue, S1true, S2true, Sinjtrue)
err_pu, zdatameas = py_95bus_meas_data.Zdatameas_creation(zdata, zdatameas)
measurements_set.meas_creation_test(err_pu)

#compare results:
#print("zdata.mval==measurements_set.measurements?: " + str((zdata.mval==measurements_set.getMeasValues()).all()))

Vest, Iest, Iinjest, S1est, S2est, Sinjest = nv_state_estimator.DsseCall(branch, node, zdatameas, system.Ymatrix, system.Adjacencies)
results = nv_state_estimator_cim.DsseCall(system, measurements_set) 
print("test==test_cim?: " + str(np.array_equal(np.around(Vest.complex,6), np.around(results.get_voltages(),6))))


"""
nodes_num = len(system.nodes)
Gmatrix = system.Ymatrix.real
Bmatrix = system.Ymatrix.imag
Yabs_matrix = np.absolute(system.Ymatrix)
Yphase_matrix = np.angle(system.Ymatrix)
Adj = system.Adjacencies
test = nv_state_estimator.DsseMixed(branch, node, zdatameas, system.Ymatrix, system.Adjacencies)
test_cim = nv_state_estimator_cim.DsseMixed(nodes_num, measurements_set, Gmatrix, Bmatrix, Yabs_matrix, Yphase_matrix, Adj)
print("test==test_cim?: " +  str(np.array_equal(np.around(test,6), np.around(test_cim,6))))
"""