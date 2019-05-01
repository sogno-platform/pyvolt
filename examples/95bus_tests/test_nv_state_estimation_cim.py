import numpy as np
import math
import py_95bus_meas_data
import py_95bus_network_data
from acs.state_estimation import nv_state_estimator
from acs.state_estimation import nv_powerflow
from acs.state_estimation import network
from acs.state_estimation import nv_powerflow_cim
from acs.state_estimation import nv_state_estimator_cim
from acs.state_estimation import measurement

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

def load_python_data(nodes, branches, type):
	"""
	to create a new network.System object from the objects 
	py_95bus_network_data.Node and py_95bus_network_data.Branch
	"""
	system = network.System()
	
	for node_idx in range(0, nodes.num):
		node = network.Node()
		node.index = node_idx
		node.type = network.BusType[type[node_idx]]
		if node.type == network.BusType.slack:
			node.voltage_pu = nodes.P2[node_idx]*np.cos(nodes.Q2[node_idx]) + 1j * nodes.P2[node_idx]*np.sin(nodes.Q2[node_idx])
		elif node.type == network.BusType.PQ:
			node.power_pu = complex(nodes.P2[node_idx], nodes.Q2[node_idx])
		elif node.type == network.BusType.PV:
			pass
		system.nodes.append(node)

	for branch_idx in range(0, branches.num):
		branch = network.Branch()
		branch.r_pu = branches.R[branch_idx]
		branch.x_pu = branches.X[branch_idx]
		branch.z_pu = complex(branch.r_pu, branch.x_pu)
		branch.y_pu = 1/branch.z_pu if (branch.z_pu != 0) else float("inf")
		branch.start_node = system.nodes[branches.start[branch_idx]-1]
		branch.end_node = system.nodes[branches.end[branch_idx]-1]
		system.branches.append(branch)
	
	system.Ymatrix_calc()
	return system

""" Insert here per unit values of the grid for power and voltage """
S = 100*(10**6)
V = (11*(10**3))/math.sqrt(3)
slackV = 1.02

#execute powerflow analysis using nv_powerflow
Base = PerUnit(S,V)
branch, node = py_95bus_network_data.Network_95_nodes(Base, slackV)
Vtrue, Itrue, Iinjtrue, S1true, S2true, Sinjtrue, num_iter = nv_powerflow.solve(branch, node)

#execute powerflow analysis using nv_powerflow_cim
system = load_python_data(node, branch, node.type)
res, num_iter = nv_powerflow_cim.solve(system)

# Show numerical comparison with 6 decimals
"""
print("res.V==V? " +  str(not np.any(np.around(res.get_voltages(),6)-np.around(Vtrue,6))))
print("res.Iinj==Iinj? " +  str(not np.any(np.around(res.get_Iinj(),6)-np.around(Iinjtrue,6))))
print("res.I==I? " +  str(not np.any(np.around(res.getI(),6)-np.around(Itrue,6))))
print("res.S1==S1? " +  str(not np.any(np.around(res.get_S1(),6)-np.around(S1true,6))))
print("res.S2==S2? " +  str(not np.any(np.around(res.get_S2(),6)-np.around(S2true,6))))
print("res.Sinj==Sinj? " +  str(not np.any(np.around(res.get_Sinj(),6)-np.around(Sinjtrue,6))))
"""

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
	measurements_set.create_measurement(PowerflowNode.topology_node, nv_state_estimator_cim.ElemType.Node, nv_state_estimator_cim.MeasType.V_mag, np.absolute(PowerflowNode.voltage_pu), V_unc)
for i in range(1,len(res.nodes)):
	PowerflowNode = res.get_node(i)
	measurements_set.create_measurement(PowerflowNode.topology_node, nv_state_estimator_cim.ElemType.Node, nv_state_estimator_cim.MeasType.Sinj_real, PowerflowNode.power_pu.real, Sinj_unc)
for i in range(1,len(res.nodes)):
	PowerflowNode = res.get_node(i)
	measurements_set.create_measurement(PowerflowNode.topology_node, nv_state_estimator_cim.ElemType.Node, nv_state_estimator_cim.MeasType.Sinj_imag, PowerflowNode.power_pu.imag, Sinj_unc)
for i in [0,10,27,54,58]:
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch, nv_state_estimator_cim.MeasType.S1_real, res.branches[i].power_pu.real, S_unc)
for i in [9,53]:
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch, nv_state_estimator_cim.MeasType.S2_real, res.branches[i].power2_pu.real, S_unc)
for i in [0,10,27,54,58]:	
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch, nv_state_estimator_cim.MeasType.S1_imag, res.branches[i].power_pu.imag, S_unc)
for i in [9,53]:
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch, nv_state_estimator_cim.MeasType.S2_imag, res.branches[i].power2_pu.imag, S_unc)
for i in [12,35]:
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch, nv_state_estimator_cim.MeasType.I_mag, np.absolute(res.branches[i].current_pu), I_unc)
for i in [0,10,54]:
	PowerflowNode = res.get_node(i)
	measurements_set.create_measurement(PowerflowNode.topology_node, nv_state_estimator_cim.ElemType.Node, nv_state_estimator_cim.MeasType.Vpmu_mag, np.absolute(PowerflowNode.voltage_pu), Pmu_mag_unc)
for i in [0,10,54]:
	PowerflowNode = res.get_node(i)
	measurements_set.create_measurement(PowerflowNode.topology_node, nv_state_estimator_cim.ElemType.Node, nv_state_estimator_cim.MeasType.Vpmu_phase, np.angle(PowerflowNode.voltage_pu), Pmu_phase_unc)
for i in [12,35]:
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch, nv_state_estimator_cim.MeasType.Ipmu_mag, np.absolute(res.branches[i].current_pu), Pmu_mag_unc)
for i in [12,35]:
	measurements_set.create_measurement(res.branches[i].topology_branch, nv_state_estimator_cim.ElemType.Branch, nv_state_estimator_cim.MeasType.Ipmu_phase, np.angle(res.branches[i].current_pu), Pmu_phase_unc)


""" create measurements set (for nv_state_estimator.py)"""
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

# Perform state estimation
zdata, zdatameas = py_95bus_meas_data.Zdatatrue_creation(zdata, zdatameas, meas, branch, Vtrue, Itrue, Iinjtrue, S1true, S2true, Sinjtrue)
err_pu, zdatameas = py_95bus_meas_data.Zdatameas_creation(zdata, zdatameas)
measurements_set.meas_creation_test(err_pu)
#print("mval:")
#print(measurements_set.getmVal_test())
Vest, Iest, Iinjest, S1est, S2est, Sinjest = nv_state_estimator.DsseCall(branch, node, zdatameas, system.Ymatrix, system.Adjacencies)
results = nv_state_estimator_cim.DsseCall(system, measurements_set) 

# Show numerical comparison
print("results.V==Vest?: " + str(np.array_equal(np.around(Vest.complex,5), np.around(results.get_voltages(),5))))
print("results.Iinj==Iinjest? " +  str(not np.any(np.around(results.get_Iinj(),5)-np.around(Iinjest.complex,5))))
print("results.I==Iest? " +  str(not np.any(np.around(results.getI(),5)-np.around(Iest.complex,5))))
print("results.S1==S1est? " +  str(not np.any(np.around(results.get_S1(),5)-np.around(S1est.complex,5))))
print("results.S2==S2est? " +  str(not np.any(np.around(results.get_S2(),5)-np.around(S2est.complex,5))))
print("results.Sinj==Sinjest? " +  str(not np.any(np.around(results.get_Sinj(),5)-np.around(Sinjest.complex,5))))