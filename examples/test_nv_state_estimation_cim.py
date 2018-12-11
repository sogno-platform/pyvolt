import numpy as np
import math
import sys
import copy 

sys.path.append("../acs/state_estimation")
import py_95bus_network_data
import py_95bus_meas_data
import nv_state_estimator
import nv_powerflow

import network
import nv_powerflow_cim
import cim_py_95bus_meas_data
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
system = network.load_python_data(node, branch, node.type)
Vtrue, Itrue, Iinjtrue, S1true, S2true, Sinjtrue, num_iter = nv_powerflow.solve(branch, node)
Vtrue_cim, Itrue_cim, Iinjtrue_cim, S1true_cim, S2true_cim, Sinjtrue_cim, num_iter_cim = nv_powerflow_cim.solve(system)

""" Write here the indexes of the nodes/branches where the measurements are"""
V_idx = np.array([1,11,55])
I_idx = np.array([13,36])
Sinj_idx = np.linspace(2,node.num,node.num-1)
S1_idx = np.array([1,11,28,55,59])
S2_idx = np.array([10,54])
Ipmu_idx = np.array([1,11,55])
Vpmu_idx = np.array([13,36])
#Ipmu_idx = np.array([])
#Vpmu_idx = np.array([])

""" Write here the percent uncertainties of the measurements""" 
V_unc = 1
I_unc = 2
Sinj_unc = 5
S_unc = 2
Pmu_mag_unc = 0.7
Pmu_phase_unc = 0.7

V = Measurements(V_idx,V_unc)
I = Measurements(I_idx,I_unc)
Sinj = Measurements(Sinj_idx,Sinj_unc)
S1 = Measurements(S1_idx,S_unc)
S2 = Measurements(S2_idx,S_unc)
Ipmu_mag = Measurements(Ipmu_idx,Pmu_mag_unc)
Ipmu_phase = Measurements(Ipmu_idx,Pmu_phase_unc)
Vpmu_mag = Measurements(Vpmu_idx,Pmu_mag_unc)
Vpmu_phase = Measurements(Vpmu_idx,Pmu_phase_unc)
meas = Measurement_set(V, I, Sinj, S1, S2, Ipmu_mag, Ipmu_phase, Vpmu_mag, Vpmu_phase)
zdata = Zdata_init(meas)
zdatameas =  Zdata_init(meas)
zdata_cim = copy.deepcopy(zdata) 
zdatameas_cim =  copy.deepcopy(zdatameas)

zdata, zdatameas = py_95bus_meas_data.Zdatatrue_creation(zdata, zdatameas, meas, branch, Vtrue, Itrue, Iinjtrue, S1true, S2true, Sinjtrue)
zdata_cim, zdatameas_cim = cim_py_95bus_meas_data.Zdatatrue_creation(zdata_cim, zdatameas_cim, meas, system, Vtrue_cim, Itrue_cim, Iinjtrue_cim, S1true_cim, S2true_cim, Sinjtrue_cim)

#compare results:
print("Vtrue==Vtrue_cim?: " + str((Vtrue==Vtrue_cim).all()))
print("Itrue==Itrue_cim?: " + str((Itrue==Itrue_cim).all()))
print("Iinjtrue==Iinjtrue_cim?: " + str((Iinjtrue==Iinjtrue_cim).all()))
print("S1true==S1true_cim?: " + str((S1true==S1true_cim).all()))
print("S2true==S2true_cim?: " + str((S2true==S2true_cim).all()))
print("Sinjtrue==Sinjtrue_cim?: " + str((Sinjtrue==Sinjtrue_cim).all()))
print("zdata.mtype==zdata_cim.mtype?: {}".format(np.array_equal(zdata.mtype, zdata_cim.mtype)))
print("zdata.mval==zdata_cim.mval?: {}".format(np.array_equal(zdata.mval, zdata_cim.mval)))
print("zdatameas.mval==zdatameas_cim.mval?: {}".format(np.array_equal(zdatameas.mval, zdatameas_cim.mval)))
#print("zdata.mbranch==zdata_cim.mbranch?: {}".format(np.array_equal(zdata.mbranch, zdata_cim.mbranch)))
print("zdata.mfrom==zdata_cim.mfrom?: {}".format(np.array_equal(zdata.mfrom, zdata_cim.mfrom)))
#print(zdata.mfrom-zdata_cim.mfrom)
print("zdata.mto==zdata_cim.mto?: {}".format(np.array_equal(zdata.mto, zdata_cim.mto)))
#print(zdata.mto-zdata_cim.mto)
print("zdata.mstddev==zdata_cim.mstddev?: {}".format(np.array_equal(zdata.mstddev, zdata_cim.mstddev)))


MC_trials = 1
iter_counter = np.zeros(MC_trials)
Ymatr, Adj = nv_powerflow.Ymatrix_calc(branch,node)
Ymatr_cim, Adj_cim = network.Ymatrix_calc(system)

for mciter in range(0,MC_trials):      
	zdatameas = py_95bus_meas_data.Zdatameas_creation(zdata, zdatameas)
	zdatameas_cim=copy.deepcopy(zdatameas) 		#to compare use the same random normal (Gaussian) distribution.
	#zdatameas_cim = cim_py_95bus_meas_data.Zdatameas_creation(zdata_cim, zdatameas_cim)
	Vest, Iest, Iinjest, S1est, S2est, Sinjest = nv_state_estimator.DsseCall(branch, node, zdatameas, Ymatr, Adj)
	Vest_cim, Iest_cim, Iinjest_cim, S1est_cim, S2est_cim, Sinjest_cim = nv_state_estimator_cim.DsseCall(system, zdatameas_cim, Ymatr_cim, Adj_cim)

print("Ymatr==Ymatr_cim?: {}".format(np.array_equal(Ymatr, Ymatr_cim)))
print("Adj==Adj_cim?: {}".format(np.array_equal(Adj, Adj_cim)))
print("Vest==Vest_cim?: {}".format(np.array_equal(Vest.complex, Vest_cim.complex)))
#print(Vest.complex - Vest_cim.complex)
print("Iest==Iest_cim?: {}".format(np.array_equal(Iest.complex, Iest_cim.complex)))
print("Iinjest==Iinjest_cim?: {}".format(np.array_equal(Iinjest.complex, Iinjest_cim.complex)))
print("S1est==S1est_cim?: {}".format(np.array_equal(S1est.complex, S1est_cim.complex)))
print("S2est==S2est_cim?: {}".format(np.array_equal(S2est.complex, S2est_cim.complex)))
print("Sinjest==Sinjest_cim?: {}".format(np.array_equal(Sinjest.complex, Sinjest_cim.complex)))
print(Vest.complex)