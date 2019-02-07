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
results = nv_powerflow_cim.PowerflowResults(system)
Vtrue, Itrue, Iinjtrue, S1true, S2true, Sinjtrue, num_iter = nv_powerflow.solve(branch, node)
Vtrue_cim, num_iter_cim = nv_powerflow_cim.solve(system)
results.load_voltages(Vtrue_cim)
results.calculateI()
results.calculateIinj()
results.calculateSinj()

print("Vtrue==Vtrue_cim?: " + str((results.get_voltages()==Vtrue).all()))
print("Iinj==Iinj_cim?: " + str((results.get_Iinj()==Iinjtrue).all()))
print("Sinj==Sinj_cim?: " + str((results.get_Sinj()==Sinjtrue).all()))