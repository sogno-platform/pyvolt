import numpy
import math
import matplotlib as plt
import Network_data
import Power_flow
import Meas_data
import NvStateEstimator
#import BcStateEstimator

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
        self.mtype = numpy.zeros(nmeas)
        self.mval = numpy.zeros(nmeas)
        self.mbranch = numpy.zeros(nmeas)
        self.mfrom = numpy.zeros(nmeas)
        self.mto = numpy.zeros(nmeas)
        self.mstddev = numpy.zeros(nmeas)
        
""" Insert here per unit values of the grid for power and voltage """
S = 100*(10**6)
V = (11*(10**3))/math.sqrt(3)
slackV = 1.02

Base = PerUnit(S,V)
branch, node = Network_data.Network_95_nodes(Base, slackV)

Vtrue, Itrue, Iinjtrue, S1true, S2true, Sinjtrue, num_iter = Power_flow.NV_power_flow(branch, node)

""" Write here the indexes of the nodes/branches where the measurements are"""
V_idx = numpy.array([1,11,55])
I_idx = numpy.array([13,36])
Sinj_idx = numpy.linspace(2,node.num,node.num-1)
S1_idx = numpy.array([1,11,28,55,59])
S2_idx = numpy.array([10,54])
Ipmu_idx = numpy.array([1,11,55])
Vpmu_idx = numpy.array([13,36])

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

zdata, zdatameas = Meas_data.Zdatatrue_creation(zdata, zdatameas, meas, branch, Vtrue, Itrue, Iinjtrue, S1true, S2true, Sinjtrue)

MC_trials = 1

#Iest_matrix = numpy.zeros((MC_trials,branch.num),dtype=numpy.complex_)
#Vest_matrix = numpy.zeros((MC_trials,node.num),dtype=numpy.complex_)
#S1est_matrix = numpy.zeros((MC_trials,branch.num),dtype=numpy.complex_)
#S2est_matrix = numpy.zeros((MC_trials,branch.num),dtype=numpy.complex_)
#Sinjest_matrix = numpy.zeros((MC_trials,node.num),dtype=numpy.complex_)
iter_counter = numpy.zeros(MC_trials)

Ymatr, Adj = Power_flow.Ymatrix_calc(branch,node)

for mciter in range(0,MC_trials):
    if mciter%10 == 0:
        print(mciter)
        
    zdatameas = Meas_data.Zdatameas_creation(zdata, zdatameas)
  
    Vest, Iest, Iinjest, S1est, S2est, Sinjest = NvStateEstimator.DsseCall(branch, node, zdatameas, Ymatr, Adj)
    
#    Iest_matrix[mciter] = Iest.complex
#    Vest_matrix[mciter] = Vest.complex
#    S1est_matrix[mciter] = S1est.complex
#    S2est_matrix[mciter] = S2est.complex
#    Sinjest_matrix[mciter] = Sinjest.complex
    