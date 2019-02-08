import numpy as np

class Measurements:
	def __init__(self, index, unc):
		self.index = index
		self.unc = unc
		self.num = len(index)
		
class Measurement_set:
	def __init__(self, V, I, Pinj, P1, P2, Qinj, Q1, Q2, Vpmu_mag, Vpmu_phase, Ipmu_mag, Ipmu_phase):
		self.V = V
		self.Pinj = Pinj
		self.P1 = P1
		self.P2 = P2
		self.Qinj = Qinj
		self.Q1 = Q1
		self.Q2 = Q2
		self.I = I
		self.Vpmu_mag = Vpmu_mag
		self.Vpmu_phase = Vpmu_phase 
		self.Ipmu_mag = Ipmu_mag
		self.Ipmu_phase = Ipmu_phase
		
class Zdata_init:
	def __init__(self, meas):
		nmeas = meas.V.num + meas.I.num + meas.Pinj.num + meas.P1.num + meas.P2.num + meas.Qinj.num + meas.Q1.num + meas.Q2.num + meas.Ipmu_mag.num + meas.Ipmu_phase.num + meas.Vpmu_mag.num + meas.Vpmu_phase.num
		self.mtype = np.zeros(nmeas)
		self.mval = np.zeros(nmeas)
		self.mbranch = np.zeros(nmeas)
		self.mfrom = np.zeros(nmeas)
		self.mto = np.zeros(nmeas)
		self.mstddev = np.zeros(nmeas)

def Zdata_structure_creation(data, system):
    """ It gives the structure of the zdata matrix used within the state estimator."""
		
    """ Writing here the indexes of the nodes/branches where the measurements are"""
    V_idx = data["Measurement"]["Vmag"]["idx"]
    I_idx = data["Measurement"]["Imag"]["idx"]
    Pinj_idx = data["Measurement"]["Pinj"]["idx"]
    P1_idx = data["Measurement"]["P1"]["idx"]
    P2_idx = data["Measurement"]["P2"]["idx"]
    Qinj_idx = data["Measurement"]["Qinj"]["idx"]
    Q1_idx = data["Measurement"]["Q1"]["idx"]
    Q2_idx = data["Measurement"]["Q2"]["idx"]
    Ipmu_idx = data["Measurement"]["Ipmu"]["idx"]
    Vpmu_idx = data["Measurement"]["Vpmu"]["idx"]

    """ Writing here the percent uncertainties of the measurements""" 
    V_unc = data["Measurement"]["Vmag"]["unc"]
    I_unc = data["Measurement"]["Imag"]["unc"]
    Pinj_unc = data["Measurement"]["Pinj"]["unc"]
    Qinj_unc = data["Measurement"]["Qinj"]["unc"]
    P1_unc = data["Measurement"]["P1"]["unc"]
    Q1_unc = data["Measurement"]["Q1"]["unc"]
    P2_unc = data["Measurement"]["P2"]["unc"]
    Q2_unc = data["Measurement"]["Q2"]["unc"]
    Ipmu_mag_unc = data["Measurement"]["Ipmu"]["unc_mag"]
    Ipmu_phase_unc = data["Measurement"]["Ipmu"]["unc_phase"]
    Vpmu_mag_unc = data["Measurement"]["Vpmu"]["unc_mag"]
    Vpmu_phase_unc = data["Measurement"]["Vpmu"]["unc_phase"]
    V = Measurements(V_idx,V_unc)
    I = Measurements(I_idx,I_unc)
    Pinj = Measurements(Pinj_idx,Pinj_unc)
    P1 = Measurements(P1_idx,P1_unc)
    P2 = Measurements(P2_idx,P2_unc)
    Qinj = Measurements(Qinj_idx,Qinj_unc)
    Q1 = Measurements(Q1_idx,Q1_unc)
    Q2 = Measurements(Q2_idx,Q2_unc)
    Ipmu_mag = Measurements(Ipmu_idx,Ipmu_mag_unc)
    Ipmu_phase = Measurements(Ipmu_idx,Ipmu_phase_unc)
    Vpmu_mag = Measurements(Vpmu_idx,Vpmu_mag_unc)
    Vpmu_phase = Measurements(Vpmu_idx,Vpmu_phase_unc)
    meas = Measurement_set(V, I, Pinj, P1, P2, Qinj, Q1, Q2, Vpmu_mag, Vpmu_phase, Ipmu_mag, Ipmu_phase)

    zdatameas = Zdata_init(meas)
    zdatameas = Zdatatrue_creation(zdatameas, meas, system)
	
    return zdatameas

def Zdatatrue_creation(zdatameas, meas, system):
    """ It gives the true values at the measurement points."""
    import numpy as np
    
    t1 = np.ones(meas.V.num)
    t2 = 2*np.ones(meas.Pinj.num)
    t3 = 3*np.ones(meas.Qinj.num)
    t4 = 4*np.ones(meas.P1.num + meas.P2.num)
    t5 = 5*np.ones(meas.Q1.num + meas.Q2.num)
    t6 = 6*np.ones(meas.I.num)
    t7 = 7*np.ones(meas.Vpmu_mag.num)
    t8 = 8*np.ones(meas.Vpmu_phase.num)
    t9 = 9*np.ones(meas.Ipmu_mag.num)
    t10 = 10*np.ones(meas.Ipmu_phase.num)
    
    zdatameas.mtype = np.concatenate((t1,t2,t3,t4,t5,t6,t7,t8,t9,t10),axis=0)
    
    # z1 = np.abs(V[meas.V.index-1])
    # z2 = Sinj.real[meas.Sinj.index-1]
    # z3 = Sinj.imag[meas.Sinj.index-1]
    # z4_1 = S1.real[meas.S1.index-1]
    # z4_2 = S2.real[meas.S2.index-1]
    # z5_1 = S1.imag[meas.S1.index-1]
    # z5_2 = S2.imag[meas.S2.index-1]
    # z6 = np.abs(I[meas.I.index-1])
    # z7 = np.abs(V[meas.Vpmu_mag.index-1])
    # z8 = np.angle(V[meas.Vpmu_phase.index-1])
    # z9 = np.abs(I[meas.Ipmu_mag.index-1])
    # z10 = np.angle(I[meas.Ipmu_phase.index-1])
    
    # zdata.mval = np.concatenate((z1,z2,z3,z4_1,z4_2,z5_1,z5_2,z6,z7,z8,z9,z10),axis=0)
    
    fr1 = meas.V.index
    fr2 = meas.Pinj.index
    fr3 = meas.Qinj.index
    fr4_1 = np.array([])
    fr4_2 = np.array([])
    fr5_1 = np.array([])
    fr5_2 = np.array([])
    for index in range(len(meas.P1.index)):
        fr4_1 = np.append(fr4_1, system.branches[meas.P1.index[index]-1].start_node.index)
    for index in range(len(meas.Q1.index)):
        fr5_1 = np.append(fr5_1, system.branches[meas.Q1.index[index]-1].start_node.index)
    for index in range(len(meas.P2.index)):
        fr4_2 = np.append(fr4_2, system.branches[meas.P2.index[index]-1].end_node.index)
    for index in range(len(meas.Q2.index)):
        fr5_2 = np.append(fr5_2, system.branches[meas.Q2.index[index]-1].end_node.index)
    fr6 = np.array([])
    for index in range(len(meas.I.index)):
        fr6 = np.append(fr6, system.branches[meas.I.index[index]-1].start_node.index)
    fr7 = meas.Vpmu_mag.index
    fr8 = meas.Vpmu_phase.index
    fr9 = np.array([])
    fr10 = np.array([])
    for index in range(len(meas.Ipmu_mag.index)):
        fr9 = np.append(fr9, system.branches[meas.Ipmu_mag.index[index]-1].start_node.index)
    for index in range(len(meas.Ipmu_phase.index)):
        fr10 = np.append(fr10, system.branches[meas.Ipmu_phase.index[index]-1].start_node.index)	
    
    zdatameas.mfrom = np.concatenate((fr1,fr2,fr3,fr4_1,fr4_2,fr5_1,fr5_2,fr6,fr7,fr8,fr9,fr10),axis=0)
    zdatameas.mfrom = zdatameas.mfrom.astype(int)
    
    to1 = np.zeros(meas.V.num)
    to2 = np.zeros(meas.Pinj.num)
    to3 = np.zeros(meas.Qinj.num)
    to4_1 = np.array([])
    to4_2 = np.array([])
    to5_1 = np.array([])
    to5_2 = np.array([])
    for index in range(len(meas.P1.index)):
        to4_1 = np.append(to4_1, system.branches[meas.P1.index[index]-1].end_node.index)
    for index in range(len(meas.Q1.index)):
        to5_1 = np.append(to5_1, system.branches[meas.Q1.index[index]-1].end_node.index)
    for index in range(len(meas.P2.index)):
        to4_2 = np.append(to4_2, system.branches[meas.P2.index[index]-1].start_node.index)
    for index in range(len(meas.Q2.index)):
        to5_2 = np.append(to5_2, system.branches[meas.Q2.index[index]-1].start_node.index)
    to6 = np.array([])
    for index in range(len(meas.I.index)):
        to6 = np.append(to6, system.branches[meas.I.index[index]-1].end_node.index)
    to7 = np.zeros(meas.Vpmu_mag.num)
    to8 = np.zeros(meas.Vpmu_phase.num)
    to9 = np.array([])
    to10 = np.array([])
    for index in range(len(meas.Ipmu_mag.index)):
        to9 = np.append(to9, system.branches[meas.Ipmu_mag.index[index]-1].end_node.index)
    for index in range(len(meas.Ipmu_phase.index)):
        to10 = np.append(to10, system.branches[meas.Ipmu_phase.index[index]-1].end_node.index)
	
    zdatameas.mto = np.concatenate((to1,to2,to3,to4_1,to4_2,to5_1,to5_2,to6,to7,to8,to9,to10),axis=0)
    zdatameas.mto = zdatameas.mto.astype(int)
    
    d1 = (meas.V.unc/300)*np.ones(meas.V.num)
    d2 = (meas.Pinj.unc/300)*np.ones(meas.Pinj.num)
    d3 = (meas.Qinj.unc/300)*np.ones(meas.Qinj.num)
    d4_1 = (meas.P1.unc/300)*np.ones(meas.P1.num)
    d4_2 = (meas.P2.unc/300)*np.ones(meas.P2.num)
    d5_1 = (meas.Q1.unc/300)*np.ones(meas.Q1.num)
    d5_2 = (meas.Q2.unc/300)*np.ones(meas.Q2.num)
    d6 = (meas.I.unc/300)*np.ones(meas.I.num)
    d7 = (meas.Vpmu_mag.unc/300)*np.ones(meas.Vpmu_mag.num)
    d8 = (meas.Vpmu_phase.unc/300)*np.ones(meas.Vpmu_phase.num)
    d9 = (meas.Ipmu_mag.unc/300)*np.ones(meas.Ipmu_mag.num)
    d10 = (meas.Ipmu_phase.unc/300)*np.ones(meas.Ipmu_phase.num)
    
    zdatameas.mstddev = np.concatenate((d1,d2,d3,d4_1,d4_2,d5_1,d5_2,d6,d7,d8,d9,d10),axis=0)
#	idx = np.where(zdata.mstddev<10**(-6))
#	zdata.mstddev[idx] = 10**(-6)
	
    return zdatameas

def Zdatameas_creation_fromPF(data, zdatameas, values):
    """ It gives the measured values (affected by uncertainty) at the measurement points."""
    import numpy as np
    
    num_data = len(values)
    index = np.linspace(0,num_data-2,int(num_data/2)).astype(int)
    
    Vmag = [values[i] for i in index]
    Vtheta = [values[i] for i in index + 1]
    Pinj = [];
    Qinj = [];
    P1 = [];
    P2 = [];
    Q1 = [];
    Q2 = [];
    Imag = [];
    Itheta = [];
    
    z1 = [Vmag[i-1] for i in data["Measurement"]["Vmag"]["idx"]]
    z2 = [Pinj[i-1] for i in data["Measurement"]["Pinj"]["idx"]]
    z3 = [Qinj[i-1] for i in data["Measurement"]["Qinj"]["idx"]]
    z4_1 = [P1[i-1] for i in data["Measurement"]["P1"]["idx"]]
    z4_2 = [P2[i-1] for i in data["Measurement"]["P2"]["idx"]]
    z5_1 = [Q1[i-1] for i in data["Measurement"]["Q1"]["idx"]]
    z5_2 = [Q2[i-1] for i in data["Measurement"]["Q2"]["idx"]]
    z6 = [Imag[i-1] for i in data["Measurement"]["Imag"]["idx"]]
    z7 = [Vmag[i-1] for i in data["Measurement"]["Vpmu"]["idx"]]
    z8 = [Vtheta[i-1] for i in data["Measurement"]["Vpmu"]["idx"]]
    z9 = [Imag[i-1] for i in data["Measurement"]["Ipmu"]["idx"]]
    z10 = [Itheta[i-1] for i in data["Measurement"]["Ipmu"]["idx"]]
    
    zmeas = np.concatenate((z1,z2,z3,z4_1,z4_2,z5_1,z5_2,z6,z7,z8,z9,z10),axis=0)
    zdatameas.mval = zmeas
    zdev = np.multiply(zmeas,zdatameas.mstddev)
    vthetaidx = np.where(zdatameas.mtype==8)
    ithetaidx = np.where(zdatameas.mtype==10)
    zdev[vthetaidx] = zdatameas.mstddev[vthetaidx]
    zdev[ithetaidx] = zdatameas.mstddev[ithetaidx]
    err_pu = np.random.normal(0,1,len(zmeas))
    zdatameas.mval = zmeas + np.multiply(zdev,err_pu)
    
    return zdatameas