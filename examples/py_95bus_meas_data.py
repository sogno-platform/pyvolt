def Zdatatrue_creation(zdata, zdatameas, meas, branch, V, I, Iinj, S1, S2, Sinj):
    """ It gives the true values at the measurement points."""
    import numpy as np
    
    t1 = np.ones(meas.V.num)
    t2 = 2*np.ones(meas.Sinj.num)
    t3 = 3*np.ones(meas.Sinj.num)
    t4 = 4*np.ones(meas.S1.num + meas.S2.num)
    t5 = 5*np.ones(meas.S1.num + meas.S2.num)
    t6 = 6*np.ones(meas.I.num)
    t7 = 7*np.ones(meas.Vpmu_mag.num)
    t8 = 8*np.ones(meas.Vpmu_phase.num)
    t9 = 9*np.ones(meas.Ipmu_mag.num)
    t10 = 10*np.ones(meas.Ipmu_phase.num)
    
    zdata.mtype = np.concatenate((t1,t2,t3,t4,t5,t6,t7,t8,t9,t10),axis=0)
    zdatameas.mtype = np.concatenate((t1,t2,t3,t4,t5,t6,t7,t8,t9,t10),axis=0)
    
    z1 = np.abs(V[meas.V.index-1])
    z2 = Sinj.real[meas.Sinj.index-1]
    z3 = Sinj.imag[meas.Sinj.index-1]
    z4_1 = S1.real[meas.S1.index-1]
    z4_2 = S2.real[meas.S2.index-1]
    z5_1 = S1.imag[meas.S1.index-1]
    z5_2 = S2.imag[meas.S2.index-1]
    z6 = np.abs(I[meas.I.index-1])
    z7 = np.abs(V[meas.Vpmu_mag.index-1])
    z8 = np.angle(V[meas.Vpmu_phase.index-1])
    z9 = np.abs(I[meas.Ipmu_mag.index-1])
    z10 = np.angle(I[meas.Ipmu_phase.index-1])
    
    zdata.mval = np.concatenate((z1,z2,z3,z4_1,z4_2,z5_1,z5_2,z6,z7,z8,z9,z10),axis=0)
    
    b1 = np.zeros(meas.V.num)
    b2 = np.zeros(meas.Sinj.num)
    b3 = np.zeros(meas.Sinj.num)
    b4_1 = branch.code[meas.S1.index-1]
    b4_2 = -branch.code[meas.S2.index-1]
    b5_1 = branch.code[meas.S1.index-1]
    b5_2 = -branch.code[meas.S2.index-1]
    b6 = branch.code[meas.I.index-1]
    b7 = np.zeros(meas.Vpmu_mag.num)
    b8 = np.zeros(meas.Vpmu_phase.num)
    b9 = branch.code[meas.Ipmu_mag.index-1]
    b10 = branch.code[meas.Ipmu_phase.index-1]
    
    zdata.mbranch = np.concatenate((b1,b2,b3,b4_1,b4_2,b5_1,b5_2,b6,b7,b8,b9,b10),axis=0)
    zdata.mbranch = zdata.mbranch.astype(int)
    zdatameas.mbranch = np.concatenate((b1,b2,b3,b4_1,b4_2,b5_1,b5_2,b6,b7,b8,b9,b10),axis=0)
    zdatameas.mbranch = zdata.mbranch.astype(int)
    
    fr1 = meas.V.index
    fr2 = meas.Sinj.index
    fr3 = meas.Sinj.index
    fr4_1 = branch.start[meas.S1.index-1]
    fr4_2 = branch.end[meas.S2.index-1]
    fr5_1 = branch.start[meas.S1.index-1]
    fr5_2 = branch.end[meas.S2.index-1]
    fr6 = branch.start[meas.I.index-1]
    fr7 = meas.Vpmu_mag.index
    fr8 = meas.Vpmu_phase.index
    fr9 = branch.start[meas.Ipmu_mag.index-1]
    fr10 = branch.start[meas.Ipmu_phase.index-1]
    
    zdata.mfrom = np.concatenate((fr1,fr2,fr3,fr4_1,fr4_2,fr5_1,fr5_2,fr6,fr7,fr8,fr9,fr10),axis=0)
    zdata.mfrom = zdata.mfrom.astype(int)
    zdatameas.mfrom = np.concatenate((fr1,fr2,fr3,fr4_1,fr4_2,fr5_1,fr5_2,fr6,fr7,fr8,fr9,fr10),axis=0)
    zdatameas.mfrom = zdata.mfrom.astype(int)
    
    to1 = np.zeros(meas.V.num)
    to2 = np.zeros(meas.Sinj.num)
    to3 = np.zeros(meas.Sinj.num)
    to4_1 = branch.end[meas.S1.index-1]
    to4_2 = branch.start[meas.S2.index-1]
    to5_1 = branch.end[meas.S1.index-1]
    to5_2 = branch.start[meas.S2.index-1]
    to6 = branch.end[meas.I.index-1]
    to7 = np.zeros(meas.Vpmu_mag.num)
    to8 = np.zeros(meas.Vpmu_phase.num)
    to9 = branch.end[meas.Ipmu_mag.index-1]
    to10 = branch.end[meas.Ipmu_phase.index-1]
    
    zdata.mto = np.concatenate((to1,to2,to3,to4_1,to4_2,to5_1,to5_2,to6,to7,to8,to9,to10),axis=0)
    zdata.mto = zdata.mto.astype(int)
    zdatameas.mto = np.concatenate((to1,to2,to3,to4_1,to4_2,to5_1,to5_2,to6,to7,to8,to9,to10),axis=0)
    zdatameas.mto = zdata.mto.astype(int)
    
    d1 = z1*(meas.V.unc/300)
    d2 = np.absolute(z2*(meas.Sinj.unc/300))
    d3 = np.absolute(z3*(meas.Sinj.unc/300))
    d4_1 = np.absolute(z4_1*(meas.S1.unc/300))
    d4_2 = np.absolute(z4_2*(meas.S2.unc/300))
    d5_1 = np.absolute(z5_1*(meas.S1.unc/300))
    d5_2 = np.absolute(z5_2*(meas.S2.unc/300))
    d6 = z6*(meas.I.unc/300)
    d7 = z7*(meas.Vpmu_mag.unc/300)
    d8 = (meas.Vpmu_phase.unc/300)*np.ones(meas.Vpmu_phase.num)
    d9 = z9*(meas.Ipmu_mag.unc/300)
    d10 = (meas.Ipmu_phase.unc/300)*np.ones(meas.Ipmu_phase.num)
    
    zdata.mstddev = np.concatenate((d1,d2,d3,d4_1,d4_2,d5_1,d5_2,d6,d7,d8,d9,d10),axis=0)
    zdatameas.mstddev = np.concatenate((d1,d2,d3,d4_1,d4_2,d5_1,d5_2,d6,d7,d8,d9,d10),axis=0)
#    idx = np.where(zdata.mstddev<10**(-6))
#    zdata.mstddev[idx] = 10**(-6)
    
    return zdata, zdatameas

def Zdatameas_creation(zdata, zdatameas):
    """ It gives the measured values (affected by uncertainty) at the measurement points."""
    import numpy as np
    
    zmeas = zdata.mval
    zdev = zdata.mstddev
    err_pu = np.random.normal(0,1,len(zmeas))
    zdatameas.mval = zmeas + np.multiply(zdev,err_pu)
    
    return zdatameas