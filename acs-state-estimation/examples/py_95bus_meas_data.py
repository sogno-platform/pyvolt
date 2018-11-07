def Zdatatrue_creation(zdata, zdatameas, meas, branch, V, I, Iinj, S1, S2, Sinj):
    """ It gives the true values at the measurement points."""
    import numpy
    
    t1 = numpy.ones(meas.V.num)
    t2 = 2*numpy.ones(meas.Sinj.num)
    t3 = 3*numpy.ones(meas.Sinj.num)
    t4 = 4*numpy.ones(meas.S1.num + meas.S2.num)
    t5 = 5*numpy.ones(meas.S1.num + meas.S2.num)
    t6 = 6*numpy.ones(meas.I.num)
    t7 = 7*numpy.ones(meas.Vpmu_mag.num)
    t8 = 8*numpy.ones(meas.Vpmu_phase.num)
    t9 = 9*numpy.ones(meas.Ipmu_mag.num)
    t10 = 10*numpy.ones(meas.Ipmu_phase.num)
    
    zdata.mtype = numpy.concatenate((t1,t2,t3,t4,t5,t6,t7,t8,t9,t10),axis=0)
    zdatameas.mtype = numpy.concatenate((t1,t2,t3,t4,t5,t6,t7,t8,t9,t10),axis=0)
    
    z1 = V.mag[meas.V.index-1]
    z2 = Sinj.real[meas.Sinj.index-1]
    z3 = Sinj.imag[meas.Sinj.index-1]
    z4_1 = S1.real[meas.S1.index-1]
    z4_2 = S2.real[meas.S2.index-1]
    z5_1 = S1.imag[meas.S1.index-1]
    z5_2 = S2.imag[meas.S2.index-1]
    z6 = I.mag[meas.I.index-1]
    z7 = V.mag[meas.Vpmu_mag.index-1]
    z8 = V.phase[meas.Vpmu_phase.index-1]
    z9 = I.mag[meas.Ipmu_mag.index-1]
    z10 = I.phase[meas.Ipmu_phase.index-1]
    
    zdata.mval = numpy.concatenate((z1,z2,z3,z4_1,z4_2,z5_1,z5_2,z6,z7,z8,z9,z10),axis=0)
    
    b1 = numpy.zeros(meas.V.num)
    b2 = numpy.zeros(meas.Sinj.num)
    b3 = numpy.zeros(meas.Sinj.num)
    b4_1 = branch.code[meas.S1.index-1]
    b4_2 = -branch.code[meas.S2.index-1]
    b5_1 = branch.code[meas.S1.index-1]
    b5_2 = -branch.code[meas.S2.index-1]
    b6 = branch.code[meas.I.index-1]
    b7 = numpy.zeros(meas.Vpmu_mag.num)
    b8 = numpy.zeros(meas.Vpmu_phase.num)
    b9 = branch.code[meas.Ipmu_mag.index-1]
    b10 = branch.code[meas.Ipmu_phase.index-1]
    
    zdata.mbranch = numpy.concatenate((b1,b2,b3,b4_1,b4_2,b5_1,b5_2,b6,b7,b8,b9,b10),axis=0)
    zdata.mbranch = zdata.mbranch.astype(int)
    zdatameas.mbranch = numpy.concatenate((b1,b2,b3,b4_1,b4_2,b5_1,b5_2,b6,b7,b8,b9,b10),axis=0)
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
    
    zdata.mfrom = numpy.concatenate((fr1,fr2,fr3,fr4_1,fr4_2,fr5_1,fr5_2,fr6,fr7,fr8,fr9,fr10),axis=0)
    zdata.mfrom = zdata.mfrom.astype(int)
    zdatameas.mfrom = numpy.concatenate((fr1,fr2,fr3,fr4_1,fr4_2,fr5_1,fr5_2,fr6,fr7,fr8,fr9,fr10),axis=0)
    zdatameas.mfrom = zdata.mfrom.astype(int)
    
    to1 = numpy.zeros(meas.V.num)
    to2 = numpy.zeros(meas.Sinj.num)
    to3 = numpy.zeros(meas.Sinj.num)
    to4_1 = branch.end[meas.S1.index-1]
    to4_2 = branch.start[meas.S2.index-1]
    to5_1 = branch.end[meas.S1.index-1]
    to5_2 = branch.start[meas.S2.index-1]
    to6 = branch.end[meas.I.index-1]
    to7 = numpy.zeros(meas.Vpmu_mag.num)
    to8 = numpy.zeros(meas.Vpmu_phase.num)
    to9 = branch.end[meas.Ipmu_mag.index-1]
    to10 = branch.end[meas.Ipmu_phase.index-1]
    
    zdata.mto = numpy.concatenate((to1,to2,to3,to4_1,to4_2,to5_1,to5_2,to6,to7,to8,to9,to10),axis=0)
    zdata.mto = zdata.mto.astype(int)
    zdatameas.mto = numpy.concatenate((to1,to2,to3,to4_1,to4_2,to5_1,to5_2,to6,to7,to8,to9,to10),axis=0)
    zdatameas.mto = zdata.mto.astype(int)
    
    d1 = z1*(meas.V.unc/300)
    d2 = numpy.absolute(z2*(meas.Sinj.unc/300))
    d3 = numpy.absolute(z3*(meas.Sinj.unc/300))
    d4_1 = numpy.absolute(z4_1*(meas.S1.unc/300))
    d4_2 = numpy.absolute(z4_2*(meas.S2.unc/300))
    d5_1 = numpy.absolute(z5_1*(meas.S1.unc/300))
    d5_2 = numpy.absolute(z5_2*(meas.S2.unc/300))
    d6 = z6*(meas.I.unc/300)
    d7 = z7*(meas.Vpmu_mag.unc/300)
    d8 = (meas.Vpmu_phase.unc/300)*numpy.ones(meas.Vpmu_phase.num)
    d9 = z9*(meas.Ipmu_mag.unc/300)
    d10 = (meas.Ipmu_phase.unc/300)*numpy.ones(meas.Ipmu_phase.num)
    
    zdata.mstddev = numpy.concatenate((d1,d2,d3,d4_1,d4_2,d5_1,d5_2,d6,d7,d8,d9,d10),axis=0)
    zdatameas.mstddev = numpy.concatenate((d1,d2,d3,d4_1,d4_2,d5_1,d5_2,d6,d7,d8,d9,d10),axis=0)
#    idx = numpy.where(zdata.mstddev<10**(-6))
#    zdata.mstddev[idx] = 10**(-6)
    
    return zdata, zdatameas

def Zdatameas_creation(zdata, zdatameas):
    """ It gives the measured values (affected by uncertainty) at the measurement points."""
    import numpy
    
    zmeas = zdata.mval
    zdev = zdata.mstddev
    err_pu = numpy.random.normal(0,1,len(zmeas))
    zdatameas.mval = zmeas + numpy.multiply(zdev,err_pu)
    
    return zdatameas