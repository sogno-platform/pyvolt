import numpy

class Real_to_all:
    def __init__(self,Ar, Ax):
        """ Converts real and imaginary parts to all the other formats"""
        self.real = Ar
        self.imag = Ax
        self.complex = Ar + 1j*Ax
        self.mag = numpy.absolute(self.complex)
        self.phase = numpy.angle(self.complex)
            
class Znod:
    def __init__(self,branch,node):
        """ Calculates a matrix of impedances for calculating the voltage in all
        the nodes of the grid, starting from the slack bus voltage and knowing
        the branch currents.
        
        In each row it gives the set of branch impedances to calculate the voltage 
        of the corresponding node through the voltage drops starting from the 
        first node of the grid """
        
        self.Z = numpy.zeros((node.num-1,branch.num),dtype=numpy.complex_)
        for k in range(1,node.num):
            i = k-1
            aa = branch.start[i]-2
            if aa == -1:
                self.Z[i] = self.Z[i]
            else:
                self.Z[i] = self.Z[aa]
            self.Z[i][i] = branch.Z[i]
        zzz = numpy.zeros((1,branch.num))
        self.Z = numpy.concatenate((zzz,self.Z), axis=0)
        self.R = self.Z.real
        self.X = self.Z.imag

def DsseCall(branch, node, zdata):
    """ It identifies the type of measurements present in the measurement set and 
    calls the appropriate estimator for dealing with them."""
    
    trad_code = 0
    PMU_code = 0
        
    Vmag_meas = numpy.where(zdata.mtype==1)
    if len(Vmag_meas[0])>0:
        trad_code = 1
        
    Vpmu_meas = numpy.where(zdata.mtype==7)
    if len(Vpmu_meas[0])>0:
        PMU_code = 2
        
    est_code = trad_code + PMU_code
        
    if est_code == 1:
        Vest, Iest = DsseTrad(branch, node, zdata)
    elif est_code == 2:
        Vest, Iest = DssePmu(branch, node, zdata)
    else:
        Vest, Iest = DsseMixed(branch, node, zdata)
        
    """ From here on, all the other quantities of the grid are calculated. """        
    Iinj_r = numpy.zeros(node.num)
    Iinj_x = numpy.zeros(node.num)
    for k in range(1,node.num+1):
        to = numpy.where(branch.end==k)
        fr = numpy.where(branch.start==k)
        Iinj_r[k-1] = numpy.sum(Iest.real[to[0]]) - numpy.sum(Iest.real[fr[0]])
        Iinj_x[k-1] = numpy.sum(Iest.imag[to[0]]) - numpy.sum(Iest.imag[fr[0]])
        
    Iinjest = Real_to_all(Iinj_r,Iinj_x)
    Sinj_rx = numpy.multiply(Vest.complex,numpy.conj(Iinjest.complex))
    Sinjest = Real_to_all(numpy.real(Sinj_rx),numpy.imag(Sinj_rx))
    S1_rx = numpy.multiply(Vest.complex[branch.start-1],numpy.conj(Iest.complex))
    S2_rx = -numpy.multiply(Vest.complex[branch.end-1],numpy.conj(Iest.complex))
    S1est = Real_to_all(numpy.real(S1_rx),numpy.imag(S1_rx))
    S2est = Real_to_all(numpy.real(S2_rx),numpy.imag(S2_rx))
        
    return Vest, Iest, Iinjest, S1est, S2est, Sinjest


def DsseTrad(branch, node, zdata):
    """ It performs state estimation using branch current state variables and it is customized
    to work without PMU measurements. """
    
    vidx = numpy.where(zdata.mtype==1)
    pidx = numpy.where(zdata.mtype==2)
    qidx = numpy.where(zdata.mtype==3)
    pfidx = numpy.where(zdata.mtype==4)
    qfidx = numpy.where(zdata.mtype==5)
    iidx = numpy.where(zdata.mtype==6)
    
    nvi = len(vidx[0])
    npi = len(pidx[0])
    nqi = len(qidx[0])
    npf = len(pfidx[0])
    nqf = len(qfidx[0])
    nii = len(iidx[0])
    
    busvi = zdata.mfrom[vidx]
    buspi = zdata.mfrom[pidx]
    busqi = zdata.mfrom[qidx]        
    buspf = zdata.mfrom[pfidx]
    busqf = zdata.mfrom[qfidx]
    branchpf = zdata.mbranch[pfidx]
    branchqf = zdata.mbranch[qfidx]
    branchii = zdata.mbranch[iidx]   

    z = zdata.mval
    idx = numpy.where(zdata.mbranch<0)
    z[idx] = -z[idx]
    
    Pinj = z[pidx]
    Qinj = z[qidx]
    Pbr = z[pfidx]
    Qbr = z[qfidx]
    
    idx = numpy.where(zdata.mstddev<10**(-6))
    zdata.mstddev[idx] = 10**(-6)
    weights = zdata.mstddev**(-2)
    W = numpy.diag(weights)
    
    znod = Znod(branch,node)
    
    """ Jacobian for Power Injection Measurements (converted to equivalent 
    rectangualar current measurements) """
    H2 = numpy.zeros((npi,2*branch.num+1))
    H3 = numpy.zeros((nqi,2*branch.num+1))
    for i in range(0,npi):
        m = buspi[i]
        to = numpy.where(branch.end==m)
        fr = numpy.where(branch.start==m)
        to1 = to[0]+1
        to2 = to1+branch.num
        fr1 = fr[0]+1
        fr2 = fr1+branch.num
        H2[i][to1] = 1
        H2[i][fr1] = -1
        H3[i][to2] = 1
        H3[i][fr2] = -1
        
    """ Jacobian for branch Power Measurements (converted to equivalent 
    rectangualar current measurements)"""
    H4 = numpy.zeros((npf,2*branch.num+1))
    H5 = numpy.zeros((nqf,2*branch.num+1))
    for i in range(0,npf):
        br = abs(branchpf[i])
        H4[i][br] = 1
        H5[i][br+branch.num] = 1
        
    epsilon = 5
    Vr = numpy.ones(node.num)
    Vx = numpy.zeros(node.num)
    V = Real_to_all(Vr, Vx)
    Ir = numpy.zeros(node.num)
    Ix = numpy.zeros(node.num)
    I = Real_to_all(Ir, Ix)
    num_iter = 0
    
    State = numpy.zeros(2*branch.num)
    State = numpy.concatenate((numpy.array([1,]),State),axis=0)
    
    while epsilon>10**(-6):
        """ Computation of equivalent current measurements in place of the power measurements """
        Irinj = (Pinj*V.real[buspi-1] + Qinj*V.imag[busqi-1])/(V.mag[buspi-1]**2)
        Ixinj = (Pinj*V.imag[buspi-1] - Qinj*V.real[busqi-1])/(V.mag[buspi-1]**2)
        z[pidx] = Irinj
        z[qidx] = Ixinj
        
        Irbr = (Pbr*V.real[buspf-1] + Qbr*V.imag[busqf-1])/(V.mag[buspf-1]**2)
        Ixbr = (Pbr*V.imag[buspf-1] - Qbr*V.real[busqf-1])/(V.mag[buspf-1]**2)
        z[pfidx] = Irbr
        z[qfidx] = Ixbr
        
        """ Voltage Magnitude Measurements """
        h1 = V.mag[busvi-1]
        H1 = numpy.zeros((nvi,2*branch.num+1))
        for i in range(0,nvi):
            m = busvi[i]-1
            if m == 0:
                H1[i][0] = 1
            else:
                idx = numpy.where(znod.Z[m]!=0)
                H1[i][0] = numpy.cos(V.phase[m])
                H1[i][idx[0]+1] = -znod.R[m][idx]*numpy.cos(V.phase[m]) - znod.X[m][idx]*numpy.sin(V.phase[m])
                H1[i][idx[0]+1+branch.num] = -znod.R[m][idx]*numpy.sin(V.phase[m]) + znod.X[m][idx]*numpy.cos(V.phase[m])
        
        """ Power Injection Measurements """
        h2 = numpy.inner(H2,State)
        h3 = numpy.inner(H3,State)
        
        """ Power Flow Measurements """
        h4 = numpy.inner(H4,State)
        h5 = numpy.inner(H5,State)
        
        """ Current Magnitude Measurements """
        h6 = I.mag[branchii]
        phi = I.phase[branchii]
        H6 = numpy.zeros((nii,2*branch.num+1))
        for i in range(0,nii):
            idx = branchii[i]
            H6[i][idx] = numpy.cos(phi[i])
            H6[i][idx+branch.num] = numpy.sin(phi[i])
            
        """ WLS computation """
        H = numpy.concatenate((H1,H2,H3,H4,H5,H6),axis=0)
        y = numpy.concatenate((h1,h2,h3,h4,h5,h6),axis=0)
        res = numpy.subtract(z,y)
        g = numpy.inner(H.transpose(),numpy.inner(W,res))
        WH = numpy.inner(W,H.transpose())
        G = numpy.inner(H.transpose(),WH.transpose())
        
        Ginv = numpy.linalg.inv(G)
        Delta_State = numpy.inner(Ginv,g)
        
        State = State + Delta_State
        epsilon = numpy.amax(numpy.absolute(Delta_State))
        
        V.real[0] = State[0]
        Ir = State[1:branch.num+1]
        Ix = State[branch.num+1:]
        Irx = Ir + 1j*Ix
        
        """ Forward sweep to calculate node voltage along the whole grid """
        Vcomplex = numpy.inner((V.real[0] + 1j*V.imag[0]),numpy.ones(node.num))
        DeltaV = numpy.inner(znod.Z,Irx)
        V.complex = Vcomplex[0] - DeltaV
        V.real = numpy.real(V.complex)
        V.imag = numpy.imag(V.complex)
        V = Real_to_all(V.real, V.imag)
        I = Real_to_all(Ir,Ix)
        
        num_iter = num_iter+1
    
    return V, I



def DssePmu(branch, node, zdata):
    """ It performs state estimation using branch current state variables and it is customized
    to work using only PMU measurements."""
    
    class Znod:
        def __init__(self,branch,node):
            self.Z = numpy.zeros((node.num-1,branch.num),dtype=numpy.complex_)
            for k in range(1,node.num):
                i = k-1
                aa = branch.start[i]-2
                if aa == -1:
                    self.Z[i] = self.Z[i]
                else:
                    self.Z[i] = self.Z[aa]
                self.Z[i][i] = branch.Z[i]
            zzz = numpy.zeros((1,branch.num))
            self.Z = numpy.concatenate((zzz,self.Z), axis=0)
            self.R = self.Z.real
            self.X = self.Z.imag
            
    znod = Znod(branch,node)
    
    

def DsseMixed(branch, node, zdata):
    """ It performs state estimation using branch current state variables and it is built
    to work in scenarios where both conventional and PMU measurements are simultaneously present."""
    
    class Znod:
        def __init__(self,branch,node):
            self.Z = numpy.zeros((node.num-1,branch.num),dtype=numpy.complex_)
            for k in range(1,node.num):
                i = k-1
                aa = branch.start[i]-2
                if aa == -1:
                    self.Z[i] = self.Z[i]
                else:
                    self.Z[i] = self.Z[aa]
                self.Z[i][i] = branch.Z[i]
            zzz = numpy.zeros((1,branch.num))
            self.Z = numpy.concatenate((zzz,self.Z), axis=0)
            self.R = self.Z.real
            self.X = self.Z.imag
            
            
    znod = Znod(branch,node)
    
    
    
    
    