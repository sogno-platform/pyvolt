import numpy

class Real_to_all:
    def __init__(self,Ar,Ax):
        """ Converts real and imaginary parts to all the other formats. """
        self.real = Ar
        self.imag = Ax
        self.complex = Ar + 1j*Ax
        self.mag = numpy.absolute(self.complex)
        self.phase = numpy.angle(self.complex)
        
class Complex_to_all:
    def __init__(self,Arx):
        """ Converts complex numbers to all the other formats. """
        self.complex = Arx
        self.real = numpy.real(self.complex)
        self.imag = numpy.imag(self.complex)      
        self.mag = numpy.absolute(self.complex)
        self.phase = numpy.angle(self.complex)
            
def DsseCall(branch, node, zdata, Ymatrix, Adj):
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
        Vest = DsseTrad(branch, node, zdata, Ymatrix, Adj)
    elif est_code == 2:
        Vest = DssePmu(branch, node, zdata, Ymatrix, Adj)
    else:
        Vest = DsseMixed(branch, node, zdata, Ymatrix, Adj)
        
    """ From here on, all the other quantities of the grid are calculated """
    Irx = numpy.zeros((branch.num),dtype=numpy.complex)
    for idx in range(branch.num):
        fr = branch.start[idx]-1
        to = branch.end[idx]-1
        Irx[idx] = - (Vest.complex[fr] - Vest.complex[to])*Ymatrix[fr][to]
    Ir = numpy.real(Irx)
    Ix = numpy.imag(Irx)
    
    Iest = Real_to_all(Ir,Ix)
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


def DsseTrad(branch, node, zdata, Ymatrix, Adj):
    """ It performs state estimation using rectangular node voltage state variables 
    and it is customized to work without PMU measurements"""
    
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
    fbuspf = zdata.mfrom[pfidx]
    tbuspf = zdata.mto[pfidx]
    fbusqf = zdata.mfrom[qfidx]
    fbusiamp = zdata.mfrom[iidx]
    tbusiamp = zdata.mto[iidx]
       
    z = zdata.mval
    
    Pinj = z[pidx]
    Qinj = z[qidx]
    Pbr = z[pfidx]
    Qbr = z[qfidx]
    
    idx = numpy.where(zdata.mstddev<10**(-6))
    zdata.mstddev[idx] = 10**(-6)
    weights = zdata.mstddev**(-2)
    W = numpy.diag(weights)
    
    Admittance = Complex_to_all(Ymatrix)
    Gmatrix = Admittance.real
    Bmatrix = Admittance.imag
    Yabs_matrix = Admittance.mag
    Yphase_matrix = Admittance.phase
        
    """ Jacobian for Power Injection Measurements (converted to equivalent 
    rectangualar current measurements) """
    H2 = numpy.zeros((npi,2*node.num-1))
    H3 = numpy.zeros((nqi,2*node.num-1))
    for i in range(npi):
        m = buspi[i]-1
        m2 = m + node.num - 1
        H2[i][m] = - Gmatrix[m][m]
        H2[i][m2] = Bmatrix[m][m]
        H3[i][m] = - Bmatrix[m][m]
        H3[i][m2] = - Gmatrix[m][m]
        idx = numpy.subtract(Adj[m],1)
        H2[i][idx] = - Gmatrix[m][idx]
        H3[i][idx] = - Bmatrix[m][idx]
        if 0 in idx:
            pos = numpy.where(idx==0)
            idx = numpy.delete(idx,pos)
        idx2 = idx + node.num-1
        H2[i][idx2] = Bmatrix[m][idx]
        H3[i][idx2] = - Gmatrix[m][idx]
        
    """ Jacobian for branch Power Measurements (converted to equivalent 
    rectangualar current measurements)"""
    H4 = numpy.zeros((npf,2*node.num-1))
    H5 = numpy.zeros((nqf,2*node.num-1))
    for i in range(npf):
        m = fbuspf[i]-1
        n = tbuspf[i]-1
        H4[i][m] = - Gmatrix[m][n]
        H4[i][n] = Gmatrix[m][n]
        H5[i][m] = - Bmatrix[m][n]
        H5[i][n] = Bmatrix[m][n]
        if m > 0:
            m2 = m + node.num-1
            H4[i][m2] = Bmatrix[m][n]
            H5[i][m2] = - Gmatrix[m][n]
        if n > 0:
            n2 = n + node.num-1
            H4[i][n2] = - Bmatrix[m][n]
            H5[i][n2] = Gmatrix[m][n]
        
    epsilon = 5
    Vr = numpy.ones(node.num)
    Vx = numpy.zeros(node.num)
    V = Real_to_all(Vr, Vx)
    num_iter = 0
    
    StateVr = numpy.ones(node.num)
    StateVx = numpy.zeros(node.num-1)
    State = numpy.concatenate((StateVr,StateVx),axis=0)
    
    while epsilon>10**(-6):
        """ Computation of equivalent current measurements in place of the power measurements """
        Irinj = (Pinj*V.real[buspi-1] + Qinj*V.imag[busqi-1])/(V.mag[buspi-1]**2)
        Ixinj = (Pinj*V.imag[buspi-1] - Qinj*V.real[busqi-1])/(V.mag[buspi-1]**2)
        z[pidx] = Irinj
        z[qidx] = Ixinj
        
        Irbr = (Pbr*V.real[fbuspf-1] + Qbr*V.imag[fbusqf-1])/(V.mag[fbuspf-1]**2)
        Ixbr = (Pbr*V.imag[fbuspf-1] - Qbr*V.real[fbusqf-1])/(V.mag[fbuspf-1]**2)
        z[pfidx] = Irbr
        z[qfidx] = Ixbr
        
        """ Voltage Magnitude Measurements """
        h1 = V.mag[busvi-1]
        H1 = numpy.zeros((nvi,2*node.num-1))
        for i in range(nvi):
            m = busvi[i]-1
            H1[i][m] = numpy.cos(V.phase[m])
            if m > 0:
                m2 = m + node.num-1
                H1[i][m2] = numpy.sin(V.phase[m])
                
        """ Power Injection Measurements """
        h2 = numpy.inner(H2,State)
        h3 = numpy.inner(H3,State)
        
        """ Power Flow Measurements """
        h4 = numpy.inner(H4,State)
        h5 = numpy.inner(H5,State)
        
        """ Current Magnitude Measurements """
        h6re = numpy.zeros((nii))
        h6im = numpy.zeros((nii))
        h6complex = numpy.zeros((nii),dtype=complex)
        h6 = numpy.ones((nii))
        H6 = numpy.zeros((nii,2*node.num-1))
        for i in range(nii):
            m = fbusiamp[i]-1
            n = tbusiamp[i]-1
            h6re[i] = Yabs_matrix[m][n]*((V.real[n]-V.real[m])*numpy.cos(Yphase_matrix[m][n]) + (V.imag[m]-V.imag[n])*numpy.sin(Yphase_matrix[m][n]))
            h6im[i] = Yabs_matrix[m][n]*((V.real[n]-V.real[m])*numpy.sin(Yphase_matrix[m][n]) + (V.imag[n]-V.imag[m])*numpy.cos(Yphase_matrix[m][n]))
            h6complex[i] = h6re[i] + 1j*h6im[i]
            if num_iter>0:
                h6[i] = numpy.absolute(h6complex[i])
            H6[i][m] = - Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6re[i] + numpy.sin(Yphase_matrix[m][n])*h6im[i])/h6[i]
            H6[i][n] = Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6re[i] + numpy.sin(Yphase_matrix[m][n])*h6im[i])/h6[i]
            if m > 0:
                m2 = m + node.num-1
                H6[i][m2] = - Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6im[i] - numpy.sin(Yphase_matrix[m][n])*h6re[i])/h6[i]
            if n > 0:
                n2 = n + node.num-1
                H6[i][n2] = Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6im[i] - numpy.sin(Yphase_matrix[m][n])*h6re[i])/h6[i]
        
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
        
        V.real = State[:node.num]
        V.imag = State[node.num:]
        V.imag = numpy.concatenate(([0],V.imag), axis=0)
        V = Real_to_all(V.real, V.imag)
        
        num_iter = num_iter+1
    
    return V



def DssePmu(branch, node, zdata, Ymatrix, Adj):
    """ It performs state estimation using rectangular node voltage state variables 
    and it is customized to work using PMU measurements."""
    
    pidx = numpy.where(zdata.mtype==2)
    qidx = numpy.where(zdata.mtype==3)
    pfidx = numpy.where(zdata.mtype==4)
    qfidx = numpy.where(zdata.mtype==5)
    vmagpmuidx = numpy.where(zdata.mtype==7)
    vphasepmuidx = numpy.where(zdata.mtype==8)
    imagpmuidx = numpy.where(zdata.mtype==9)
    iphasepmuidx = numpy.where(zdata.mtype==10)
    
    npi = len(pidx[0])
    nqi = len(qidx[0])
    npf = len(pfidx[0])
    nqf = len(qfidx[0])
    nvpmu = len(vmagpmuidx[0])
    nipmu = len(imagpmuidx[0])
    
    buspi = zdata.mfrom[pidx]
    busqi = zdata.mfrom[qidx]
    fbuspf = zdata.mfrom[pfidx]
    tbuspf = zdata.mto[pfidx]
    fbusqf = zdata.mfrom[qfidx]
    busvpmu = zdata.mfrom[vmagpmuidx]
    fbusipmu = zdata.mfrom[imagpmuidx]
    tbusipmu = zdata.mto[imagpmuidx]
       
    z = zdata.mval
    
    Pinj = z[pidx]
    Qinj = z[qidx]
    Pbr = z[pfidx]
    Qbr = z[qfidx]
    
    idx = numpy.where(zdata.mstddev<10**(-6))
    zdata.mstddev[idx] = 10**(-6)
    weights = zdata.mstddev**(-2)
    W = numpy.diag(weights)
    
    Admittance = Complex_to_all(Ymatrix)
    Gmatrix = Admittance.real
    Bmatrix = Admittance.imag
        
    """ Jacobian for Power Injection Measurements (converted to equivalent 
    rectangualar current measurements) """
    H2 = numpy.zeros((npi,2*node.num))
    H3 = numpy.zeros((nqi,2*node.num))
    for i in range(npi):
        m = buspi[i] - 1
        m2 = m + node.num
        H2[i][m] = - Gmatrix[m][m]
        H2[i][m2] = Bmatrix[m][m]
        H3[i][m] = - Bmatrix[m][m]
        H3[i][m2] = - Gmatrix[m][m]
        idx = numpy.subtract(Adj[m],1)
        H2[i][idx] = - Gmatrix[m][idx]
        H3[i][idx] = - Bmatrix[m][idx]
        idx2 = idx + node.num
        H2[i][idx2] = Bmatrix[m][idx]
        H3[i][idx2] = - Gmatrix[m][idx]
        
    """ Jacobian for branch Power Measurements (converted to equivalent 
    rectangualar current measurements)"""
    H4 = numpy.zeros((npf,2*node.num))
    H5 = numpy.zeros((nqf,2*node.num))
    for i in range(npf):
        m = fbuspf[i]-1
        n = tbuspf[i]-1
        H4[i][m] = - Gmatrix[m][n]
        H4[i][n] = Gmatrix[m][n]
        H5[i][m] = - Bmatrix[m][n]
        H5[i][n] = Bmatrix[m][n]
        m2 = m + node.num
        H4[i][m2] = Bmatrix[m][n]
        H5[i][m2] = - Gmatrix[m][n]
        n2 = n + node.num
        H4[i][n2] = - Bmatrix[m][n]
        H5[i][n2] = Gmatrix[m][n]
        
    """ Jacobian for Voltage Pmu Measurements (converted into rectangular) """
    H7 = numpy.zeros((nvpmu,2*node.num))
    H8 = numpy.zeros((nvpmu,2*node.num))
    for i in range(nvpmu):
        idx1 = vmagpmuidx[0][i]
        idx2 = vphasepmuidx[0][i]
        vamp = z[idx1]
        vtheta = z[idx2]
        z[idx1] = vamp*numpy.cos(vtheta)
        z[idx2] = vamp*numpy.sin(vtheta)
        rot_mat = numpy.array([[numpy.cos(vtheta), - vamp*numpy.sin(vtheta)], [numpy.sin(vtheta), vamp*numpy.cos(vtheta)]])
        starting_cov = numpy.array([[weights[idx1], 0], [0, weights[idx2]]])
        final_cov = numpy.inner(rot_mat,numpy.inner(starting_cov,rot_mat.transpose()))
        W[idx1][idx1] = final_cov[0][0]
        W[idx2][idx2] = final_cov[1][1]
        W[idx1][idx2] = final_cov[0][1]
        W[idx2][idx1] = final_cov[1][0]
        m = busvpmu[i]-1
        H7[i][m] = 1
        m2 = m + node.num
        H8[i][m2] = 1
        
    """ Jacobian for Current Pmu Measurements (converted into rectangular) """
    H9 = numpy.zeros((nipmu,2*node.num))
    H10 = numpy.zeros((nipmu,2*node.num))
    for i in range(nipmu):
        idx1 = imagpmuidx[0][i]
        idx2 = iphasepmuidx[0][i]
        iamp = z[idx1]
        itheta = z[idx2]
        z[idx1] = iamp*numpy.cos(itheta)
        z[idx2] = iamp*numpy.sin(itheta)
        rot_mat = numpy.array([[numpy.cos(itheta), - iamp*numpy.sin(itheta)], [numpy.sin(itheta), iamp*numpy.cos(itheta)]])
        starting_cov = numpy.array([[weights[idx1], 0], [0, weights[idx2]]])
        final_cov = numpy.inner(rot_mat,numpy.inner(starting_cov,rot_mat.transpose()))
        W[idx1][idx1] = final_cov[0][0]
        W[idx2][idx2] = final_cov[1][1]
        W[idx1][idx2] = final_cov[0][1]
        W[idx2][idx1] = final_cov[1][0]
        m = fbusipmu[i]-1
        n = tbusipmu[i]-1
        H9[i][m] = - Gmatrix[m][n]
        H9[i][n] = Gmatrix[m][n]
        H10[i][m] = - Bmatrix[m][n]
        H10[i][n] = Bmatrix[m][n]
        m2 = m + node.num
        n2 = n + node.num
        H9[i][m2] = Bmatrix[m][n]
        H9[i][n2] = - Bmatrix[m][n]
        H10[i][m2] = - Gmatrix[m][n]
        H10[i][n2] = Gmatrix[m][n]
        
    epsilon = 5
    Vr = numpy.ones(node.num)
    Vx = numpy.zeros(node.num)
    V = Real_to_all(Vr, Vx)
    num_iter = 0
    
    StateVr = numpy.ones(node.num)
    StateVx = numpy.zeros(node.num)
    State = numpy.concatenate((StateVr,StateVx),axis=0)
    
    H = numpy.concatenate((H2,H3,H4,H5,H7,H8,H9,H10),axis=0)
    WH = numpy.inner(W,H.transpose())
    G = numpy.inner(H.transpose(),WH.transpose())
    Ginv = numpy.linalg.inv(G)
    
    
    while epsilon>10**(-6):
        """ Computation of equivalent current measurements in place of the power measurements """
        Irinj = (Pinj*V.real[buspi-1] + Qinj*V.imag[busqi-1])/(V.mag[buspi-1]**2)
        Ixinj = (Pinj*V.imag[buspi-1] - Qinj*V.real[busqi-1])/(V.mag[buspi-1]**2)
        z[pidx] = Irinj
        z[qidx] = Ixinj
        
        Irbr = (Pbr*V.real[fbuspf-1] + Qbr*V.imag[fbusqf-1])/(V.mag[fbuspf-1]**2)
        Ixbr = (Pbr*V.imag[fbuspf-1] - Qbr*V.real[fbusqf-1])/(V.mag[fbuspf-1]**2)
        z[pfidx] = Irbr
        z[qfidx] = Ixbr
        
        
        """ WLS computation """
        y = numpy.inner(H,State)
        res = numpy.subtract(z,y)
        g = numpy.inner(H.transpose(),numpy.inner(W,res))

        Delta_State = numpy.inner(Ginv,g)
        
        State = State + Delta_State
        epsilon = numpy.amax(numpy.absolute(Delta_State))
        
        V.real = State[:node.num]
        V.imag = State[node.num:]
        V = Real_to_all(V.real, V.imag)
        
        num_iter = num_iter+1
        
    return V
    
    

def DsseMixed(branch, node, zdata, Ymatrix, Adj):
    """ It performs state estimation using rectangular node voltage state variables
    and it is built to work in scenarios where both conventional and PMU measurements 
    are simultaneously present."""
    
    vidx = numpy.where(zdata.mtype==1)
    pidx = numpy.where(zdata.mtype==2)
    qidx = numpy.where(zdata.mtype==3)
    pfidx = numpy.where(zdata.mtype==4)
    qfidx = numpy.where(zdata.mtype==5)
    iidx = numpy.where(zdata.mtype==6)
    vmagpmuidx = numpy.where(zdata.mtype==7)
    vphasepmuidx = numpy.where(zdata.mtype==8)
    imagpmuidx = numpy.where(zdata.mtype==9)
    iphasepmuidx = numpy.where(zdata.mtype==10)
    
    nvi = len(vidx[0])
    npi = len(pidx[0])
    nqi = len(qidx[0])
    npf = len(pfidx[0])
    nqf = len(qfidx[0])
    nii = len(iidx[0])
    nvpmu = len(vmagpmuidx[0])
    nipmu = len(imagpmuidx[0])
    
    busvi = zdata.mfrom[vidx]
    buspi = zdata.mfrom[pidx]
    busqi = zdata.mfrom[qidx]
    fbuspf = zdata.mfrom[pfidx]
    tbuspf = zdata.mto[pfidx]
    fbusqf = zdata.mfrom[qfidx]
    fbusiamp = zdata.mfrom[iidx]
    tbusiamp = zdata.mto[iidx]
    busvpmu = zdata.mfrom[vmagpmuidx]
    fbusipmu = zdata.mfrom[imagpmuidx]
    tbusipmu = zdata.mto[imagpmuidx]
       
    z = zdata.mval
    
    Pinj = z[pidx]
    Qinj = z[qidx]
    Pbr = z[pfidx]
    Qbr = z[qfidx]
    
    idx = numpy.where(zdata.mstddev<10**(-6))
    zdata.mstddev[idx] = 10**(-6)
    weights = zdata.mstddev**(-2)
    W = numpy.diag(weights)
    
    Admittance = Complex_to_all(Ymatrix)
    Gmatrix = Admittance.real
    Bmatrix = Admittance.imag
    Yabs_matrix = Admittance.mag
    Yphase_matrix = Admittance.phase
        
    """ Jacobian for Power Injection Measurements (converted to equivalent 
    rectangualar current measurements) """
    H2 = numpy.zeros((npi,2*node.num))
    H3 = numpy.zeros((nqi,2*node.num))
    for i in range(npi):
        m = buspi[i] - 1
        m2 = m + node.num
        H2[i][m] = - Gmatrix[m][m]
        H2[i][m2] = Bmatrix[m][m]
        H3[i][m] = - Bmatrix[m][m]
        H3[i][m2] = - Gmatrix[m][m]
        idx = numpy.subtract(Adj[m],1)
        H2[i][idx] = - Gmatrix[m][idx]
        H3[i][idx] = - Bmatrix[m][idx]
        idx2 = idx + node.num
        H2[i][idx2] = Bmatrix[m][idx]
        H3[i][idx2] = - Gmatrix[m][idx]
        
    """ Jacobian for branch Power Measurements (converted to equivalent 
    rectangualar current measurements)"""
    H4 = numpy.zeros((npf,2*node.num))
    H5 = numpy.zeros((nqf,2*node.num))
    for i in range(npf):
        m = fbuspf[i]-1
        n = tbuspf[i]-1
        H4[i][m] = - Gmatrix[m][n]
        H4[i][n] = Gmatrix[m][n]
        H5[i][m] = - Bmatrix[m][n]
        H5[i][n] = Bmatrix[m][n]
        m2 = m + node.num
        H4[i][m2] = Bmatrix[m][n]
        H5[i][m2] = - Gmatrix[m][n]
        n2 = n + node.num
        H4[i][n2] = - Bmatrix[m][n]
        H5[i][n2] = Gmatrix[m][n]
        
    """ Jacobian for Voltage Pmu Measurements (converted into rectangular) """
    H7 = numpy.zeros((nvpmu,2*node.num))
    H8 = numpy.zeros((nvpmu,2*node.num))
    for i in range(nvpmu):
        idx1 = vmagpmuidx[0][i]
        idx2 = vphasepmuidx[0][i]
        vamp = z[idx1]
        vtheta = z[idx2]
        z[idx1] = vamp*numpy.cos(vtheta)
        z[idx2] = vamp*numpy.sin(vtheta)
        rot_mat = numpy.array([[numpy.cos(vtheta), - vamp*numpy.sin(vtheta)], [numpy.sin(vtheta), vamp*numpy.cos(vtheta)]])
        starting_cov = numpy.array([[weights[idx1], 0], [0, weights[idx2]]])
        final_cov = numpy.inner(rot_mat,numpy.inner(starting_cov,rot_mat.transpose()))
        W[idx1][idx1] = final_cov[0][0]
        W[idx2][idx2] = final_cov[1][1]
        W[idx1][idx2] = final_cov[0][1]
        W[idx2][idx1] = final_cov[1][0]
        m = busvpmu[i]-1
        H7[i][m] = 1
        m2 = m + node.num
        H8[i][m2] = 1
        
    """ Jacobian for Current Pmu Measurements (converted into rectangular) """
    H9 = numpy.zeros((nipmu,2*node.num))
    H10 = numpy.zeros((nipmu,2*node.num))
    for i in range(nipmu):
        idx1 = imagpmuidx[0][i]
        idx2 = iphasepmuidx[0][i]
        iamp = z[idx1]
        itheta = z[idx2]
        z[idx1] = iamp*numpy.cos(itheta)
        z[idx2] = iamp*numpy.sin(itheta)
        rot_mat = numpy.array([[numpy.cos(itheta), - iamp*numpy.sin(itheta)], [numpy.sin(itheta), iamp*numpy.cos(itheta)]])
        starting_cov = numpy.array([[weights[idx1], 0], [0, weights[idx2]]])
        final_cov = numpy.inner(rot_mat,numpy.inner(starting_cov,rot_mat.transpose()))
        W[idx1][idx1] = final_cov[0][0]
        W[idx2][idx2] = final_cov[1][1]
        W[idx1][idx2] = final_cov[0][1]
        W[idx2][idx1] = final_cov[1][0]
        m = fbusipmu[i]-1
        n = tbusipmu[i]-1
        H9[i][m] = - Gmatrix[m][n]
        H9[i][n] = Gmatrix[m][n]
        H10[i][m] = - Bmatrix[m][n]
        H10[i][n] = Bmatrix[m][n]
        m2 = m + node.num
        n2 = n + node.num
        H9[i][m2] = Bmatrix[m][n]
        H9[i][n2] = - Bmatrix[m][n]
        H10[i][m2] = - Gmatrix[m][n]
        H10[i][n2] = Gmatrix[m][n]
        
    epsilon = 5
    Vr = numpy.ones(node.num)
    Vx = numpy.zeros(node.num)
    V = Real_to_all(Vr, Vx)
    num_iter = 0
    
    StateVr = numpy.ones(node.num)
    StateVx = numpy.zeros(node.num)
    State = numpy.concatenate((StateVr,StateVx),axis=0)
    
    while epsilon>10**(-6):
        """ Computation of equivalent current measurements in place of the power measurements """
        Irinj = (Pinj*V.real[buspi-1] + Qinj*V.imag[busqi-1])/(V.mag[buspi-1]**2)
        Ixinj = (Pinj*V.imag[buspi-1] - Qinj*V.real[busqi-1])/(V.mag[buspi-1]**2)
        z[pidx] = Irinj
        z[qidx] = Ixinj
        
        Irbr = (Pbr*V.real[fbuspf-1] + Qbr*V.imag[fbusqf-1])/(V.mag[fbuspf-1]**2)
        Ixbr = (Pbr*V.imag[fbuspf-1] - Qbr*V.real[fbusqf-1])/(V.mag[fbuspf-1]**2)
        z[pfidx] = Irbr
        z[qfidx] = Ixbr
        
        """ Voltage Magnitude Measurements """
        h1 = V.mag[busvi-1]
        H1 = numpy.zeros((nvi,2*node.num))
        for i in range(nvi):
            m = busvi[i]-1
            H1[i][m] = numpy.cos(V.phase[m])
            m2 = m + node.num
            H1[i][m2] = numpy.sin(V.phase[m])
                
        """ Power Injection Measurements """
        h2 = numpy.inner(H2,State)
        h3 = numpy.inner(H3,State)
        
        """ Power Flow Measurements """
        h4 = numpy.inner(H4,State)
        h5 = numpy.inner(H5,State)
        
        """ Current Magnitude Measurements """
        h6re = numpy.zeros((nii))
        h6im = numpy.zeros((nii))
        h6complex = numpy.zeros((nii),dtype=complex)
        h6 = numpy.ones((nii))
        H6 = numpy.zeros((nii,2*node.num))
        for i in range(nii):
            m = fbusiamp[i]-1
            n = tbusiamp[i]-1
            h6re[i] = Yabs_matrix[m][n]*((V.real[n]-V.real[m])*numpy.cos(Yphase_matrix[m][n]) + (V.imag[m]-V.imag[n])*numpy.sin(Yphase_matrix[m][n]))
            h6im[i] = Yabs_matrix[m][n]*((V.real[n]-V.real[m])*numpy.sin(Yphase_matrix[m][n]) + (V.imag[n]-V.imag[m])*numpy.cos(Yphase_matrix[m][n]))
            h6complex[i] = h6re[i] + 1j*h6im[i]
            if num_iter>0:
                h6[i] = numpy.absolute(h6complex[i])
            H6[i][m] = - Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6re[i] + numpy.sin(Yphase_matrix[m][n])*h6im[i])/h6[i]
            H6[i][n] = Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6re[i] + numpy.sin(Yphase_matrix[m][n])*h6im[i])/h6[i]
            m2 = m + node.num
            H6[i][m2] = - Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6im[i] - numpy.sin(Yphase_matrix[m][n])*h6re[i])/h6[i]
            n2 = n + node.num
            H6[i][n2] = Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6im[i] - numpy.sin(Yphase_matrix[m][n])*h6re[i])/h6[i]
        
        """ PMU Voltage Measurements """
        h7 = numpy.inner(H7,State)
        h8 = numpy.inner(H8,State)
        
        """ PMU Current Measurements """
        h9 = numpy.inner(H9,State)
        h10 = numpy.inner(H10,State)
        
        """ WLS computation """
        H = numpy.concatenate((H1,H2,H3,H4,H5,H6,H7,H8,H9,H10),axis=0)
        y = numpy.concatenate((h1,h2,h3,h4,h5,h6,h7,h8,h9,h10),axis=0)
        res = numpy.subtract(z,y)
        g = numpy.inner(H.transpose(),numpy.inner(W,res))
        WH = numpy.inner(W,H.transpose())
        G = numpy.inner(H.transpose(),WH.transpose())
        
        Ginv = numpy.linalg.inv(G)
        Delta_State = numpy.inner(Ginv,g)
        
        State = State + Delta_State
        epsilon = numpy.amax(numpy.absolute(Delta_State))
        
        V.real = State[:node.num]
        V.imag = State[node.num:]
        V = Real_to_all(V.real, V.imag)
        
        num_iter = num_iter+1
    
    return V
    
    
    
    