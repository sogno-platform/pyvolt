def BC_power_flow(branch, node):
    """It performs Power Flow by using branch current state variables."""
    import numpy
    import math
    
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
            
    class Real_to_all:
        def __init__(self,Ar, Ax):
            self.real = Ar
            self.imag = Ax
            self.complex = Ar + 1j*Ax
            self.mag = numpy.absolute(self.complex)
            self.phase = numpy.angle(self.complex)
            
    
    znod = Znod(branch,node)
    z = numpy.zeros(2*(branch.num+1))
    h = numpy.zeros(2*(branch.num+1))
    H = numpy.zeros((2*(branch.num+1),2*(branch.num+1)))
    
    for k in range(1,node.num+1):
        i = k-1
        m = 2*i
        t = node.type[i]
        if t == 'slack':
            V = node.pwr_flow_values[1][i]
            theta = node.pwr_flow_values[2][i]
            z[m] = V*math.cos(theta)
            z[m+1] = V*math.sin(theta)
            H[m][0] = 1;
            H[m][2:] = numpy.concatenate((-znod.R[i],znod.X[i]),axis=0)
            H[m+1][1] = 1
            H[m+1][2:] = numpy.concatenate((-znod.X[i],-znod.R[i]),axis=0)
        elif t == 'PQ':
            to = numpy.where(branch.end==k)
            fr = numpy.where(branch.start==k)
            to1 = to[0]+2
            to2 = to1+branch.num
            fr1 = fr[0]+2
            fr2 = fr1+branch.num
            H[m][to1] = 1
            H[m+1][to2] = 1
            H[m][fr1] = -1
            H[m+1][fr2] = -1
        elif t == 'PV':
            z[m+1] = node.pwr_flow_values[2][i]
            to = numpy.where(branch.end==k)
            fr = numpy.where(branch.start==k)
            to1 = to[0]+2
            fr1 = fr[0]+2
            H[m][to1] = 1
            H[m][fr1] = -1
    
    epsilon = 5
    Vr = numpy.ones(node.num)
    Vx = numpy.zeros(node.num)
    V = Real_to_all(Vr, Vx)
    num_iter = 0
    
    State = numpy.zeros(2*branch.num)
    State = numpy.concatenate((numpy.array([1,0]),State),axis=0)
    
    while epsilon>10**(-10):
        for k in range(1,node.num+1):
            i = k-1
            m = 2*i
            t = node.type[i]
            if t == 'slack':
                h[m] = numpy.inner(H[m],State)
                h[m+1] = numpy.inner(H[m+1],State)
            elif t == 'PQ':
                z[m] = (node.pwr_flow_values[1][i]*V.real[i] + node.pwr_flow_values[2][i]*V.imag[i])/(V.mag[i]**2)
                z[m+1] = (node.pwr_flow_values[1][i]*V.imag[i] - node.pwr_flow_values[2][i]*V.real[i])/(V.mag[i]**2)
                h[m] = numpy.inner(H[m],State)
                h[m+1] = numpy.inner(H[m+1],State)
            elif t == 'PV':
                z[m] = (node.pwr_flow_values[1][i]*V.real[i] + node.pwr_flow_values[2][i]*V.imag[i])/(V.mag[i]**2)
                h[m] = numpy.inner(H[m],State)
                h[m+1] = V.mag[i]
                H[m+1][0] = numpy.cos(V.theta[i])
                H[m+1][1] = numpy.sin(V.theta[i])
                idx = numpy.where(znod.Z[i] != 0)
                idx1 = idx[0]+2
                idx2 = idx1[0]+branch.num
                H[m+1][idx1] = - znod.R[i][idx]*numpy.cos(V.theta[i]) - znod.X[i][idx]*numpy.sin(V.theta[i])
                H[m+1][idx2] = - znod.R[i][idx]*numpy.sin(V.theta[i]) - znod.X[i][idx]*numpy.cos(V.theta[i])
        
        r = numpy.subtract(z,h)
        Hinv = numpy.linalg.inv(H)
        Delta_State = numpy.inner(Hinv,r)
        State = State + Delta_State
        epsilon = numpy.amax(numpy.absolute(Delta_State))
        
        V.real[0] = State[0]
        V.imag[0] = State[1]
        Ir = State[2:branch.num+2]
        Ix = State[branch.num+2:]
        Irx = Ir + 1j*Ix
        
        Vcomplex = numpy.inner((V.real[0] + 1j*V.imag[0]),numpy.ones(node.num))
        DeltaV = numpy.inner(znod.Z,Irx)
        V.complex = Vcomplex[0] - DeltaV
        V.real = numpy.real(V.complex)
        V.imag = numpy.imag(V.complex)
        V = Real_to_all(V.real, V.imag)
        
        num_iter = num_iter+1
        
    I = Real_to_all(Ir,Ix)
    Iinj_r = numpy.zeros(node.num)
    Iinj_x = numpy.zeros(node.num)
    for k in range(1,node.num+1):
        to = numpy.where(branch.end==k)
        fr = numpy.where(branch.start==k)
        Iinj_r[k-1] = numpy.sum(I.real[to[0]]) - numpy.sum(I.real[fr[0]])
        Iinj_x[k-1] = numpy.sum(I.imag[to[0]]) - numpy.sum(I.imag[fr[0]])
        
    Iinj = Real_to_all(Iinj_r,Iinj_x)
    Sinj_rx = numpy.multiply(V.complex,numpy.conj(Iinj.complex))
    Sinj = Real_to_all(numpy.real(Sinj_rx),numpy.imag(Sinj_rx))
    S1_rx = numpy.multiply(V.complex[branch.start-1],numpy.conj(I.complex))
    S2_rx = numpy.multiply(V.complex[branch.end-1],numpy.conj(I.complex))
    S1 = Real_to_all(numpy.real(S1_rx),numpy.imag(S1_rx))
    S2 = Real_to_all(numpy.real(S2_rx),numpy.imag(S2_rx))
    
    return V, I, Iinj, S1, S2, Sinj, num_iter