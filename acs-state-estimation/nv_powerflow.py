import numpy
import math

class Real_to_all:
    def __init__(self,Ar, Ax):
        self.real = Ar
        self.imag = Ax
        self.complex = Ar + 1j*Ax
        self.mag = numpy.absolute(self.complex)
        self.phase = numpy.angle(self.complex)
        
def Ymatrix_calc(branch, node):
    Ymatrix = numpy.zeros((node.num,node.num),dtype=numpy.complex)
    Adjacencies = [[] for _ in range(node.num)]
    for index in range(branch.num):
        fr = branch.start[index] - 1
        to = branch.end[index] - 1
        Ymatrix[fr][to] = - branch.Y[index]
        Ymatrix[to][fr] = - branch.Y[index]
        Ymatrix[fr][fr] += branch.Y[index]
        Ymatrix[to][to] += branch.Y[index]
        Adjacencies[fr].append(to+1)
        Adjacencies[to].append(fr+1)
    return Ymatrix, Adjacencies
        
def NV_power_flow(branch, node):
    """It performs Power Flow by using rectangular node voltage state variables."""
    
    Ymatrix, Adj = Ymatrix_calc(branch,node)
    
    z = numpy.zeros(2*(node.num))
    h = numpy.zeros(2*(node.num))
    H = numpy.zeros((2*(node.num),2*(node.num)))
    
    for k in range(1,node.num+1):
        i = k-1
        m = 2*i
        i2 = i + node.num
        t = node.type[i]
        if t == 'slack':
            V = node.pwr_flow_values[1][i]
            theta = node.pwr_flow_values[2][i]
            z[m] = V*math.cos(theta)
            z[m+1] = V*math.sin(theta)
            H[m][i] = 1
            H[m+1][i2] = 1
        elif t == 'PQ':
            H[m][i] = - numpy.real(Ymatrix[i][i])
            H[m][i2] = numpy.imag(Ymatrix[i][i])
            H[m+1][i] = - numpy.imag(Ymatrix[i][i])
            H[m+1][i2] = - numpy.real(Ymatrix[i][i])
            idx1 = numpy.subtract(Adj[i],1)
            idx2 = idx1 + node.num
            H[m][idx1] = - numpy.real(Ymatrix[i][idx1])
            H[m][idx2] = numpy.imag(Ymatrix[i][idx1])
            H[m+1][idx1] = - numpy.imag(Ymatrix[i][idx1])
            H[m+1][idx2] = - numpy.real(Ymatrix[i][idx1])
        elif t == 'PV':
            z[m+1] = node.pwr_flow_values[2][i]
            H[m][i] = - numpy.real(Ymatrix[i][i])
            H[m][i2] = numpy.imag(Ymatrix[i][i])
            idx1 = numpy.subtract(Adj[i],1)
            idx2 = idx1 + node.num
            H[m][idx1] = - numpy.real(Ymatrix[i][idx1])
            H[m][idx2] = numpy.imag(Ymatrix[i][idx1])
    
    epsilon = 5
    Vr = numpy.ones(node.num)
    Vx = numpy.zeros(node.num)
    V = Real_to_all(Vr, Vx)
    num_iter = 0
    
    StateVr = numpy.ones(node.num)
    StateVx = numpy.zeros(node.num)
    State = numpy.concatenate((StateVr,StateVx),axis=0)
    
    while epsilon>10**(-10):
        for k in range(1,node.num+1):
            i = k-1
            m = 2*i
            i2 = i + node.num
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
                H[m+1][i] = numpy.cos(V.phase[i])
                H[m+1][i2] = numpy.sin(V.phase[i])
        
        r = numpy.subtract(z,h)
        Hinv = numpy.linalg.inv(H)
        Delta_State = numpy.inner(Hinv,r)
        State = State + Delta_State
        epsilon = numpy.amax(numpy.absolute(Delta_State))
        
        V.real = State[:node.num]
        V.imag = State[node.num:]
        V = Real_to_all(V.real, V.imag)
                
        num_iter = num_iter+1
        
    Irx = numpy.zeros((branch.num),dtype=numpy.complex)
    for idx in range(branch.num):
        fr = branch.start[idx]-1
        to = branch.end[idx]-1
        Irx[idx] = - (V.complex[fr] - V.complex[to])*Ymatrix[fr][to]
    Ir = numpy.real(Irx)
    Ix = numpy.imag(Irx)
    
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
    S2_rx = - numpy.multiply(V.complex[branch.end-1],numpy.conj(I.complex))
    S1 = Real_to_all(numpy.real(S1_rx),numpy.imag(S1_rx))
    S2 = Real_to_all(numpy.real(S2_rx),numpy.imag(S2_rx))
    
    return V, I, Iinj, S1, S2, Sinj, num_iter
    