import numpy as np
import math

class Znod:
    def __init__(self, branches, nodes):
        # Complex impedance matrix
        self.Z = np.zeros((nodes.num-1,branches.num), dtype=np.complex)
        for k in range(1,nodes.num):
            i = k-1
            aa = branches.start[i]-2
            if aa == -1:
                self.Z[i] = self.Z[i]
            else:
                self.Z[i] = self.Z[aa]
            self.Z[i][i] = branches.Z[i]

        # Add one row of zeros
        zzz = np.zeros((1,branches.num))
        self.Z = np.concatenate((zzz,self.Z), axis=0)
        self.R = self.Z.real
        self.X = self.Z.imag

def BC_power_flow(branches, nodes):
    """It performs Power Flow by using rectangular branches current state variables."""

    znod = Znod(branches,nodes)
    # real matrices storing complex numbers 
    # real in position n and imag in position n+1
    z = np.zeros(2*(branches.num+1))
    h = np.zeros(2*(branches.num+1))
    H = np.zeros((2*(branches.num+1),2*(branches.num+1)))
    
    for k in range(1,nodes.num+1):
        i = k-1
        m = 2*i
        if nodes.type[i] == 'slack':
            V = nodes.pwr_flow_values[1][i]
            theta = nodes.pwr_flow_values[2][i]
            z[m] = V*math.cos(theta)
            z[m+1] = V*math.sin(theta)
            H[m][0] = 1;
            H[m][2:] = np.concatenate((-znod.R[i],znod.X[i]),axis=0)
            H[m+1][1] = 1
            H[m+1][2:] = np.concatenate((-znod.X[i],-znod.R[i]),axis=0)
        elif nodes.type[i] == 'PQ':
            to = np.where(branches.end==k)
            fr = np.where(branches.start==k)
            to1 = to[0]+2
            to2 = to1+branches.num
            fr1 = fr[0]+2
            fr2 = fr1+branches.num
            H[m][to1] = 1
            H[m+1][to2] = 1
            H[m][fr1] = -1
            H[m+1][fr2] = -1
        elif nodes.type[i] == 'PV':
            z[m+1] = nodes.pwr_flow_values[2][i]
            to = np.where(branches.end==k)
            fr = np.where(branches.start==k)
            to1 = to[0]+2
            fr1 = fr[0]+2
            H[m][to1] = 1
            H[m][fr1] = -1
    
    diff = 5
    epsilon = 10**(-10)
    V = np.ones(nodes.num) + 1j* np.zeros(nodes.num)
    num_iter = 0
    
    State = np.ones(2*branches.num)
    State = np.concatenate((np.array([1,0]),State),axis=0)
    
    while diff > epsilon:
        for k in range(1,nodes.num+1):
            i = k-1
            m = 2*i
            if nodes.type[i] == 'slack':
                h[m] = np.inner(H[m],State)
                h[m+1] = np.inner(H[m+1],State)
            elif nodes.type[i] == 'PQ':
                z[m] = (nodes.pwr_flow_values[1][i]*V[i].real + nodes.pwr_flow_values[2][i]*V[i].imag)/(np.abs(V[i])**2)
                z[m+1] = (nodes.pwr_flow_values[1][i]*V[i].imag - nodes.pwr_flow_values[2][i]*V[i].real)/(np.abs(V[i])**2)
                h[m] = np.inner(H[m],State)
                h[m+1] = np.inner(H[m+1],State)
            elif nodes.type[i] == 'PV':
                z[m] = (nodes.pwr_flow_values[1][i]*V[i].real + nodes.pwr_flow_values[2][i]*V[i].imag)/(np.abs(V[i])**2)
                h[m] = np.inner(H[m],State)
                h[m+1] = np.abs(V[i])
                H[m+1][0] = np.cos(np.angle(V[i]))
                H[m+1][1] = np.sin(np.angle(V[i]))
                idx = np.where(znod.Z[i] != 0)
                idx1 = idx[0]+2
                idx2 = idx1[0]+branches.num
                H[m+1][idx1] = - znod.R[i][idx]*np.cos(np.angle(V[i])) - znod.X[i][idx]*np.sin(np.angle(V[i]))
                H[m+1][idx2] = - znod.R[i][idx]*np.sin(np.angle(V[i])) - znod.X[i][idx]*np.cos(np.angle(V[i]))
        
        
        r = np.subtract(z,h)
        Hinv = np.linalg.inv(H)
        Delta_State = np.inner(Hinv,r)
        State = State + Delta_State
        diff = np.amax(np.absolute(Delta_State))
        
        V[0] = State[0] + 1j * State[1]
        Ir = State[2:branches.num+2]
        Ix = State[branches.num+2:]
        Irx = Ir + 1j*Ix
        
        Vcomplex = np.inner((V[0].real + 1j*V[0].imag), np.ones(nodes.num))
        DeltaV = np.inner(znod.Z,Irx)
        V = Vcomplex[0] - DeltaV
        
        num_iter = num_iter+1
        
    I = Ir + 1j*Ix
    Iinj_r = np.zeros(nodes.num)
    Iinj_x = np.zeros(nodes.num)
    for k in range(1,nodes.num+1):
        to = np.where(branches.end==k)
        fr = np.where(branches.start==k)
        Iinj_r[k-1] = np.sum(I[to[0]].real) - np.sum(I[fr[0]].real)
        Iinj_x[k-1] = np.sum(I[to[0]].imag) - np.sum(I[fr[0]].imag)
        
    Iinj = Iinj_r + 1j * Iinj_x
    Sinj_rx = np.multiply(V, np.conj(Iinj))
    
    Sinj = np.real(Sinj_rx) + 1j * np.imag(Sinj_rx)
    S1 = np.multiply(V[branchs.start-1], np.conj(I))
    S2 = np.multiply(V[branchs.end-1], np.conj(I))
    
    return V, I, Iinj, S1, S2, Sinj, num_iter