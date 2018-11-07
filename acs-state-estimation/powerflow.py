import numpy as np
import math

class Znod:
    def __init__(self, branches, nodes):
        # Complex impedance matrix
        self.Z = np.zeros( (nodes.num-1,branches.num), dtype=np.complex_ )
        # Add impedance of each branch
        for k in range(1, nodes.num):
            i = k-1
            aa = branches.start[i]-2
            if aa == -1:
                self.Z[i] = self.Z[i]
            else:
                self.Z[i] = self.Z[aa]
            self.Z[i][i] = branches.Z[i]
        
        # Add one row of zeros
        zzz = np.zeros((1, branches.num))
        self.Z = np.concatenate( (zzz,self.Z), axis=0 )
        # Resistance matrix
        self.R = self.Z.real
        # Reactance matrix
        self.X = self.Z.imag
            
def BC_power_flow(branch, node):
    """It performs Power Flow by using branch current state variables."""
      
    znod = Znod(branch,node)
    # real matrices storing complex numbers 
    # real in pos n and imag in pos n+1
    z = np.zeros(2*(branch.num+1))
    h = np.zeros(2*(branch.num+1))
    H = np.zeros((2*(branch.num+1),2*(branch.num+1)))
    
    for k in range(1, node.num+1):
        # i is matrix index
        # the size of z,h,H is twice the number of branches
        # but here we are using the node number to work on the 
        # matrix. Is this not dangerous?
        i = k-1
        m = 2*i
        if node.type[i] == 'slack':
            V = node.pwr_flow_values[1][i]
            theta = node.pwr_flow_values[2][i]
            # real node voltage
            z[m] = V*math.cos(theta)
            # imag node voltage
            z[m+1] = V*math.sin(theta)
            H[m][0] = 1;
            H[m][2:] = np.concatenate((-znod.R[i],znod.X[i]), axis=0)
            H[m+1][1] = 1
            H[m+1][2:] = np.concatenate((-znod.X[i],-znod.R[i]), axis=0)
        elif node.type[i] == 'PQ':
            to = np.where(branch.end==k)
            fr = np.where(branch.start==k)
            to1 = to[0]+2
            to2 = to1+branch.num
            fr1 = fr[0]+2
            fr2 = fr1+branch.num
            H[m][to1] = 1
            H[m+1][to2] = 1
            H[m][fr1] = -1
            H[m+1][fr2] = -1
        elif node.type[i] == 'PV':
            z[m+1] = node.pwr_flow_values[2][i]
            to = np.where(branch.end==k)
            fr = np.where(branch.start==k)
            to1 = to[0]+2
            fr1 = fr[0]+2
            H[m][to1] = 1
            H[m][fr1] = -1
    
    print(z)
    print(h)
    print(H)

    epsilon = 5
    V = np.ones(node.num) + 1j* np.zeros(node.num)
    num_iter = 0
    
    State = np.zeros(2*branch.num)
    State = np.concatenate( (np.array([1,0]),State), axis=0 )
    
    while epsilon > 10**(-10):
        for k in range(1,node.num+1):
            i = k-1
            m = 2*i
            t = node.type[i]
            if t == 'slack':
                h[m] = np.inner(H[m],State)
                h[m+1] = np.inner(H[m+1],State)
            elif t == 'PQ':
                z[m] = (node.pwr_flow_values[1][i]*V[i].real + node.pwr_flow_values[2][i]*V[i].imag)/(np.abs(V[i])**2)
                z[m+1] = (node.pwr_flow_values[1][i]*V[i].imag - node.pwr_flow_values[2][i]*V[i].real)/(np.abs(V[i])**2)
                h[m] = np.inner(H[m],State)
                h[m+1] = np.inner(H[m+1],State)
            elif t == 'PV':
                z[m] = (node.pwr_flow_values[1][i]*V[i].real + node.pwr_flow_values[2][i]*V[i].imag)/(np.abs(V[i])**2)
                h[m] = np.inner(H[m],State)
                h[m+1] = np.abs(V[i])
                H[m+1][0] = np.cos(np.angle(V[i]))
                H[m+1][1] = np.sin(np.angle(V[i]))
                idx = np.where(znod.Z[i] != 0)
                idx1 = idx[0]+2
                idx2 = idx1[0]+branch.num
                H[m+1][idx1] = - znod.R[i][idx]*np.cos(np.angle(V[i])) - znod.X[i][idx]*np.sin(np.angle(V[i]))
                H[m+1][idx2] = - znod.R[i][idx]*np.sin(np.angle(V[i])) - znod.X[i][idx]*np.cos(np.angle(V[i]))
        
        r = np.subtract(z,h)
        Hinv = np.linalg.inv(H)
        Delta_State = np.inner(Hinv,r)
        State = State + Delta_State
        epsilon = np.amax(np.absolute(Delta_State))
        
        V[0] = State[0] + 1j * State[1]
        Ir = State[2:branch.num+2]
        Ix = State[branch.num+2:]
        Irx = Ir + 1j*Ix
        
        Vcomplex = np.inner((V[0].real + 1j*V[0].imag), np.ones(node.num))
        DeltaV = np.inner(znod.Z,Irx)
        V = Vcomplex[0] - DeltaV
        
        num_iter = num_iter+1
        
    I = Ir + 1j*Ix
    Iinj_r = np.zeros(node.num)
    Iinj_x = np.zeros(node.num)

    for k in range(1,node.num+1):
        to = np.where(branch.end==k)
        fr = np.where(branch.start==k)
        Iinj_r[k-1] = np.sum(I[to[0]].real) - np.sum(I[fr[0]].real)
        Iinj_x[k-1] = np.sum(I[to[0]].imag) - np.sum(I[fr[0]].imag)
        
    Iinj = Iinj_r + 1j * Iinj_x
    Sinj_rx = np.multiply(V, np.conj(Iinj))
    
    Sinj = np.real(Sinj_rx) + 1j * np.imag(Sinj_rx)
    S1_rx = np.multiply(V[branch.start-1], np.conj(I))
    S2_rx = np.multiply(V[branch.end-1], np.conj(I))
    S1 = np.real(S1_rx) + 1j * np.imag(S1_rx)
    S2 = np.real(S2_rx) + 1j * np.imag(S2_rx)
    
    return V, I, Iinj, S1, S2, Sinj, num_iter