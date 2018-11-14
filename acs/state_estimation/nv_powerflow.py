import numpy as np
import math
        
def Ymatrix_calc(branch, node):
    Ymatrix = np.zeros((node.num,node.num),dtype=np.complex)
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
        
def solve(branch, node):
    """It performs Power Flow by using rectangular node voltage state variables."""
    
    Ymatrix, Adj = Ymatrix_calc(branch,node)
    
    z = np.zeros(2*(node.num))
    h = np.zeros(2*(node.num))
    H = np.zeros((2*(node.num),2*(node.num)))
    
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
            H[m][i] = - np.real(Ymatrix[i][i])
            H[m][i2] = np.imag(Ymatrix[i][i])
            H[m+1][i] = - np.imag(Ymatrix[i][i])
            H[m+1][i2] = - np.real(Ymatrix[i][i])
            idx1 = np.subtract(Adj[i],1)
            idx2 = idx1 + node.num
            H[m][idx1] = - np.real(Ymatrix[i][idx1])
            H[m][idx2] = np.imag(Ymatrix[i][idx1])
            H[m+1][idx1] = - np.imag(Ymatrix[i][idx1])
            H[m+1][idx2] = - np.real(Ymatrix[i][idx1])
        elif t == 'PV':
            z[m+1] = node.pwr_flow_values[2][i]
            H[m][i] = - np.real(Ymatrix[i][i])
            H[m][i2] = np.imag(Ymatrix[i][i])
            idx1 = np.subtract(Adj[i],1)
            idx2 = idx1 + node.num
            H[m][idx1] = - np.real(Ymatrix[i][idx1])
            H[m][idx2] = np.imag(Ymatrix[i][idx1])
    
    epsilon = 10**(-10)
    diff = 5
    V = np.ones(node.num) + 1j* np.zeros(node.num)
    num_iter = 0
    
    State = np.ones(2*branch.num)
    State = np.concatenate((np.array([1,0]),State),axis=0)
    
    while diff > epsilon:
        for k in range(1,node.num+1):
            i = k-1
            m = 2*i
            i2 = i + node.num
            if node.type[i] == 'slack':
                h[m] = np.inner(H[m],State)
                h[m+1] = np.inner(H[m+1],State)
            elif node.type[i] == 'PQ':
                z[m] = (node.pwr_flow_values[1][i]*V[i].real + node.pwr_flow_values[2][i]*V[i].imag)/(np.abs(V[i])**2)
                z[m+1] = (node.pwr_flow_values[1][i]*V[i].imag - node.pwr_flow_values[2][i]*V[i].real)/(np.abs(V[i])**2)
                h[m] = np.inner(H[m],State)
                h[m+1] = np.inner(H[m+1],State)
            elif node.type[i] == 'PV':
                z[m] = (node.pwr_flow_values[1][i]*V[i].real + node.pwr_flow_values[2][i]*V[i].imag)/(np.abs(V[i])**2)
                h[m] = np.inner(H[m],State)
                h[m+1] = np.abs(V[i])
                H[m+1][i] = np.cos(np.angle(V[i]))
                H[m+1][i2] = np.sin(np.angle(V[i]))
        
        r = np.subtract(z,h)
        Hinv = np.linalg.inv(H)
        Delta_State = np.inner(Hinv,r)
        State = State + Delta_State
        diff = np.amax(np.absolute(Delta_State))
        
        V = State[:node.num] + 1j * State[node.num:]
                
        num_iter = num_iter+1
        
    Irx = np.zeros((branch.num),dtype=np.complex)
    for idx in range(branch.num):
        fr = branch.start[idx]-1
        to = branch.end[idx]-1
        Irx[idx] = - (V[fr] - V[to])*Ymatrix[fr][to]
    Ir = np.real(Irx)
    Ix = np.imag(Irx)
    
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
    S1 = np.multiply(V[branch.start-1], np.conj(I))
    S2 = - np.multiply(V[branch.end-1], np.conj(I))
    
    return V, I, Iinj, S1, S2, Sinj, num_iter
    