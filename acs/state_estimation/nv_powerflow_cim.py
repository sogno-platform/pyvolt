import numpy as np
import math
from network import BusType
		
def Ymatrix_calc(system):
	nodes_num = len(system.nodes)
	branches_num = len(system.branches)
	Ymatrix = np.zeros((nodes_num, nodes_num),dtype=np.complex)
	Adjacencies = [[] for _ in range(nodes_num)]
	for index in range(branches_num):
		fr = system.branches[index].start_node.index
		to = system.branches[index].end_node.index
		Ymatrix[fr][to] -= system.branches[index].y
		Ymatrix[to][fr] -= system.branches[index].y
		Ymatrix[fr][fr] += system.branches[index].y
		Ymatrix[to][to] += system.branches[index].y
		Adjacencies[fr].append(to+1)
		Adjacencies[to].append(fr+1)
	return Ymatrix, Adjacencies
		
def solve(system):
	"""It performs Power Flow by using rectangular node voltage state variables."""
	
	Ymatrix, Adj = Ymatrix_calc(system)
	
	nodes_num = len(system.nodes)
	branches_num = len(system.branches)
	
	z = np.zeros(2*(nodes_num))
	h = np.zeros(2*(nodes_num))
	H = np.zeros((2*(nodes_num),2*(nodes_num)))
	
	for i in range(0,nodes_num):
		m = 2*i
		i2 = i + nodes_num
		type = system.nodes[i].type 
		if  type is BusType.SLACK:
			z[m] = np.real(system.nodes[i].voltage)
			z[m+1] = np.imag(system.nodes[i].voltage)
			H[m][i] = 1
			H[m+1][i2] = 1
		elif type is BusType.PQ:
			H[m][i] = - np.real(Ymatrix[i][i])
			H[m][i2] = np.imag(Ymatrix[i][i])
			H[m+1][i] = - np.imag(Ymatrix[i][i])
			H[m+1][i2] = - np.real(Ymatrix[i][i])
			idx1 = np.subtract(Adj[i],1)
			idx2 = idx1 + nodes_num
			H[m][idx1] = - np.real(Ymatrix[i][idx1])
			H[m][idx2] = np.imag(Ymatrix[i][idx1])
			H[m+1][idx1] = - np.imag(Ymatrix[i][idx1])
			H[m+1][idx2] = - np.real(Ymatrix[i][idx1])
		elif type is BusType.PV:
			z[m+1] = np.real(system.nodes[i].power)
			H[m][i] = - np.real(Ymatrix[i][i])
			H[m][i2] = np.imag(Ymatrix[i][i])
			idx1 = np.subtract(Adj[i],1)
			idx2 = idx1 + nodes_num
			H[m][idx1] = - np.real(Ymatrix[i][idx1])
			H[m][idx2] = np.imag(Ymatrix[i][idx1])
	
	epsilon = 10**(-10)
	diff = 5
	V = np.ones(nodes_num) + 1j* np.zeros(nodes_num)
	num_iter = 0
	
	State = np.ones(2*branches_num)
	State = np.concatenate((np.array([1,0]),State),axis=0)
	
	while diff > epsilon:
		for i in range(0, nodes_num):
			m = 2*i
			i2 = i + nodes_num
			type = system.nodes[i].type 
			if type is BusType.SLACK:
				h[m] = np.inner(H[m],State)
				h[m+1] = np.inner(H[m+1],State)
			elif type is BusType.PQ:
				z[m] = (np.real(system.nodes[i].power)*V[i].real + np.imag(system.nodes[i].power)*V[i].imag)/(np.abs(V[i])**2)
				z[m+1] = (np.real(system.nodes[i].power)*V[i].imag - np.imag(system.nodes[i].power)*V[i].real)/(np.abs(V[i])**2)
				h[m] = np.inner(H[m],State)
				h[m+1] = np.inner(H[m+1],State)
			elif type is BusType.PV:
				z[m] = (np.real(system.nodes[i].power)*V[i].real + np.imag(system.nodes[i].power)*V[i].imag)/(np.abs(V[i])**2)
				h[m] = np.inner(H[m],State)
				h[m+1] = np.abs(V[i])
				H[m+1][i] = np.cos(np.angle(V[i]))
				H[m+1][i2] = np.sin(np.angle(V[i]))
		
		r = np.subtract(z,h)
		Hinv = np.linalg.inv(H)
		Delta_State = np.inner(Hinv,r)
		State = State + Delta_State
		diff = np.amax(np.absolute(Delta_State))
		
		V = State[:nodes_num] + 1j * State[nodes_num:]
				
		num_iter = num_iter+1
		
	Irx = np.zeros((branches_num),dtype=np.complex)
	for idx in range(branches_num):
		fr = system.branches[idx].start_node.index
		to = system.branches[idx].end_node.index
		Irx[idx] = - (V[fr] - V[to])*Ymatrix[fr][to]
	Ir = np.real(Irx)
	Ix = np.imag(Irx)
	
	I = Ir + 1j*Ix
	Iinj_r = np.zeros(nodes_num)
	Iinj_x = np.zeros(nodes_num)
	
	for k in range(0, nodes_num):
		to=[]
		fr=[]
		for m in range(len(system.branches)):
			if k==system.branches[m].start_node.index:
				fr.append(m)
			if k==system.branches[m].end_node.index:
				to.append(m)

		Iinj_r[k] = np.sum(I[to].real) - np.sum(I[fr].real)
		Iinj_x[k] = np.sum(I[to].imag) - np.sum(I[fr].imag)
		
	Iinj = Iinj_r + 1j * Iinj_x
	Sinj_rx = np.multiply(V, np.conj(Iinj))
	Sinj = np.real(Sinj_rx) + 1j * np.imag(Sinj_rx)
	
	S1=np.array([])
	S2=np.array([])
	I_conj=np.conj(I)
	for i in range(0, branches_num):
		S1=np.append(S1, V[system.branches[i].start_node.index]*(I_conj[i]))
		S2=np.append(S2, -V[system.branches[i].end_node.index]*(I_conj[i]))

	return V, I, Iinj, S1, S2, Sinj, num_iter
	