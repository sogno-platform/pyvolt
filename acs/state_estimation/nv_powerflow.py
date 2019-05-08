import numpy as np
from .network import BusType
from .results import Results

def solve(system):
	"""It performs Power Flow by using rectangular node voltage state variables."""
	
	nodes_num = len(system.nodes)
	branches_num = len(system.branches)
	
	z = np.zeros(2*nodes_num)
	h = np.zeros(2*nodes_num)
	H = np.zeros((2*nodes_num,2*nodes_num))
	
	for i in range(0,nodes_num):
		m = 2*i
		i2 = i + nodes_num
		node_type = system.nodes[i].type 
		if  node_type == BusType.SLACK:
			z[m] = np.real(system.nodes[i].voltage_pu)
			z[m+1] = np.imag(system.nodes[i].voltage_pu)
			H[m][i] = 1
			H[m+1][i2] = 1
		elif node_type is BusType.PQ:
			H[m][i] = - np.real(system.Ymatrix[i][i])
			H[m][i2] = np.imag(system.Ymatrix[i][i])
			H[m+1][i] = - np.imag(system.Ymatrix[i][i])
			H[m+1][i2] = - np.real(system.Ymatrix[i][i])
			idx1 = np.subtract(system.Adjacencies[i],1)
			idx2 = idx1 + nodes_num
			H[m][idx1] = - np.real(system.Ymatrix[i][idx1])
			H[m][idx2] = np.imag(system.Ymatrix[i][idx1])
			H[m+1][idx1] = - np.imag(system.Ymatrix[i][idx1])
			H[m+1][idx2] = - np.real(system.Ymatrix[i][idx1])
		elif node_type is BusType.PV:
			z[m+1] = np.real(system.nodes[i].power)
			H[m][i] = - np.real(system.Ymatrix[i][i])
			H[m][i2] = np.imag(system.Ymatrix[i][i])
			idx1 = np.subtract(system.Adjacencies[i],1)
			idx2 = idx1 + nodes_num
			H[m][idx1] = - np.real(system.Ymatrix[i][idx1])
			H[m][idx2] = np.imag(system.Ymatrix[i][idx1])
	
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
			node_type = system.nodes[i].type 
			if node_type is BusType.SLACK:
				h[m] = np.inner(H[m],State)
				h[m+1] = np.inner(H[m+1],State)
			elif node_type is BusType.PQ:
				z[m] = (np.real(system.nodes[i].power_pu)*np.real(V[i]) + np.imag(system.nodes[i].power_pu)*np.imag(V[i]))/(np.abs(V[i])**2)
				z[m+1] = (np.real(system.nodes[i].power_pu)*np.imag(V[i]) - np.imag(system.nodes[i].power_pu)*np.real(V[i]))/(np.abs(V[i])**2)
				h[m] = np.inner(H[m],State)
				h[m+1] = np.inner(H[m+1],State)
			elif node_type is BusType.PV:
				z[m] = (np.real(system.nodes[i].power_pu)*np.real(V[i])+ np.imag(system.nodes[i].power_pu)*np.imag(V[i]))(np.abs(V[i])**2)
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
		
	# calculate all the other quantities of the grid
	powerflow_results = Results(system)
	powerflow_results.load_voltages(V)
	powerflow_results.calculate_all()

	return powerflow_results, num_iter