import sys
import math
import numpy as np
from network import Ymatrix_calc
from network import BusType

sys.path.append("../../../dataprocessing")
from villas.dataprocessing.readtools import read_timeseries_dpsim

class PowerflowNode():
	def __init__(self, topo_node):		
		self.topology_node = topo_node		
		self.voltage = complex(0, 0)		
		self.current = complex(0, 0)
		self.power = complex(0, 0)

class PowerflowBranch():
	def __init__(self, topo_branch):
		self.topology_branch = topo_branch
		self.current = complex(0, 0)
		self.power = complex(0, 0)			#complex power flow at branch, measured at initial node
		self.power2 = complex(0, 0)			#complex power flow at branch, measured at final node 
		
class PowerflowResults():	
	def __init__(self, system):
		self.nodes=[]
		self.branches=[]
		self.Ymatrix=system.Ymatrix
		self.Adjacencies=system.Adjacencies
		for node in system.nodes:
			self.nodes.append(PowerflowNode(topo_node=node))
		for branch in system.branches:
			self.branches.append(PowerflowBranch(topo_branch=branch))
			
	def read_data_dpsim(self, file_name):
		"""
		read the voltages from a dpsim input file
		"""
		loadflow_results = read_timeseries_dpsim(file_name, print_status=False)
		for node in self.nodes:
			node.V = loadflow_results[node.topology_node.uuid].values[0]

	def load_voltages(self, V):
		"""
		load the voltages of V-array (result of powerflow_cim.solve)
		"""
		for index in range(len(V)):
			for elem in self.nodes:
				if elem.topology_node.index == index:
					elem.voltage = V[index]
					
	def calculateI(self):
		"""
		To calculate the branch currents
		"""	
		for branch in self.branches:
			fr = branch.topology_branch.start_node.index
			to = branch.topology_branch.end_node.index
			branch.current = - (self.nodes[fr].voltage - self.nodes[to].voltage)*self.Ymatrix[fr][to]
	
	def calculateIinj(self):
		"""
		calculate current injections at a node
		"""
		for node in self.nodes:
			to=complex(0, 0)	#sum of the currents flowing to the node
			fr=complex(0, 0)	#sum of the currents flowing from the node
			for branch in self.branches:
				if node.topology_node.index==branch.topology_branch.start_node.index:
					fr = fr + branch.current
				if node.topology_node.index==branch.topology_branch.end_node.index:
					to = to + branch.current
			node.current = to - fr
	
	def calculateSinj(self):
		"""
		calculate power injection at a node
		"""
		for node in self.nodes:
			node.power = node.voltage*np.conj(node.current)
	
	def calculateS1(self):
		"""
		calculate complex power flow at branch, measured at initial node
		"""
		for branch in self.branches:
			branch_index = branch.topology_branch.start_node.index
			for node in self.nodes:
				if branch_index == node.topology_node.index:
					branch.power = node.voltage*(np.conj(branch.current))
	
	def calculateS2(self):
		"""
		calculate complex ower flow at branch, measured at final node 
		"""
		for branch in self.branches:
			branch_index = branch.topology_branch.end_node.index
			for node in self.nodes:
				if branch_index == node.topology_node.index:
					branch.power2 = -node.voltage*(np.conj(branch.current))
					
	def get_voltages(self):
		"""
		get complex Power Injection at nodes
		for a test purpose
		"""
		voltages = np.zeros(len(self.nodes), dtype=np.complex_)
		for node in self.nodes:
			voltages[node.topology_node.index] = node.voltage
		return voltages
	
	def get_Iinj(self):
		"""
		get node currents 
		for a test purpose
		"""
		Iinj = np.zeros(len(self.nodes), dtype=np.complex_)
		for node in self.nodes:
			Iinj[node.topology_node.index] = node.current
		return Iinj
		
	def get_Sinj(self):
		"""
		get node power 
		for a test purpose
		"""
		Sinj = np.zeros(len(self.nodes), dtype=np.complex_)
		for node in self.nodes:
			Sinj[node.topology_node.index] = node.power
		return Sinj
		
	
	def getI(self):
		"""
		get branch currents 
		for a test purpose
		"""
		I = np.array([])
		for branch in self.branches:
			I = np.append(I, branch.current)
		return I
		
	def get_S1(self):
		"""
		get complex Power flow at branch, measured at initial node
		for a test purpose
		"""
		S1 = np.array([])
		for branch in self.branches:
			S1 = np.append(S, branch.power)
		return S1
		
	def get_S2(self):
		"""
		get complex Power flow at branch, measured at final node  
		for a test purpose
		"""
		S2 = np.array([])
		for branch in self.branches:
			S2 = np.append(S, branch.power2)
		return S2
		
def solve(system):
	"""It performs Power Flow by using rectangular node voltage state variables."""
	
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
		elif type is BusType.PV:
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
		
	results = PowerflowResults(system)
	results.load_voltages(V)
	results.calculateI()
	results.calculateIinj()
	results.calculateSinj()
	results.calculateI()
	results.calculateS()
	
	return results, num_iter
	
def calculateI(system, V):
	"""
	To calculate the branch currents
	"""
	Ymatrix, Adj = Ymatrix_calc(system)
	I = np.zeros((len(system.branches)), dtype=np.complex)
	for idx in range(len(system.branches)):
		fr = system.branches[idx].start_node.index
		to = system.branches[idx].end_node.index
		I[idx] = - (V[fr] - V[to])*Ymatrix[fr][to]
	
	return I
		
def calculateInj(system, I):
	"""
	To calculate current injections at a node
	"""
	Iinj = np.zeros((len(system.nodes)), dtype=np.complex)
	for k in range(0, (len(system.nodes))):
		to=[]
		fr=[]
		for m in range(len(system.branches)):
			if k==system.branches[m].start_node.index:
				fr.append(m)
			if k==system.branches[m].end_node.index:
				to.append(m)
		Iinj[k] = np.sum(I[to]) - np.sum(I[fr])
	
	return Iinj

def calculateS1(system, V, I):
	"""
	To calculate powerflow on branches
	"""
	S1 = np.zeros((len(system.branches)), dtype=np.complex)
	for i in range(0, len(system.branches)):
		S1[i] = V[system.branches[i].start_node.index]*(np.conj(I[i]))
	return S1
	
def calculateS2(system, V, I):
	"""
	To calculate powerflow on branches
	"""
	S2 = np.zeros((len(system.branches)), dtype=np.complex)
	for i in range(0, len(system.branches)):
		S2[i] = -V[system.branches[i].end_node.index]*(np.conj(I[i]))
	return S2
	
def calculateSinj(V, Iinj):
	"""
	To calculate power injection at a node
	"""
	Sinj = np.multiply(V, np.conj(Iinj))
	return Sinj