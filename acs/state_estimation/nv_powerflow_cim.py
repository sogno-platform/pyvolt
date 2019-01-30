import sys
import math
import numpy as np
from network import Ymatrix_calc
from network import BusType

sys.path.append("../../../dataprocessing")
from dpsim_reader import read_data

class PF_Results():
	def __init__(self):
		self.V = [] #Node Voltage
		self.I = [] #Branch Current
		self.Iinj = [] #Injection Current
		self.S1 = [] #Branch Complex Power measured at 1st node of the line
		self.S2 = [] #Branch Complex Power measured at 2nd node of the line
		self.SInj = [] #Complex Power Injection
		self.num_iter = 0 #Number of Iterations of Newton Rapson
		
	def read_data(self, file_name, system):
		"""
		To read the voltages from a input file
		"""
		loadflow_results = read_data(loadflow_results_file)
		#order readed data according to system.nodes
		self.V = np.zeros(len(loadflow_results), dtype=np.complex_)
		for elem in range(len(system.nodes)): 
			self.V[elem] = loadflow_results[system.nodes[elem].uuid]

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
		self.power = complex(0, 0)

class PowerflowResults():	
	def __init__(self):
		self.nodes=[]
		self.branches=[]

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
		
	return V, num_iter
	
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