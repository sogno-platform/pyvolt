import numpy as np
from enum import Enum

class BusType(Enum):
	SLACK = 1
	slack = 1
	PV = 2
	pv = 2
	PQ = 3
	pq = 3


class Node():
	def __init__(self, name='', uuid='', v_mag=0.0, v_phase=0.0, p=0.0, q=0.0, index=0,  bus_type='PQ'):
		self.index = index
		self.name = name
		self.uuid = uuid
		self.power = complex(p, q)
		self.voltage = v_mag*np.cos(v_phase) + 1j * v_mag*np.sin(v_phase)
		self.type = BusType[bus_type]

class Branch():
	def __init__(self, r, x, start_node, end_node):
		self.r = r
		self.x = x
		self.start_node = start_node
		self.end_node = end_node
		self.z = self.r + 1j*self.x
		self.y = 1/self.z if (self.z != 0) else float("inf")

class System():	
	def __init__(self):
		self.nodes=[]
		self.branches=[]
		self.bR=[]
		self.bX=[]
		self.P=[]
		self.Q=[]
		self.Ymatrix = np.zeros((1, 1),dtype=np.complex)
		self.Adjacencies = np.array([])

	#def load_cim_data(self, res, Sb, Vb, Zb):
	def load_cim_data(self, res):
		#this function is used to fill the vectors node, branch, bR, bX, P and Q
		for key, value in res.items():
			if value.__class__.__name__=="TopologicalNode":
				self.P.append(value.pInjection)
				self.Q.append(value.qInjection)
				index=len(self.P)-1
				self.nodes.append(Node(name=value.name, uuid=value.mRID, p=value.pInjection, q=value.qInjection, index=index))
		for key, value in res.items():
			if value.__class__.__name__=="ACLineSegment":
				length=value.length
				if length==0.0:
					length=1.0
				self.bR.append(value.r*length)
				self.bX.append(value.x*length)
				for i in range(len(self.nodes)):
					if value.startNodeID==self.nodes[i].uuid:
						startNode=self.nodes[i]
						break
				for i in range(len(self.nodes)):
					if value.endNodeID==self.nodes[i].uuid:
						endNode=self.nodes[i]
						break
				self.branches.append(Branch(value.r, value.x, startNode, endNode))	
			elif value.__class__.__name__=="PowerTransformer":
				self.bR.append(value.primaryConnection.r)
				self.bX.append(value.primaryConnection.x)
				for i in range(len(self.nodes)):
					if value.startNodeID==self.nodes[i].uuid:
						startNode=self.nodes[i]
						break
				for i in range(len(self.nodes)):
					if value.endNodeID==self.nodes[i].uuid:
						endNode=self.nodes[i]
						break
				self.branches.append(Branch(value.primaryConnection.r, value.primaryConnection.x, startNode, endNode))
			else:
				continue

		#determine the impedance matrix
		#self.Ymatrix_calc(Zb)
		self.Ymatrix_calc()

	#def Ymatrix_calc(self, Zb):
	def Ymatrix_calc(self):
		"""
		@param Zb: base value of impedance
		"""
		nodes_num = len(self.nodes)
		self.Ymatrix = np.zeros((nodes_num, nodes_num),dtype=np.complex)
		self.Adjacencies = [[] for _ in range(nodes_num)]
		for branch in self.branches:
			fr = branch.start_node.index
			to = branch.end_node.index
			self.Ymatrix[fr][to] -= branch.y
			self.Ymatrix[to][fr] -= branch.y
			self.Ymatrix[fr][fr] += branch.y
			self.Ymatrix[to][to] += branch.y
			self.Adjacencies[fr].append(to+1)	#to + 1???
			self.Adjacencies[to].append(fr+1)	#fr + 1???

		#self.Ymatrix = self.Ymatrix*Zb

def load_python_data(nodes, branches, type):
	system = System()
	
	for node_idx in range(0, nodes.num):
		if BusType[type[node_idx]] == BusType.slack:		
			system.nodes.append(Node(v_mag=nodes.P2[0], v_phase=nodes.Q2[0], p=0, q=0, index=node_idx, bus_type=type[node_idx]))
		else:
			system.nodes.append(Node(v_mag=0, v_phase=0, p=nodes.P2[node_idx], q=nodes.Q2[node_idx], index=node_idx, bus_type=type[node_idx]))
	
	for branch_idx in range(0, branches.num):
		system.branches.append(Branch(branches.R[branch_idx], branches.X[branch_idx], 
		system.nodes[branches.start[branch_idx]-1], system.nodes[branches.end[branch_idx]-1]))
	
	system.Ymatrix_calc()
	return system