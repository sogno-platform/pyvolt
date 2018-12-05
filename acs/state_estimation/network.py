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
		self.y = 1/self.z if (self.z != 0) else 0

class System():	
	def __init__(self):
		self.nodes=[]
		self.branches=[]
		self.bR=[]
		self.bX=[]
		self.P=[]
		self.Q=[]

	def load_cim_data(self, res):
		#this function is used to fill the vectors node, branch, bR, bX, P and Q
		for key, value in res.items():
			if value.__class__.__name__=="TopologicalNode":
				self.P.append(value.pInjection)
				self.Q.append(value.qInjection)
				index=len(self.P)-1
				self.nodes.append(Node(name=value.name, uuid=value.mRID, p=value.pInjection, q=value.qInjection, index=index))
			elif value.__class__.__name__=="ACLineSegment":
				length=value.length
				if length==0.0:
					length=1.0
				self.bR.append(value.r*length)
				self.bX.append(value.x*length)
				startNode=res[value.startNodeID]
				endNode=res[value.endNodeID]
				self.branches.append(Branch(value.r, value.x, startNode, endNode))	
			elif value.__class__.__name__=="PowerTransformer":
				self.bR.append(value.primaryConnection.r)
				self.bX.append(value.primaryConnection.x)
				startNode=res[value.startNodeID]
				endNode=res[value.endNodeID]
				self.branches.append(Branch(value.primaryConnection.r, value.primaryConnection.x, startNode, endNode))
			else:
				continue


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

	return system