import numpy as np

class Node():
	def __init__(self, name, uuid, v_mag, v_phase, p, q, index):
		self.index = index
		self.name = name
		self.uuid = uuid
		self.power = complex(p, q)
		self.voltage = v_mag*np.cos(v_phase) + 1j * v_mag*np.sin(v_phase)
				
class Branch():
	def __init__(self, r, x, start_node, end_node):
		self.r = r
		self.x = x
		self.start_node = start_node
		self.end_node = end_node

class System():	
	def __init__(self):
		self.nodes=[]
		self.branches=[]
		self.bR=[]
		self.bX=[]
		self.P=[]
		self.Q=[]
	
	def load_cim_data(res):
		#this function is used to fill the vectors node, branch, bR, bX, P and Q
		for key, value in res.items():
			if value.__class__.__name__=="TopologicalNode":
				self.P.append(value.pInjection)
				self.Q.append(value.qInjection)
				index=len(self.P)-1
				self.nodes.append(Node(value.name, value.mRID, value.pInjection, value.qInjection, index))
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


def load_python_data(nodes, branches):
	system = System()
	
	for node_idx in range(0, nodes.num):
		if node_idx == 0:			
			system.nodes.append(Node('', '', nodes.P2[0], nodes.Q2[0], 0, 0, node_idx))
		else:
			system.nodes.append(Node('', '', 0, 0, nodes.P2[node_idx], nodes.Q2[node_idx], node_idx))
	
	for branch_idx in range(0, branches.num):
		system.branches.append(Branch(branches.R[branch_idx], branches.X[branch_idx], 
		system.nodes[branches.start[branch_idx]-1], system.nodes[branches.end[branch_idx]-1]))

	return system