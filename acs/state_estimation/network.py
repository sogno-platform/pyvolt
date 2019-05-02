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
	def __init__(self, uuid='', base_voltage = 1.0, base_apparent_power=1.0, v_mag=0.0, 
				 v_phase=0.0, p=0.0, q=0.0, index=0,  bus_type='PQ', pu=False):
		self.uuid = uuid
		self.index = index
		self.baseVoltage = base_voltage
		self.base_apparent_power = base_apparent_power
		self.base_current = self.base_apparent_power/self.baseVoltage
		self.type = BusType[bus_type]
		self.voltage = v_mag*np.cos(v_phase) + 1j * v_mag*np.sin(v_phase)
		self.power = complex(p, q)
		self.power_pu = complex(p, q)/self.base_apparent_power
		self.voltage_pu = self.voltage/self.baseVoltage
		
class Branch():
	def __init__(self, uuid='', r=0.0, x=0.0, start_node=None, end_node=None, 
					base_voltage=1.0, base_apparent_power=1.0):
		self.uuid = uuid
		self.baseVoltage = base_voltage
		self.base_apparent_power = base_apparent_power
		self.base_current = self.base_apparent_power/self.baseVoltage
		self.base_impedance = base_voltage**2/self.base_apparent_power
		self.start_node = start_node
		self.end_node = end_node
		self.r = r
		self.x = x
		self.z = self.r + 1j*self.x
		self.y = 1/self.z if (self.z != 0) else float("inf")
		self.r_pu = r/self.base_impedance
		self.x_pu = x/self.base_impedance
		self.z_pu = self.r_pu + 1j*self.x_pu
		self.y_pu = 1/self.z_pu if (self.z_pu != 0) else float("inf")

class System():	
	def __init__(self):
		self.nodes=[]
		self.branches=[]
		self.Ymatrix = np.zeros([], dtype=np.complex)
		self.Adjacencies = np.array([])

	def load_cim_data(self, res, base_apparent_power):
		"""
		To fill the vectors node, branch, bR, bX, P and Q
		"""
		index = 0
		for uuid, element in res.items():
			if element.__class__.__name__=="TopologicalNode":
				vmag = 0.0
				vphase = 0.0
				pInj = 0.0
				qinj = 0.0
				for uuid2, element2 in res.items():
					if element2.__class__.__name__=="SvVoltage":
						if element2.getNodeUUID() == uuid:
							vmag = element2.v
							vphase = element2.angle
							break
				for uuid2, element2 in res.items():
					if element2.__class__.__name__=="SvPowerFlow":
						if element2.getNodeUUID() == uuid:
							pInj += element2.p
							qinj += element2.q
				node_type = self._getNodeType(element)
				base_voltage = element.BaseVoltage.nominalVoltage
				self.nodes.append(Node(uuid=uuid, base_voltage=base_voltage, v_mag=vmag, 
										base_apparent_power=base_apparent_power, v_phase=vphase, 
										p=pInj, q=qinj, index=index, bus_type=node_type))
				index = index+1
		
		for uuid, element in res.items():
			if element.__class__.__name__=="ACLineSegment":
				for node in self.nodes:
					if element.startNodeID==node.uuid:
						startNode = node
						break
				for node in self.nodes:
					if element.endNodeID==node.uuid:
						endNode=node
						break
				base_voltage = element.BaseVoltage.nominalVoltage
				self.branches.append(Branch(uuid=uuid, r=element.r, x=element.x, start_node=startNode, 
											end_node=endNode, base_voltage=base_voltage, 
											base_apparent_power=base_apparent_power))
				continue
			elif element.__class__.__name__=="PowerTransformer":
				for i in range(len(self.nodes)):
					if element.startNodeID==self.nodes[i].uuid:
						startNode=self.nodes[i]
						break
				for i in range(len(self.nodes)):
					if element.endNodeID==self.nodes[i].uuid:
						endNode=self.nodes[i]
						break
				#base voltage = high voltage side (=primaryConnection)
				base_voltage = element.primaryConnection.BaseVoltage.nominalVoltage
				self.branches.append(Branch(uuid=uuid, r=element.primaryConnection.r, x=element.primaryConnection.x, 
											start_node=startNode, end_node=endNode, base_voltage=base_voltage, 
											base_apparent_power=base_apparent_power))
				continue
			else:
				continue

		#calculate admitance matrix
		self.Ymatrix_calc()

	def _getNodeType(self, node):
		"""
		return the type of a node: PQ, PV or SLACK

		@param node: element of class cimpy.Topology.TopologicalNode
		"""
		for terminal in node.Terminal:
			if terminal.ConductingEquipment.__class__.__name__=="ExternalNetworkInjection":
				return "SLACK"
			elif terminal.ConductingEquipment.__class__.__name__=="SynchronousMachine":
				return "PV"
		return "PQ"

	def Ymatrix_calc(self):
		nodes_num = len(self.nodes)
		self.Ymatrix = np.zeros((nodes_num, nodes_num), dtype=np.complex)
		self.Adjacencies = [[] for _ in range(nodes_num)]
		for branch in self.branches:
			fr = branch.start_node.index
			to = branch.end_node.index
			self.Ymatrix[fr][to] -= branch.y_pu
			self.Ymatrix[to][fr] -= branch.y_pu
			self.Ymatrix[fr][fr] += branch.y_pu
			self.Ymatrix[to][to] += branch.y_pu
			self.Adjacencies[fr].append(to+1)	#to + 1???
			self.Adjacencies[to].append(fr+1)	#fr + 1???