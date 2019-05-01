import math
import numpy as np
from acs.state_estimation import network
from acs.state_estimation import nv_powerflow
from acs.state_estimation import nv_powerflow_cim
import py_95bus_network_data

class PerUnit:
	def __init__(self, S, V):
		self.S = S
		self.V = V
		self.I = S/V
		self.Z = S/(V**2)
		
def load_python_data(nodes, branches, type):
	"""
	to create a new network.System object from the objects 
	py_95bus_network_data.Node and py_95bus_network_data.Branch
	"""
	system = network.System()
	
	for node_idx in range(0, nodes.num):
		node = network.Node()
		node.index = node_idx
		node.type = network.BusType[type[node_idx]]
		if node.type == network.BusType.slack:
			node.voltage_pu = nodes.P2[node_idx]*np.cos(nodes.Q2[node_idx]) + 1j * nodes.P2[node_idx]*np.sin(nodes.Q2[node_idx])
		elif node.type == network.BusType.PQ:
			node.power_pu = complex(nodes.P2[node_idx], nodes.Q2[node_idx])
		elif node.type == network.BusType.PV:
			pass
		system.nodes.append(node)

	for branch_idx in range(0, branches.num):
		branch = network.Branch()
		branch.r_pu = branches.R[branch_idx]
		branch.x_pu = branches.X[branch_idx]
		branch.z_pu = complex(branch.r_pu, branch.x_pu)
		branch.y_pu = 1/branch.z_pu if (branch.z_pu != 0) else float("inf")
		branch.start_node = system.nodes[branches.start[branch_idx]-1]
		branch.end_node = system.nodes[branches.end[branch_idx]-1]
		system.branches.append(branch)
	
	system.Ymatrix_calc()
	return system

""" Insert here per unit values of the grid for power and voltage """
S = 100*(10**6)
V = (11*(10**3))/math.sqrt(3)
slackV = 1.02

#execute powerflow analysis using nv_powerflow
Base = PerUnit(S,V)
branch, node = py_95bus_network_data.Network_95_nodes(Base, slackV)
V, I, Iinj, S1, S2, Sinj, num_iter = nv_powerflow.solve(branch, node)

#execute powerflow analysis using nv_powerflow_cim
system = load_python_data(node, branch, node.type)
res, num_iter = nv_powerflow_cim.solve(system)

# Show numerical comparison with 6 decimals
print("res.V==V? " +  str(not np.any(np.around(res.get_voltages(),6)-np.around(V,6))))
print("res.Iinj==Iinj? " +  str(not np.any(np.around(res.get_Iinj(),6)-np.around(Iinj,6))))
print("res.I==I? " +  str(not np.any(np.around(res.getI(),6)-np.around(I,6))))
print("res.S1==S1? " +  str(not np.any(np.around(res.get_S1(),6)-np.around(S1,6))))
print("res.S2==S2? " +  str(not np.any(np.around(res.get_S2(),6)-np.around(S2,6))))
print("res.Sinj==Sinj? " +  str(not np.any(np.around(res.get_Sinj(),6)-np.around(Sinj,6))))