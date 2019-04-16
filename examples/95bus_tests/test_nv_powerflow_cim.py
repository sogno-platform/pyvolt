import math
import sys

sys.path.append("../../acs/state_estimation")
import network

import py_95bus_network_data
import nv_powerflow_cim
from nv_powerflow import *

class PerUnit:
	def __init__(self, S, V):
		self.S = S
		self.V = V
		self.I = S/V
		self.Z = S/(V**2)
		
		
""" Insert here per unit values of the grid for power and voltage """
S = 100*(10**6)
V = (11*(10**3))/math.sqrt(3)
slackV = 1.02

Base = PerUnit(S,V)
branch, node = py_95bus_network_data.Network_95_nodes(Base, slackV)
V, I, Iinj, S1, S2, Sinj, num_iter = solve(branch, node)

system = network.load_python_data(node, branch, node.type)
res, num_iter = nv_powerflow_cim.solve(system)

# Show numerical comparison
print("res.V==V? " +  str(not np.any(res.get_voltages()-V)))
print("res.I==I? " +  str(not np.any(res.getI()-I)))
print("res.Iinj==Iinj? " +  str(not np.any(res.get_Iinj()-Iinj)))
print("res.S1==S1? " +  str(not np.any(res.get_S1()-S1)))
print("res.S2==S2? " +  str(not np.any(res.get_S2()-S2)))
print("res.Sinj==Sinj? " +  str(not np.any(res.get_Sinj()-Sinj)))