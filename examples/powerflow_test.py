import numpy
import math
import matplotlib as plt

from acs.state_estimation import *  
import py_95bus_network_data

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

Vtrue, Itrue, Iinjtrue, S1true, S2true, Sinjtrue, num_iter = bc_powerflow.solve(branch, node)

system = network.load_python_data(node, branch)
#Vtrue_cim, Itrue_cim, Iinjtrue_cim, S1true_cim, S2true_cim, Sinjtrue_cim, num_iter_cim = nv_powerflow_cim.solve(system)

print('done')
