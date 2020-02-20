import cmath
import numpy as np
from villas.dataprocessing.readtools import read_timeseries_dpsim


class ResultsNode():
    def __init__(self, topo_node):
        self.topology_node = topo_node
        self.voltage = complex(0, 0)
        self.current = complex(0, 0)
        self.power = complex(0, 0)
        self.voltage_pu = complex(0, 0)
        self.current_pu = complex(0, 0)
        self.power_pu = complex(0, 0)

    def __str__(self):
        str = 'class=Node\n'
        attributes = self.__dict__
        for key in attributes.keys():
            str = str + key + '={}\n'.format(attributes[key])
        return str


class ResultsBranch():
    def __init__(self, topo_branch):
        self.topology_branch = topo_branch
        self.current = complex(0, 0)
        self.power = complex(0, 0)  # complex power flow at branch, measured at initial node
        self.power2 = complex(0, 0)  # complex power flow at branch, measured at final node
        self.current_pu = complex(0, 0)
        self.power_pu = complex(0, 0)  # complex power flow at branch, measured at initial node
        self.power2_pu = complex(0, 0)  # complex power flow at branch, measured at final node

    def __str__(self):
        str = 'class=Branch\n'
        attributes = self.__dict__
        for key in attributes.keys():
            str = str + key + '={}\n'.format(attributes[key])
        return str


class Results():
    def __init__(self, system):
        self.nodes = []
        self.branches = []
        self.Ymatrix = system.Ymatrix
        self.Adjacencies = system.Adjacencies
        for node in system.nodes:
            if node.ideal_connected_with == '':
                self.nodes.append(ResultsNode(topo_node=node))
        for branch in system.branches:
            self.branches.append(ResultsBranch(topo_branch=branch))

    def get_node_by_index(self, index):
        """
        return the node with node.index==index
        """
        for node in self.nodes:
            if (node.topology_node.index==index):
                return node
        
        return None
    
    def read_data_dpsim(self, file_name, pu=False):
        """
        read the voltages from a dpsim input file

        @param file_name
        @param pu: - True if voltages are expressed in per unit system
        """
        loadflow_results = read_timeseries_dpsim(file_name, print_status=False)

        if pu == True:
            for node in self.nodes:
                node.voltage_pu = loadflow_results[node.topology_node.uuid].values[0]
                node.voltage = node.voltage_pu * node.topology_node.baseVoltage
        elif pu == False:
            for node in self.nodes:
                node.voltage = loadflow_results[node.topology_node.uuid].values[0]
                node.voltage_pu = node.voltage / node.topology_node.baseVoltage

    def load_voltages(self, V):
        """
        load the voltages of V-array (result of powerflow_cim.solve)
        """
        for index in range(len(V)):
            node = self.get_node_by_index(index)
            node.voltage_pu = V[index]
            node.voltage = node.voltage_pu * node.topology_node.baseVoltage
    
    def calculate_all(self):
        """
        calculate all quantities of the grid
        """
        self.calculateI()
        self.calculateIinj()
        self.calculateSinj()
        self.calculateIinj()
        self.calculateS1()
        self.calculateS2()

    def calculateI(self):
        """
        To calculate the branch currents
        """
        for branch in self.branches:
            fr = branch.topology_branch.start_node.index
            to = branch.topology_branch.end_node.index
            branch.current_pu = - (self.nodes[fr].voltage_pu - self.nodes[to].voltage_pu) * self.Ymatrix[fr][to]
            branch.current = branch.current_pu * branch.topology_branch.base_current

    def calculateIinj(self):
        """
        calculate current injections at a node
        """
        for node in self.nodes:
            to = complex(0, 0)  # sum of the currents flowing to the node
            fr = complex(0, 0)  # sum of the currents flowing from the node
            for branch in self.branches:
                if node.topology_node.index == branch.topology_branch.start_node.index:
                    fr = fr + branch.current_pu
                if node.topology_node.index == branch.topology_branch.end_node.index:
                    to = to + branch.current_pu
            node.current_pu = to - fr
            node.current = node.current_pu * node.topology_node.base_current

    def calculateSinj(self):
        """
        calculate power injection at a node
        """
        for node in self.nodes:
            node.power_pu = node.voltage_pu * np.conj(node.current_pu)
            node.power = node.power_pu * node.topology_node.base_apparent_power

    def calculateS1(self):
        """
        calculate complex power flow at branch, measured at initial node
        """
        for branch in self.branches:
            branch_index = branch.topology_branch.start_node.index
            for node in self.nodes:
                if branch_index == node.topology_node.index:
                    branch.power_pu = node.voltage_pu * (np.conj(branch.current_pu))
                    branch.power = branch.power_pu * branch.topology_branch.base_apparent_power

    def calculateS2(self):
        """
        calculate complex ower flow at branch, measured at final node
        """
        for branch in self.branches:
            branch_index = branch.topology_branch.end_node.index
            for node in self.nodes:
                if branch_index == node.topology_node.index:
                    branch.power2_pu = -node.voltage_pu * (np.conj(branch.current_pu))
                    branch.power = branch.power2_pu * branch.topology_branch.base_apparent_power

    def get_node(self, index=None, uuid=None):
        """
        returns a PowerflowNode with a certain uuid or a certain index (not both!):
        - if index in not None --> return the PowerflowNode with PowerflowNode.topology_node.index == index
        - if uuid in not None --> return the PowerflowNode with PowerflowNode.topology_node.uuid == uuid
        """
        if index is not None:
            for node in self.nodes:
                if index == node.topology_node.index:
                    return node
        elif uuid is not None:
            for node in self.nodes:
                if uuid == node.topology_node.uuid:
                    return node

    def get_branch(self, uuid):
        """
        returns a PowerflowBranch with a certain uuid
        """
        for branch in self.branches:
            if uuid == branch.topology_branch.uuid:
                return branch

    def get_voltages(self, pu=True):
        """
        get node voltages
        - if pu==True --> voltages are expressed as per-unit
        """
        voltages = np.zeros(len(self.nodes), dtype=np.complex_)
        if pu == True:
            for node in self.nodes:
                voltages[node.topology_node.index] = node.voltage_pu
        elif pu == False:
            for node in self.nodes:
                voltages[node.topology_node.index] = node.voltage

        return voltages

    def get_Iinj(self, pu=True):
        """
        get node currents
        - if pu==True --> voltages are expressed as per-unit
        """
        Iinj = np.zeros(len(self.nodes), dtype=np.complex_)
        if pu == True:
            for node in self.nodes:
                Iinj[node.topology_node.index] = node.current_pu
        elif pu == False:
            for node in self.nodes:
                Iinj[node.topology_node.index] = node.current

        return Iinj

    def get_Sinj(self, pu=True):
        """
        get node power
        - if pu==True --> voltages are expressed as per-unit
        """
        Sinj = np.zeros(len(self.nodes), dtype=np.complex_)
        if pu == True:
            for node in self.nodes:
                Sinj[node.topology_node.index] = node.power_pu
        elif pu == False:
            for node in self.nodes:
                Sinj[node.topology_node.index] = node.power

        return Sinj

    def getI(self, pu=True):
        """
        get branch currents
        - if pu==True --> voltages are expressed as per-unit
        """
        I = np.zeros(len(self.branches), dtype=np.complex_)
        if pu == True:
            for branch_idx in range(len(self.branches)):
                I[branch_idx] = self.branches[branch_idx].current_pu
        elif pu == False:
            for branch_idx in range(len(self.branches)):
                I[branch_idx] = self.branches[branch_idx].current

        return I

    def get_S1(self, pu=True):
        """
        get complex Power flow at branch, measured at initial node
        - if pu==True --> voltages are expressed as per-unit
        """
        S1 = np.zeros(len(self.branches), dtype=np.complex_)
        if pu == True:
            for branch_idx in range(len(self.branches)):
                S1[branch_idx] = self.branches[branch_idx].power_pu
        elif pu == False:
            for branch_idx in range(len(self.branches)):
                S1[branch_idx] = self.branches[branch_idx].power

        return S1

    def get_S2(self, pu=True):
        """
        get complex Power flow at branch, measured at final node
        - if pu==True --> voltages are expressed as per-unit
        """
        S2 = np.zeros(len(self.branches), dtype=np.complex_)
        if pu == True:
            for branch_idx in range(len(self.branches)):
                S2[branch_idx] = self.branches[branch_idx].power2_pu
        elif pu == False:
            for branch_idx in range(len(self.branches)):
                S2[branch_idx] = self.branches[branch_idx].power2

        return S2

    def print_voltages_polar(self):
        """
        for test purposes
        """
        for node in self.nodes:
            print(node.topology_node.uuid + " = " + str(cmath.polar(node.voltage)))
