import logging
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
    def __init__(self, uuid='', name='', base_voltage=1.0, base_apparent_power=1.0, v_mag=0.0,
                 v_phase=0.0, p=0.0, q=0.0, index=0, bus_type='PQ', pu=False):
        self.uuid = uuid
        self.name = name
        self.index = index
        self.baseVoltage = base_voltage
        self.base_apparent_power = base_apparent_power
        self.base_current = self.base_apparent_power / self.baseVoltage / np.sqrt(3)
        self.type = BusType["PQ"]
        self.voltage = v_mag * np.cos(np.radians(v_phase)) + 1j * v_mag * np.sin(np.radians(v_phase))
        self.power = complex(p, q)
        self.power_pu = complex(p, q) / self.base_apparent_power
        self.voltage_pu = self.voltage / self.baseVoltage

    def __str__(self):
        str = 'class=Node\n'
        attributes = self.__dict__
        for key in attributes.keys():
            str = str + key + '={}\n'.format(attributes[key])
        return str


class Branch():
    def __init__(self, uuid='', r=0.0, x=0.0, start_node=None, end_node=None,
                 base_voltage=1.0, base_apparent_power=1.0):
        self.uuid = uuid
        self.baseVoltage = base_voltage
        self.base_apparent_power = base_apparent_power
        self.base_current = self.base_apparent_power / self.baseVoltage / np.sqrt(3)
        self.base_impedance = base_voltage ** 2 / self.base_apparent_power
        self.start_node = start_node
        self.end_node = end_node
        self.r = r
        self.x = x
        self.z = self.r + 1j * self.x
        self.y = 1 / self.z if (self.z != 0) else float("inf")
        self.r_pu = r / self.base_impedance
        self.x_pu = x / self.base_impedance
        self.z_pu = self.r_pu + 1j * self.x_pu
        self.y_pu = 1 / self.z_pu if (self.z_pu != 0) else float("inf")

    def __str__(self):
        str = 'class=Branch\n'
        attributes = self.__dict__
        for key in attributes.keys():
            str = str + key + '={}\n'.format(attributes[key])
        return str


class Breaker():
    def __init__(self, node_load, node_breaker, branch_connected_with_node_breaker, connected=True):
        """
        :param node_load: node which is connected with pq loads
        :param node_breaker:
        :param branch_connected_with_node_breaker: ACLineSegment which is connected with node_breaker
        :param connected: True if the breaker is considered closed and False if the broker is open 
        """
        #node_breaker and node_load are used when the breaker is open
        #only node_breaker is used when the breaker is closed
        self.node_breaker = node_breaker
        self.node_load = node_load
        self.branch_connected_with_node_breaker = branch_connected_with_node_breaker
        self.connected = connected

    def __str__(self):
        str = 'class=Breaker\n'
        attributes = self.__dict__
        for key in attributes.keys():
            str = str + key + '={}\n'.format(attributes[key])
        return str

    def open_breaker(self, system):
        self.connected == False
        node_breaker_uuid = self.node_breaker.uuid
        
        #check if node_breaker_uuid is in system.nodes
        is_node_in_list = system.get_node_by_uuid(node_breaker_uuid)
        
        #add node_breaker to list system.nodes
        if not is_node_in_list:
            system.nodes.append(self.node_breaker)
        
        #connect branch_connected_with_node_breaker with node_load
        if (self.branch_connected_with_node_breaker.start_node.uuid==self.node_load.uuid):
            self.branch_connected_with_node_breaker.start_node = self.node_breaker
        elif (self.branch_connected_with_node_breaker.end_node.uuid==self.node_load.uuid):
            self.branch_connected_with_node_breaker.end_node = self.node_breaker

        #reenumerate nodes
        system.reindex_nodes_list()

    def close_breaker(self, system):
        self.connected == True
        
        #remove node_breaker from list system.nodes
        node_breaker_uuid = self.node_breaker.uuid
        node_to_remove = system.get_node_by_uuid(node_breaker_uuid)
        if node_to_remove:
            system.nodes.remove(node_to_remove)

        #connect branch_connected_with_node_breaker with node_load
        if (self.branch_connected_with_node_breaker.start_node.uuid==self.node_breaker.uuid):
            self.branch_connected_with_node_breaker.start_node = self.node_load
        elif (self.branch_connected_with_node_breaker.end_node.uuid==self.node_breaker.uuid):
            self.branch_connected_with_node_breaker.end_node = self.node_load

        #reenumerate nodes
        system.reindex_nodes_list()


class System():
    def __init__(self):
        self.nodes = []
        self.branches = []
        self.breakers = []
        self.Ymatrix = np.zeros([], dtype=np.complex)
        self.Adjacencies = np.array([])

    def print_nodes_names(self):
        for node in self.nodes:
            print('{} {}'.format(node.name, node.index))

    def print_branchs_names(self):
        for branch in self.branches:
            print('{} {}'.format(branch.start_node.name, branch.start_node.index))
            print('{} {}'.format(branch.end_node.name, branch.end_node.index))
    
    def get_node_by_uuid(self, node_uuid):
        for node in self.nodes:
            if node.uuid == node_uuid:
                return node
        
        return False

    def reindex_nodes_list(self):
        index = 0
        for node in self.nodes:
            node.index=index
            index+=1

    def load_cim_data(self, res, base_apparent_power):
        """
        fill the vectors node and branch
        """
        index = 0
        list_TPNode = [elem for elem in res.values() if elem.__class__.__name__ == "TopologicalNode"]
        list_SvVoltage = [elem for elem in res.values() if elem.__class__.__name__ == "SvVoltage"]
        list_SvPowerFlow = [elem for elem in res.values() if elem.__class__.__name__ == "SvPowerFlow"]
        list_ACLineSegment = [elem for elem in res.values() if elem.__class__.__name__ == "ACLineSegment"]
        list_PowerTransformer = [elem for elem in res.values() if elem.__class__.__name__ == "PowerTransformer"]
        list_Terminals = [elem for elem in res.values() if elem.__class__.__name__ == "Terminal"]
        list_PowerTransformerEnds = [elem for elem in res.values() if elem.__class__.__name__ == "PowerTransformerEnd"]
        list_Breakers = [elem for elem in res.values() if elem.__class__.__name__ == "Breaker"]

        #create nodes
        for TPNode in list_TPNode:
            uuid_TPNode = TPNode.mRID
            name = TPNode.name
            vmag = 0.0
            vphase = 0.0
            pInj = 0.0
            qInj = 0.0
                
            for obj_SvVoltage in list_SvVoltage:
                if obj_SvVoltage.TopologicalNode[0].mRID == uuid_TPNode:
                    vmag = obj_SvVoltage.v
                    vphase = obj_SvVoltage.angle
                    break

            for obj_SvPowerFlow in list_SvPowerFlow:
                if obj_SvPowerFlow.Terminal[0].TopologicalNode[0].mRID == uuid_TPNode:
                    pInj += obj_SvPowerFlow.p
                    qInj += obj_SvPowerFlow.q     
        
            base_voltage = TPNode.BaseVoltage[0].nominalVoltage
            self.nodes.append(Node(name=name, uuid=uuid_TPNode, base_voltage=base_voltage, v_mag=vmag,
                                   base_apparent_power=base_apparent_power, v_phase=vphase,
                                   p=pInj, q=qInj, index=index))
            index = index + 1
        
        self._setNodeType(list_Terminals)   

        #create branches ACLineSegment
        for ACLineSegment in list_ACLineSegment:
            uuid_ACLineSegment = ACLineSegment.mRID
            nodes = self._get_nodes(list_Terminals, uuid_ACLineSegment)
            start_node = nodes[0]
            end_node = nodes[1]

            base_voltage = ACLineSegment.BaseVoltage[0].nominalVoltage
            self.branches.append(Branch(uuid=uuid_ACLineSegment, r=ACLineSegment.r, x=ACLineSegment.x, 
                                        start_node=start_node, end_node=end_node, 
                                        base_voltage=base_voltage, base_apparent_power=base_apparent_power))

        #create branches power transformer
        for power_transformer in list_PowerTransformer:
            uuid_power_transformer = power_transformer.mRID
            nodes = self._get_nodes(list_Terminals, uuid_power_transformer)
            start_node = nodes[0]
            end_node = nodes[1]
            
            # base voltage = high voltage side (=primaryConnection)
            primary_connection = self._get_primary_connection(list_PowerTransformerEnds, uuid_power_transformer)
            base_voltage = primary_connection.BaseVoltage[0].nominalVoltage
            self.branches.append(Branch(uuid=uuid_power_transformer, r=primary_connection.r, x=primary_connection.x,
                                        start_node=start_node, end_node=end_node, base_voltage=base_voltage,
                                        base_apparent_power=base_apparent_power))

        for obj_Breaker in list_Breakers:
            connected = not(obj_Breaker.normalOpen)
            nodes = self._get_nodes(list_Terminals, obj_Breaker.mRID)
            #check which node has loads
            if (nodes[0].power.real==0) and (nodes[0].power.imag == 0):
                node_breaker = nodes[0]
                node_load = nodes[1]
            elif (nodes[1].power.real==0) and (nodes[1].power.imag == 0):
                node_breaker = nodes[0]
                node_load = nodes[1]

            #search for branch_connected_with_node_breaker
            branch_connected_with_node_breaker=None
            for branch in self.branches:
                if branch.start_node.uuid == node_breaker.uuid:
                    branch_connected_with_node_breaker = branch
                    break
                if branch.end_node.uuid == node_breaker.uuid:
                    branch_connected_with_node_breaker = branch
                    break
            
            self.breakers.append(Breaker(node_breaker=node_breaker, node_load=node_load, 
                                         connected=connected, 
                                         branch_connected_with_node_breaker=branch_connected_with_node_breaker))

            #if the breaker.connected == closed --> close broker
            if connected == False:
                self.breakers[-1].close_breaker(self)
            elif connected==True:
                self.breakers[-1].open_breaker(self)

        #calculate admitance matrix
        self.Ymatrix_calc()

    def _get_nodes(self, list_Terminals, elem_uuid):
        """
        get the the start and end node of the element with uuid=elem_uuid
        This function can used only with element which are connected 
        to 2 topologicalNodes, for example: ACLineSegment, PowerTransformer and Breaker 
        :param list_Terminals: list of all elements of type Terminal
        :param elem_uuid: uuid of the element for which the start and end node ID are searched
        :return list: [startNodeID, endNodeID]
        """
        start_node_uuid = None
        end_node_uuid = None
        
        for terminal in list_Terminals:
            #TODO: has always the list ConductingEquipment only one element?
            if (len(terminal.ConductingEquipment)!=1):
                print('WARNING: len(terminal.ConductingEquipment)>1 for the element with uuid={} '.format(elem_uuid))
            conductingEquipment = terminal.ConductingEquipment[0]

            if (conductingEquipment.mRID != elem_uuid):
                continue
            sequence_number = terminal.sequenceNumber
            #TODO: has always the list Terminal.TopologicalNode only one element?
            if (len(terminal.TopologicalNode)!=1):
                print('WARNING: len(terminal.TopologicalNode)>1 for the element with uuid={}'.format(elem_uuid))
            topological_node = terminal.TopologicalNode[0]
            if sequence_number == 1:
                start_node_uuid = topological_node.mRID
            elif sequence_number == 2:
                end_node_uuid = topological_node.mRID
        
        start_node = None
        end_node = None
        if (start_node_uuid==None):
            print('WARNING: It could not find a start node for the element with uuid={}'.format(elem_uuid))
        else:
            start_node = self.get_node_by_uuid(start_node_uuid)
        if (end_node_uuid==None):
            print('WARNING: It could not find a end node for the element with uuid={}'.format(elem_uuid))
        else:
            end_node = self.get_node_by_uuid(end_node_uuid)

        return [start_node, end_node]

    def _get_primary_connection(self, list_PowerTransformerEnds, elem_uuid):
        """
        get primaryConnection of the powertransformer with uuid = elem_uuid
        :param list_PowerTransformerEnd: list of all elements of type PowerTransformerEnd
        :param elem_uuid: uuid of the power transformer for which the primary connection is searched
        :return: primary_connection
        """
        primary_connection = None
        power_transformer_ends = []

        #search for two elements of class powertransformerend that point to the powertransformer with ID = elem_uuid
        for power_transformer_end in list_PowerTransformerEnds:
            if power_transformer_end.PowerTransformer[0].mRID == elem_uuid:
                power_transformer_ends.append(power_transformer_end)

        if power_transformer_ends[0].BaseVoltage[0].nominalVoltage>=power_transformer_ends[1].BaseVoltage[0].nominalVoltage:
            primary_connection=power_transformer_ends[0]
        elif power_transformer_ends[1].BaseVoltage[0].nominalVoltage>=power_transformer_ends[0].BaseVoltage[0].nominalVoltage:        
            primary_connection=power_transformer_ends[1]

        return primary_connection

    def _setNodeType(self, list_Terminals):
        """
        set the parameter type of all elements of the list self.nodes
        :param list_PowerTransformerEnd: list of all elements of type Terminal
        :return None
        """
        #get a list of Terminals for which the ConductingEquipment is a element of class ExternalNetworkInjection
        list_Terminals_ENI = [elem for elem in list_Terminals if elem.ConductingEquipment[0].__class__.__name__ == "ExternalNetworkInjection"]
        for terminal in list_Terminals_ENI:
            #TODO: is it correct to use the element 0 for it? Is len(terminal.TopologicalNode) greater than one?
            node_uuid = terminal.TopologicalNode[0].mRID
            for node in self.nodes:
                if node.uuid == node_uuid:
                    node.type = BusType["SLACK"]
            
        #TODO the search for PV nodes has not been tested yet
        #get a list of Terminals for which the ConductingEquipment is a element of class SynchronousMachine
        list_Terminals_SM = [elem for elem in list_Terminals if elem.ConductingEquipment[0].__class__.__name__ == "SynchronousMachine"]
        for terminal in list_Terminals_SM:
            #TODO: is it correct to use the element 0 for it? Is len(terminal.TopologicalNode) greater than one?
            node_uuid = terminal.TopologicalNode[0].uuid
            for node in self.nodes:
                if node.uuid == node_uuid:
                    node.type = BusType["PV"]

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
            self.Adjacencies[fr].append(to + 1)  # to + 1???
            self.Adjacencies[to].append(fr + 1)  # fr + 1???
