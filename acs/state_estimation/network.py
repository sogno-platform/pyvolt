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
                 v_phase=0.0, p=0.0, q=0.0, index=0, ideal_connected_with=''):
        self.uuid = uuid
        self.name = name
        self.index = index
        self.baseVoltage = base_voltage
        self.base_apparent_power = base_apparent_power
        self.base_current = self.base_apparent_power / self.baseVoltage / np.sqrt(3)
        self.type = BusType["PQ"]
        self.voltage = complex(v_mag * np.cos(np.radians(v_phase)), v_mag * np.sin(np.radians(v_phase)))
        self.power = complex(p, q)
        self.power_pu = complex(p, q) / self.base_apparent_power
        self.voltage_pu = self.voltage / self.baseVoltage
        self.ideal_connected_with = ideal_connected_with

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
    def __init__(self, from_node, to_node, is_open=True):
        """
        :param from_node:
        :param to_node:
        :param is_open: True if the breaker is considered open and False if the broker is closed 
        """
        self.from_node = from_node
        self.to_node = to_node
        self.is_open = is_open

    def __str__(self):
        str = 'class=Breaker\n'
        attributes = self.__dict__
        for key in attributes.keys():
            str = str + key + '={}\n'.format(attributes[key])
        return str

    def open_breaker(self):
        self.is_open == True
        self.to_node.ideal_connected_with = ''

    def close_breaker(self):
        self.is_open == False
        self.to_node.ideal_connected_with = self.from_node.uuid


class System():
    def __init__(self):
        self.nodes = []
        self.branches = []
        self.breakers = []
        self.Ymatrix = np.zeros([], dtype=np.complex)
        self.Adjacencies = np.array([])

    def get_node_by_uuid(self, node_uuid):
        for node in self.nodes:
            if node.uuid == node_uuid:
                return node
        
        return False

    def get_node_by_index(self, index):
        """
        return the node with node.index==index 
        """
        for node in self.nodes:
            if (node.index==index) and (node.ideal_connected_with=='') :
                return node
        
        return None
           
    def get_nodes_num(self):
        """
        return the number of nodes in the list system.nodes
        Warning: if any node is ideally connected to another node, 
        the counter is increased only one time
        """
        nodes_num=0
        for node in self.nodes:
            if node.ideal_connected_with=='':
                nodes_num+=1

        return nodes_num

    def reindex_nodes_list(self):
        """
        Reenumerate the nodes in system.nodes
        If any node is ideally connected to another node, 
        both receive the same index
        """
        index = 0
        remaining_nodes_list = []
        for node in self.nodes:
            if node.ideal_connected_with=='':
                node.index=index
                index+=1
            else:
                remaining_nodes_list.append(node)

        for node in remaining_nodes_list:
            node.index = self.get_node_by_uuid(node.ideal_connected_with).index
             
    def load_cim_data(self, res, base_apparent_power):
        """
        fill the vectors node, branch and breakers
        """
        index = 0
        list_TPNode = [elem for elem in res.values() if elem.__class__.__name__ == "TopologicalNode"]
        list_SvVoltage = [elem for elem in res.values() if elem.__class__.__name__ == "SvVoltage"]
        list_SvPowerFlow = [elem for elem in res.values() if elem.__class__.__name__ == "SvPowerFlow"]
        list_EnergySources = [elem for elem in res.values() if elem.__class__.__name__ == "EnergySource"]
        list_EnergyConsumer = [elem for elem in res.values() if elem.__class__.__name__ == "EnergyConsumer"]
        list_ACLineSegment = [elem for elem in res.values() if elem.__class__.__name__ == "ACLineSegment"]
        list_PowerTransformer = [elem for elem in res.values() if elem.__class__.__name__ == "PowerTransformer"]
        list_Terminals = [elem for elem in res.values() if elem.__class__.__name__ == "Terminal"]
        list_Terminals_ES = [elem for elem in list_Terminals if elem.ConductingEquipment.__class__.__name__ == "EnergySource"]
        list_Terminals_EC = [elem for elem in list_Terminals if elem.ConductingEquipment.__class__.__name__ == "EnergyConsumer"]
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
                if obj_SvVoltage.TopologicalNode.mRID == uuid_TPNode:
                    vmag = obj_SvVoltage.v
                    vphase = obj_SvVoltage.angle
                    break
            for obj_SvPowerFlow in list_SvPowerFlow:
                if obj_SvPowerFlow.Terminal.TopologicalNode.mRID == uuid_TPNode:
                    pInj += obj_SvPowerFlow.p
                    qInj += obj_SvPowerFlow.q           
            for obj_Terminal in list_Terminals_ES:
                if obj_Terminal.TopologicalNode.mRID == uuid_TPNode:
                    for obj_EnergySource in list_EnergySources:
                        if obj_EnergySource.mRID == obj_Terminal.ConductingEquipment.mRID:
                            pInj += obj_EnergySource.activePower
                            qInj += obj_EnergySource.reactivePower
            for obj_Terminal in list_Terminals_EC:
                if obj_Terminal.TopologicalNode.mRID == uuid_TPNode:
                    for obj_EnergyConsumer in list_EnergyConsumer:
                        if obj_EnergyConsumer.mRID == obj_Terminal.ConductingEquipment.mRID:
                            pInj += obj_EnergyConsumer.p
                            qInj += obj_EnergyConsumer.q
            
            base_voltage = TPNode.BaseVoltage.nominalVoltage
            self.nodes.append(Node(name=name, uuid=uuid_TPNode, base_voltage=base_voltage, v_mag=vmag,
                                   base_apparent_power=base_apparent_power, v_phase=vphase,
                                   p=pInj, q=qInj, index=index))
            index = index + 1
        
        self._setNodeType(list_Terminals)   

        #create branches type ACLineSegment
        for ACLineSegment in list_ACLineSegment:
            uuid_ACLineSegment = ACLineSegment.mRID
            nodes = self._get_nodes(list_Terminals, uuid_ACLineSegment)
            start_node = nodes[0]
            end_node = nodes[1]

            base_voltage = ACLineSegment.BaseVoltage.nominalVoltage
            self.branches.append(Branch(uuid=uuid_ACLineSegment, r=ACLineSegment.r, x=ACLineSegment.x, 
                                        start_node=start_node, end_node=end_node, 
                                        base_voltage=base_voltage, base_apparent_power=base_apparent_power))

        #create branches type powerTransformer
        for power_transformer in list_PowerTransformer:
            uuid_power_transformer = power_transformer.mRID
            nodes = self._get_nodes(list_Terminals, uuid_power_transformer)
            start_node = nodes[0]
            end_node = nodes[1]
            
            # base voltage = high voltage side (=primaryConnection)
            primary_connection = self._get_primary_connection(list_PowerTransformerEnds, uuid_power_transformer)
            base_voltage = primary_connection.BaseVoltage.nominalVoltage
            self.branches.append(Branch(uuid=uuid_power_transformer, r=primary_connection.r, x=primary_connection.x,
                                        start_node=start_node, end_node=end_node, base_voltage=base_voltage,
                                        base_apparent_power=base_apparent_power))

        #create breakers
        for obj_Breaker in list_Breakers:
            is_open = obj_Breaker.normalOpen
            nodes = self._get_nodes(list_Terminals, obj_Breaker.mRID)
            self.breakers.append(Breaker(from_node=nodes[0], to_node=nodes[1], is_open=is_open))

            #if the breaker is open == closed --> close broker
            if is_open == False:
                self.breakers[-1].close_breaker(self)
            else:
                self.breakers[-1].ideal_connected_with = ''
            
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
            if (terminal.ConductingEquipment.mRID != elem_uuid):
                continue
            sequence_number = terminal.sequenceNumber            
            if sequence_number == 1:
                start_node_uuid = terminal.TopologicalNode.mRID
            elif sequence_number == 2:
                end_node_uuid = terminal.TopologicalNode.mRID
        
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
            power_transformer = None
            if isinstance(power_transformer_end.PowerTransformer, list):
                if (len(power_transformer_end.PowerTransformer)!=1):
                    print('WARNING: len(power_transformer_end.PowerTransformer)!=1 for the element with uuid={}. \
                        The first element will be used'.format(power_transformer_end.mRID))
                power_transformer = power_transformer_end.PowerTransformer[0]
            else:
                power_transformer = power_transformer_end.PowerTransformer
        
            if power_transformer.mRID == elem_uuid:
                power_transformer_ends.append(power_transformer_end)

        if power_transformer_ends[0].BaseVoltage.nominalVoltage>=power_transformer_ends[1].BaseVoltage.nominalVoltage:
            primary_connection=power_transformer_ends[0]
        elif power_transformer_ends[1].BaseVoltage.nominalVoltage>=power_transformer_ends[0].BaseVoltage.nominalVoltage:
            primary_connection=power_transformer_ends[1]

        return primary_connection

    def _setNodeType(self, list_Terminals):
        """
        set the parameter "type" of all elements of the list self.nodes
        :param list_PowerTransformerEnd: list of all elements of type Terminal
        :return None
        """
        #get a list of Terminals for which the ConductingEquipment is a element of class ExternalNetworkInjection
        list_Terminals_ENI = [elem for elem in list_Terminals if elem.ConductingEquipment.__class__.__name__ == "ExternalNetworkInjection"]
        for terminal in list_Terminals_ENI:
            node_uuid = terminal.TopologicalNode.mRID
            for node in self.nodes:
                if node.uuid == node_uuid:
                    node.type = BusType["SLACK"]
            
        #TODO the search for PV nodes has not been tested yet
        #get a list of Terminals for which the ConductingEquipment is a element of class SynchronousMachine
        list_Terminals_SM = [elem for elem in list_Terminals if elem.ConductingEquipment.__class__.__name__ == "SynchronousMachine"]
        for terminal in list_Terminals_SM:
            node_uuid = terminal.TopologicalNode.mRID
            for node in self.nodes:
                if node.uuid == node_uuid:
                    node.type = BusType["PV"]

    def Ymatrix_calc(self):
        self.reindex_nodes_list()
        nodes_num = self.get_nodes_num()
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
    
    #Testing functions
    def print_nodes_names(self):
        for node in self.nodes:
            print('{} {}'.format(node.name, node.index))

    def print_node_types(self):
        for node in self.nodes:
            print('{} {}'.format(node.name, node.type))
    
    def print_power(self):
        for node in self.nodes:
            print('{} {}'.format(node.name, node.power))
            
