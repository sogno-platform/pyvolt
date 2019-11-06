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
    def __init__(self, uuid='', base_voltage=1.0, base_apparent_power=1.0, v_mag=0.0,
                 v_phase=0.0, p=0.0, q=0.0, index=0, bus_type='PQ', pu=False):
        self.uuid = uuid
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

class System():
    def __init__(self):
        self.nodes = []
        self.branches = []
        self.Ymatrix = np.zeros([], dtype=np.complex)
        self.Adjacencies = np.array([])

    def load_cim_data(self, res, base_apparent_power):
        """
        To fill the vectors node, branch, bR, bX, P and Q
        """
        index = 0
        list_TPNode = [elem for elem in res.values() if elem.__class__.__name__ == "TopologicalNode"]
        list_SvVoltage = [elem for elem in res.values() if elem.__class__.__name__ == "SvVoltage"]
        list_SvPowerFlow = [elem for elem in res.values() if elem.__class__.__name__ == "SvPowerFlow"]
        #list_EnergyConsumer = [elem for elem in res.values() if elem.__class__.__name__ == "EnergyConsumer"]
        list_ACLineSegment = [elem for elem in res.values() if elem.__class__.__name__ == "ACLineSegment"]
        list_PowerTransformer = [elem for elem in res.values() if elem.__class__.__name__ == "PowerTransformer"]
        list_Terminals = [elem for elem in res.values() if elem.__class__.__name__ == "Terminal"]
        #list_Terminals_EC = [elem for elem in list_Terminals if elem.ConductingEquipment[0].__class__.__name__ == "EnergyConsumer"]
        list_PowerTransformerEnds = [elem for elem in res.values() if elem.__class__.__name__ == "PowerTransformerEnd"]

        #create nodes
        for TPNode in list_TPNode:
            name_TPNode = TPNode.name
            vmag = 0.0
            vphase = 0.0
            pInj = 0.0
            qinj = 0.0
                
            for obj_SvVoltage in list_SvVoltage:
                if obj_SvVoltage.TopologicalNode[0].name == name_TPNode:
                    vmag = obj_SvVoltage.v
                    vphase = obj_SvVoltage.angle
                    break

            for obj_SvPowerFlow in list_SvPowerFlow:
                #TODO: gibt es nur ein SvPowerFlow pro Knoten??
                if obj_SvPowerFlow.Terminal[0].TopologicalNode[0].name == name_TPNode:
                    pInj += obj_SvPowerFlow.p
                    qinj += obj_SvPowerFlow.q                
        
            #for obj_Terminal in list_Terminals_EC:
            #    #TODO: is it correct?
            #    if obj_Terminal.TopologicalNode[0].name == name_TPNode:
            #        energy_consumer_name = obj_Terminal.ConductingEquipment[0].name


            base_voltage = TPNode.BaseVoltage[0].nominalVoltage
            self.nodes.append(Node(uuid=name_TPNode, base_voltage=base_voltage, v_mag=vmag,
                                   base_apparent_power=base_apparent_power, v_phase=vphase,
                                   p=pInj, q=qinj, index=index))
            index = index + 1
        
        #create branches ACLineSegment
        for ACLineSegment in list_ACLineSegment:
            name_ACLineSegment = ACLineSegment.name
            start_node_id, end_node_id = self._get_node_names(list_Terminals, name_ACLineSegment)
            start_node = None
            end_node = None

            for node in self.nodes:
                if start_node_id == node.uuid:
                    start_node = node
                elif end_node_id == node.uuid:
                    end_node = node

            base_voltage = ACLineSegment.BaseVoltage[0].nominalVoltage
            self.branches.append(Branch(uuid=name_ACLineSegment, r=ACLineSegment.r, x=ACLineSegment.x, 
                                        start_node=start_node, end_node=end_node, 
                                        base_voltage=base_voltage, base_apparent_power=base_apparent_power))

        #create branches power transformer
        for power_transformer in list_PowerTransformer:
            name_power_transformer = power_transformer.name
            start_node_id, end_node_id = self._get_node_names(list_Terminals, name_power_transformer)
            start_node = None
            end_node = None

            for node in self.nodes:
                if start_node_id == node.uuid:
                    start_node = node
                elif end_node_id == node.uuid:
                    end_node = node

            # base voltage = high voltage side (=primaryConnection)
            primary_connection = self._get_primary_connection(list_PowerTransformerEnds, name_power_transformer)
            base_voltage = primary_connection.BaseVoltage[0].nominalVoltage
            self.branches.append(Branch(uuid=name_power_transformer, r=primary_connection.r, x=primary_connection.x,
                                        start_node=start_node, end_node=end_node, base_voltage=base_voltage,
                                        base_apparent_power=base_apparent_power))

        self._setNodeType(list_Terminals)
        # calculate admitance matrix
        self.Ymatrix_calc()

    def _get_node_names(self, list_Terminals, elem_name):
        """
        get the names of the start and end nodes of the element with name=elem_name
        This function can only used with the elements ACLineSegment, PowerTransformer, Switch and EnergyConsumer
        :param list_Terminals: list of all elements of type Terminal
        :param elem_name: name of the element for which the start and end node ID are searched
        :return tuple: (startNodeID, endNodeID)
        """
        start_node_name = ""
        end_node_name = ""
        list_elements = ["ACLineSegment", "PowerTransformer", "Switch", "EnergyConsumer"]
        
        for terminal in list_Terminals:
            #TODO: has always the list ConductingEquipment only one element?????
            conductingEquipment = terminal.ConductingEquipment[0]
            if (conductingEquipment.name != elem_name):
                continue
            elif conductingEquipment.__class__.__name__ in list_elements:
                sequenceNumber = terminal.sequenceNumber
                #TODO: has always the list TopologicalNode only one element?????
                topologicalNode = terminal.TopologicalNode[0]
                if sequenceNumber == 1:
                    start_node_name = topologicalNode.name
                elif sequenceNumber == 2:
                    end_node_name = topologicalNode.name
        
        return (start_node_name, end_node_name)

    def _get_primary_connection(self, list_PowerTransformerEnds, elem_name):
        """
        get primaryConnection of the powertransformer with name = elem_name
        :param list_PowerTransformerEnd: list of all elements of type PowerTransformerEnd
        :param elem_name: name of the power transformer for which the primary connection is searched
        :return: primary_connection
	    """
        primary_connection = None
        power_transformer_ends = []

        #search for two elements of class powertransformerend that point to the powertransformer with ID = elem_id
        for power_transformer_end in list_PowerTransformerEnds:
            if power_transformer_end.PowerTransformer[0].name == elem_name:
                power_transformer_ends.append(power_transformer_end)

        #TODO: ERROR HANDLING
        if len(power_transformer_ends) !=2:
            print("ERROR!!!")
            print(len(power_transformer_ends))
            return -1

        if power_transformer_ends[0].BaseVoltage[0].nominalVoltage>=power_transformer_ends[1].BaseVoltage[0].nominalVoltage:
            primary_connection=power_transformer_ends[0]
        elif power_transformer_ends[1].BaseVoltage[0].nominalVoltage>=power_transformer_ends[0].BaseVoltage[0].nominalVoltage:		
            primary_connection=power_transformer_ends[1]

        return primary_connection

    def _setNodeType(self, list_Terminals):
        """
        return the type of a node: PQ, PV or SLACK

        @param node: element of class cimpy.Topology.TopologicalNode
        """
        #get a list of Terminals for which the ConductingEquipment is a element of class ExternalNetworkInjection
        list_Terminals_ENI = [elem for elem in list_Terminals if elem.ConductingEquipment[0].__class__.__name__ == "ExternalNetworkInjection"]
        for terminal in list_Terminals_ENI:
            #TODO: is it correct to use the element 0 for it? Is len(terminal.TopologicalNode) greater than one?
            node_name = terminal.TopologicalNode[0].name
            for node in self.nodes:
                if node.uuid == node_name:
                    node.type = BusType["SLACK"]
            
        #get a list of Terminals for which the ConductingEquipment is a element of class ExternalNetworkInjection
        list_Terminals_SM = [elem for elem in list_Terminals if elem.ConductingEquipment[0].__class__.__name__ == "SynchronousMachine"]
        for terminal in list_Terminals_SM:
            #TODO: is it correct to use the element 0 for it? Is len(terminal.TopologicalNode) greater than one?
            node_name = terminal.TopologicalNode[0].name
            for node in self.nodes:
                if node.uuid == node_name:
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
