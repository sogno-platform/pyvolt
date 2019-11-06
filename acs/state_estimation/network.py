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
        self.type = BusType[bus_type]
        self.voltage = v_mag * np.cos(v_phase) + 1j * v_mag * np.sin(v_phase)
        self.power = complex(p, q)
        self.power_pu = complex(p, q) / self.base_apparent_power
        self.voltage_pu = self.voltage / self.baseVoltage


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
        for uuid_TPNode, TPNode in res.items():
            if TPNode.__class__.__name__ == "TopologicalNode":
                vmag = 0.0
                vphase = 0.0
                pInj = 0.0
                qinj = 0.0
                for uuid_obj_SvVoltage, obj_SvVoltage in res.items():
                    if obj_SvVoltage.__class__.__name__ == "SvVoltage":
                        if obj_SvVoltage.TopologicalNode[0].name == uuid_TPNode:
                            vmag = obj_SvVoltage.v
                            vphase = obj_SvVoltage.angle
                            break
                for uuid_obj_SvPowerFlow, obj_SvPowerFlow in res.items():
                    if obj_SvPowerFlow.__class__.__name__ == "SvPowerFlow":
                        if obj_SvPowerFlow.Terminal[0].TopologicalNode[0].name == uuid_TPNode:
                            for node in obj_SvPowerFlow.Terminal[0].TopologicalNode:
                                pInj += obj_SvPowerFlow.p
                                qinj += obj_SvPowerFlow.q
                            break
                node_type = self._getNodeType(TPNode)
                base_voltage = TPNode.BaseVoltage[0].nominalVoltage
                self.nodes.append(Node(uuid=uuid_TPNode, base_voltage=base_voltage, v_mag=vmag,
                                       base_apparent_power=base_apparent_power, v_phase=vphase,
                                       p=pInj, q=qinj, index=index, bus_type=node_type))
                index = index + 1
        
        for uuid, element in res.items():
            if (element.__class__.__name__ == "ACLineSegment") or (element.__class__.__name__ == "PowerTransformer"):
                start_node_id, end_node_id = self._get_NodeIDs(res, element.name)
                start_node = None
                end_node = None

                for node in self.nodes:
                    if start_node_id == node.uuid:
                        start_node = node
                    elif end_node_id == node.uuid:
                        end_node = node

                if element.__class__.__name__ == "ACLineSegment":
                    base_voltage = element.BaseVoltage[0].nominalVoltage
                    self.branches.append(Branch(uuid=uuid, r=element.r, x=element.x, start_node=start_node,
                                                end_node=end_node, base_voltage=base_voltage,
                                                base_apparent_power=base_apparent_power))
                elif element.__class__.__name__ == "PowerTransformer":
                    # base voltage = high voltage side (=primaryConnection)
                    primary_connection = self._get_primary_connection(res, element.name)
                    base_voltage = primary_connection.BaseVoltage[0].nominalVoltage
                    self.branches.append(Branch(uuid=uuid, r=primary_connection.r, x=primary_connection.x,
                                                start_node=start_node, end_node=end_node, base_voltage=base_voltage,
                                                base_apparent_power=base_apparent_power))

        # calculate admitance matrix
        self.Ymatrix_calc()

    def _get_NodeIDs(self, start_dict, elem_name):
        """
        get the startNodeId and endNodeID of the element with ID=mRID
        This function can only used with the elements ACLineSegment, PowerTransformer, Switch and EnergyConsumer
        :param start_dict: result of cimpy.cim_import
        :param elem_name: 
        :return tuple: (startNodeID, endNodeID)
        """
        startNodeID = ""
        endNodeID = ""
        list_elements = ["ACLineSegment", "PowerTransformer", "Switch", "EnergyConsumer"]
        
        #search for elements type 'Terminal'
        for value, key in start_dict.items():
            if key.__class__.__name__=='Terminal':
                #TODO: has always the list ConductingEquipment only one element?????
                conductingEquipment = key.ConductingEquipment[0]
                if (conductingEquipment.name != elem_name):
                    continue
                elif conductingEquipment.__class__.__name__ in list_elements:
                    sequenceNumber = key.sequenceNumber
                    #TODO: has always the list TopologicalNode only one element?????
                    TopologicalNode = key.TopologicalNode[0]
                    if sequenceNumber == 1:
                        startNodeID = TopologicalNode.name
                    elif sequenceNumber == 2:
                        endNodeID = TopologicalNode.name
        
        return (startNodeID, endNodeID)

    def _get_primary_connection(self, start_dict, elem_name):
        """
        get primaryConnection of the powertransformer with ID = elem_id
        :param start_dict: result of cimpy.cim_import
        :param elem_name: 
        :return tuple: (primary_connection, secondary_connection)
	    """
        primary_connection = None
        power_transformer_end = []
        for uuid, element in start_dict.items():
            #search for two elements of class powertransformerend that point to the powertransformer with ID = elem_id
            if element.__class__.__name__=='PowerTransformerEnd':
                if element.PowerTransformer[0].name == elem_name:
                    power_transformer_end.append(element)

        #TODO: ERROR HANDLING
        if len(power_transformer_end) !=2:
            print("ERROR!!!")
            print(len(power_transformer_end))
            return -1

        if power_transformer_end[0].BaseVoltage[0].nominalVoltage>=power_transformer_end[1].BaseVoltage[0].nominalVoltage:
            primary_connection=power_transformer_end[0]
        elif power_transformer_end[1].BaseVoltage[0].nominalVoltage>=power_transformer_end[0].BaseVoltage[0].nominalVoltage:		
            primary_connection=power_transformer_end[1]

        return primary_connection

    def _getNodeType(self, node):
        """
        return the type of a node: PQ, PV or SLACK

        @param node: element of class cimpy.Topology.TopologicalNode
        """
        for terminal in node.Terminal:
            if terminal.ConductingEquipment.__class__.__name__ == "ExternalNetworkInjection":
                return "SLACK"
            elif terminal.ConductingEquipment.__class__.__name__ == "SynchronousMachine":
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
            self.Adjacencies[fr].append(to + 1)  # to + 1???
            self.Adjacencies[to].append(fr + 1)  # fr + 1???
