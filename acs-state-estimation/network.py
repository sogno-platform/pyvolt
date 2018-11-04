class Node():
	def __init__(self, name, uuid, pInjection, qInjection, index):
		self.name=''
		self.uuid=''
		self.pInjection=0.0
		self.qInjection=0.0
		self.index=None
		
class Branch():
	def __init__(self, r, x, bstart, bend):
		self.r=0.0
		self.x=0.0
		self.startNode=bstart
		self.endNode=bend

class System():
	def __init__(self, res):
		self.node=[]
		self.branch=[]
		self.bR=[]
		self.bX=[]
		self.P=[]
		self.Q=[]
	
		#this function is used to fill the vectors node, branch, bR, bX, P and Q
		for key, value in res.items():
			if value.__class__.__name__=="TopologicalNode":
				self.P.append(value.pInjection)
				self.Q.append(value.qInjection)
				index=len(self.P)-1
				self.node.append(Node(value.name, value.mRID, value.pInjection, value.qInjection, index))
			elif value.__class__.__name__=="ACLineSegment":
				length=value.length
				if length==0.0:
					length=1.0
				self.bR.append(value.r*length)
				self.bX.append(value.x*length)
				startNode=res[value.startNodeID]
				endNode=res[value.endNodeID]
				self.branch.append(Branch(value.r, value.x, startNode, endNode))	
			elif value.__class__.__name__=="PowerTransformer":
				self.bR.append(value.primaryConnection.r)
				self.bX.append(value.primaryConnection.x)
				startNode=res[value.startNodeID]
				endNode=res[value.endNodeID]
				self.branch.append(Branch(value.primaryConnection.r, value.primaryConnection.x, startNode, endNode))
			else:
				continue
	

import sys
sys.path.append("..\..\cimpy")
import cimpy
import logging
logging.basicConfig(level=logging.INFO)

xml_files=[r"..\..\cim-grid-data\WSCC-09\WSCC-09_Neplan_EQ.xml", 
		   r"..\..\cim-grid-data\WSCC-09\WSCC-09_Neplan_SV.xml",
		   r"..\..\cim-grid-data\WSCC-09\WSCC-09_Neplan_TP.xml"]

res=cimpy.cimread(xml_files)
cimpy.setNodes(res)
cimpy.setPowerTransformerEnd(res)
network=System(res)

print("Vector bR:")
print(network.bR)
print()
print("Vector bX:")
print(network.bX)
print()
print("Vector P:")
print(network.P)
print()
print("Vector Q:")
print(network.Q)
print()
