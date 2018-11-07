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