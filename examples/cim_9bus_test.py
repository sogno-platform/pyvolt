import sys
import logging

sys.path.append("..\..\cimpy")
sys.path.append("..")
import cimpy
from acs.state_estimation.network import System  

logging.basicConfig(level=logging.INFO)

xml_files=[r"..\..\cim-grid-data\WSCC-09\WSCC-09_Neplan_EQ.xml", 
		   r"..\..\cim-grid-data\WSCC-09\WSCC-09_Neplan_SV.xml",
		   r"..\..\cim-grid-data\WSCC-09\WSCC-09_Neplan_TP.xml"]

res=cimpy.cimread(xml_files)
cimpy.setNodes(res)
cimpy.setPowerTransformerEnd(res)
network=System()
network.load_cim_data(res)

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