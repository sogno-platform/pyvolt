#!/usr/bin/env python

################
# Dependancies:
# paho, mqtt library for python (sudo pip install paho-mqtt)
import sys
sys.path.append(r"C:\Users\acs-opal\git\state-estimation\acs\state_estimation")
sys.path.append(r"C:\Users\acs-opal\git\cimpy")
sys.path.append(r"C:\Users\acs-opal\git")

import paho.mqtt.client as mqtt
import json
import numpy as np
import cimpy
import network
import nv_state_estimator_cim
import measurement_generator_online
import mqtt_client


#class PerUnit:
#    def __init__(self, S, V):
#        self.S = S
#        self.V = V
#        self.I = S/V
#        self.Z = S/(V**2)

def on_message(mqttc, userdata, msg):
    message = json.loads(msg.payload)[0]
    #print("Ecco il messaggio: ",message["data"])
    #flag = 1
    values = message["data"]
    
#    for elem in range(len(system.nodes)): 
#        values[elem] = values[system.nodes[elem].uuid].values[0]
    
    zdatameas =  measurement_generator_online.Zdatameas_creation_fromPF(meas_config, zdata, values)
    Vest, Iest, Iinjest, S1est, S2est, Sinjest = nv_state_estimator_cim.DsseCall(system, zdatameas, Ymatr, Adj)
    print("AAAAAAAAAAAAAA",Vest.mag)
    num_data = len(values)
    index = np.linspace(0,num_data-2,int(num_data/2)).astype(int)
    payload["ts"]["origin"] = message["ts"]["origin"]
    payload["sequence"] = message["sequence"]
    payload["data"] = np.append(values,values)
    payload["data"][index] = Vest.mag
    payload["data"][index+1] = Vest.phase
    #mqttc.publish("sogno-estimator",json.dumps(payload),0)
    

xml_files=[r"C:\Users\acs-opal\git\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_EQ.xml", 
		   r"C:\Users\acs-opal\git\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_SV.xml",
		   r"C:\Users\acs-opal\git\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_TP.xml"]


res=cimpy.cimread(xml_files)
cimpy.setNodes(res)
cimpy.setPowerTransformerEnd(res)
system = network.System()
system.load_cim_data(res)


Ymatr, Adj = network.Ymatrix_calc(system)


with open('Payload_config.json') as f1:
    payload = json.load(f1)

with open('Measurement_config.json') as f2:
    meas_config = json.load(f2)
zdata = measurement_generator_online.Zdata_structure_creation(meas_config, system)

mqttc = mqtt_client.connect_subscribe("dpsim-powerflow")
mqttc.on_message = on_message
#flag = 0

#while 1:
#    if flag == 1:
#        values = message["data"]
#        zdatameas =  measurement_generator_online.Zdatameas_creation_fromPF(meas_config, zdatameas, values)
#        Vest, Iest, Iinjest, S1est, S2est, Sinjest = nv_state_estimator_cim.DsseCall(system, zdatameas, Ymatr, Adj)
