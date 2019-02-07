#!/usr/bin/env python

################
# Dependancies:
# paho, mqtt library for python (sudo pip install paho-mqtt)
import sys
sys.path.append("D:/Sciebo_folder/Temporary_folder/SOGNO/Demo/acs/state_estimation")
sys.path.append("D:/Sciebo_folder/Temporary_folder/SOGNO/Demo/cimpy")

import paho.mqtt.client as mqtt
import time
import json
import numpy as np
import math
import copy
import cimpy
import network
import nv_state_estimator_cim
import measurement_generator_online
import mqtt_client


class PerUnit:
    def __init__(self, S, V):
        self.S = S
        self.V = V
        self.I = S/V
        self.Z = S/(V**2)

def on_message(mqttc, userdata, msg):
    message = json.loads(msg.payload)
    print("Ecco il messaggio: ",message["data"])
    flag = 1
    values = message["data"]
    zdatameas =  measurement_generator_online.Zdatameas_creation_fromPF(meas_config, zdata, values)
    Vest, Iest, Iinjest, S1est, S2est, Sinjest = nv_state_estimator_cim.DsseCall(system, zdatameas, Ymatr, Adj)
    index = np.linspace(0,num_data-2,num_data/2)
    payload["data"] = np.append(values,values)
    payload["data"][index] = Vest.mag
    payload["data"][index+1] = Vest.phase
    mqttc.publish("SEresults",json.dumps(payload),0)
    

xml_files=[r"D:\Sciebo_folder\Temporary_folder\SOGNO\Demo\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_EQ.xml", 
		   r"D:\Sciebo_folder\Temporary_folder\SOGNO\Demo\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_SV.xml",
		   r"D:\Sciebo_folder\Temporary_folder\SOGNO\Demo\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_TP.xml"]


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

mqttc = mqtt_client.connect_subscribe("DPsimResults")
mqttc.on_message = on_message
flag = 0

#while 1:
#    if flag == 1:
#        values = message["data"]
#        zdatameas =  measurement_generator_online.Zdatameas_creation_fromPF(meas_config, zdatameas, values)
#        Vest, Iest, Iinjest, S1est, S2est, Sinjest = nv_state_estimator_cim.DsseCall(system, zdatameas, Ymatr, Adj)
