#!/usr/bin/env python

################
# Dependancies:
# paho, mqtt library for python (sudo pip install paho-mqtt)
import paho.mqtt.client as mqtt
import json
import numpy as np
import time
import logging

logging.basicConfig(filename='recv_client.log', level=logging.INFO)

import sys
#add to path the path to the state-estimation commit 23facb6b60d3f7993eb05ba7703f137b588343ec
sys.path.append("../../../state-estimation-client/state-estimation/acs/state_estimation")
from network import System, Ymatrix_calc
from nv_state_estimator_cim import DsseCall
from nv_state_estimator_cim import Complex_to_all
from nv_state_estimator_cim import Real_to_all
from measurement_generator_online import Zdata_structure_creation, Zdatameas_creation_fromPF
import cimpy

Sb = 25e6
Vb1 = 110e3
Vb2 = 20e3
Zb = (Vb2**2)/Sb

def on_message(mqttc, userdata, msg):
    global end
    #stop simulation when message "end" is received
    if msg.payload=="end":
        end = True
    else:   
        message = json.loads(msg.payload)[0]
        values = message["data"]
        values_se = np.zeros(len(values))
        Vmag_est = np.zeros(len(system.nodes))
        Vphase_est = np.zeros(len(system.nodes))
        Vmag_true = np.zeros(len(system.nodes))
        for elem in range(len(system.nodes)):
            #if int(system.nodes[elem].uuid[1:]) == 0:
            if system.nodes[elem].uuid == "N0":
                values_se[elem*2] = values[int(system.nodes[elem].uuid[1:])*2]/Vb1
            else:
                values_se[elem * 2] = values[int(system.nodes[elem].uuid[1:])*2]/Vb2
            values_se[elem*2+1] = values[int(system.nodes[elem].uuid[1:])*2+1]
    
        if int(message["sequence"]) > 0:
            # To delete later
            nodes_num = len(system.nodes)
            branches_num = len(system.branches)
            Vt = np.zeros(len(system.nodes), dtype=np.complex)
            #index = np.linspace(0, len(values) - 2, int(len(values) / 2)).astype(int)     
            for elem in range(len(system.nodes)):
                index = 2*system.nodes[elem].index
                #Vt[elem] = values_se[index[elem]]*(np.cos(values_se[index[elem]+1]) + 1j*np.sin(values_se[index[elem]+1]))
                Vt[elem] = values_se[index]*(np.cos(values_se[index+1]) + 1j*np.sin(values_se[index+1]))
            Vtrue = Complex_to_all(Vt)
            #print(Vtrue.mag)

            """ From here on, all the other quantities of the grid are calculated """
            Irx = np.zeros(branches_num, dtype=np.complex)
            for index in range(branches_num):
                fr = system.branches[index].start_node.index  # branch.start[idx]-1
                to = system.branches[index].end_node.index  # branch.end[idx]-1
                Irx[index] = - (Vtrue.complex[fr] - Vtrue.complex[to]) * Ymatr[fr][to]
            Ir = np.real(Irx)
            Ix = np.imag(Irx)

            Itrue = Real_to_all(Ir, Ix)
            Iinj_r = np.zeros(nodes_num)
            Iinj_x = np.zeros(nodes_num)
            for k in range(nodes_num):
                to = []
                fr = []
                for m in range(branches_num):
                    if k == system.branches[m].start_node.index:
                        fr.append(m)
                    if k == system.branches[m].end_node.index:
                        to.append(m)

                Iinj_r[k] = np.sum(Itrue.real[to]) - np.sum(Itrue.real[fr])
                Iinj_x[k] = np.sum(Itrue.imag[to]) - np.sum(Itrue.imag[fr])

            Iinjtrue = Real_to_all(Iinj_r, Iinj_x)
            Sinj_rx = np.multiply(Vtrue.complex, np.conj(Iinjtrue.complex))
            Sinjtrue = Real_to_all(np.real(Sinj_rx), np.imag(Sinj_rx))
            values_se = np.append(values_se, Sinjtrue.real)
            values_se = np.append(values_se, Sinjtrue.imag)

            # till here

            if message["sequence"] < 90:
                zdatameas = Zdatameas_creation_fromPF(meas_config1, zdata1, values_se, message["sequence"])
            else:
                zdatameas = Zdatameas_creation_fromPF(meas_config2, zdata2, values_se, message["sequence"])

            Vest, Iest, Iinjest, S1est, S2est, Sinjest = DsseCall(system, zdatameas, Ymatr, Adj)
            print("state_estimation_res finish!")
            
            Vmag_est = []
            Vphase_est = []
            for elem in range(len(system.nodes)):
                Vmag_est.append(Vest.mag[elem])
                Vphase_est.append(Vest.phase[elem])

            #send results to dummy_simulator
            payload = {}
            payload["client"] = "Sogno_cigre_se"
            payload["sequence"] = message["sequence"]
            payload["V_est_mag"] = Vmag_est
            payload["V_est_phase"] = Vphase_est

            mqttc.publish(topic_dummy_simulator, "[" + json.dumps(payload) + "]", 0)

def on_connect(mqttc, userdata, flags, rc):
    if rc == 0:
        mqttc.connected_flag = False  # set flag
        print("connected OK with returned code=", rc)
    else:
        print("Bad connection, error code= ", rc)

xml_files = [
    r"..\..\..\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_EQ.xml",
    r"..\..\..\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_SV.xml",
    r"..\..\..\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_TP.xml"]

res = cimpy.cimread(xml_files)
cimpy.setNodes(res)
cimpy.setPowerTransformerEnd(res)
system = System()
system.load_cim_data(res)
system.branches[12].y = system.branches[12].y*(Vb1**2)/(Vb2**2)
system.branches[13].y = system.branches[13].y*(Vb1**2)/(Vb2**2)
Ymatrix, Adj = Ymatrix_calc(system)
Ymatr = Ymatrix*Zb

meas_configfile = r"C:\Users\Martin\Desktop\hiwi\git\state-estimation-client\sogno_demo\Measurement_config.json"
meas_configfile2 = r"C:\Users\Martin\Desktop\hiwi\git\state-estimation-client\sogno_demo\Measurement_config2.json"
meas_configfile3 = r"C:\Users\Martin\Desktop\hiwi\git\state-estimation-client\sogno_demo\Measurement_config3.json"

with open(meas_configfile) as f1:
    payload = json.load(f1)

with open(meas_configfile2) as f2:
    meas_config1 = json.load(f2)

with open(meas_configfile3) as f3:
    meas_config2 = json.load(f3)

zdata1 = Zdata_structure_creation(meas_config1, system)
zdata2 = Zdata_structure_creation(meas_config2, system)

#global variable used in callback on_message
end = False

#topic names
topic = "dpsim-powerflow"
topic_dummy_simulator = "dpsim-powerflow_results"

broker_adress = "m16.cloudmqtt.com"
mqtt.Client.connected_flag=False                        #create flag in class
mqttc = mqtt.Client("SognoDemo_Client2", True)          #create new instance
mqttc.username_pw_set("ilgtdaqk", "UbNQQjmcUdqq")
mqttc.on_connect = on_connect                           #attach function to callback
mqttc.on_message = on_message                           #attach function to callback
mqttc.connect(broker_adress, 14543)	                    #connect to broker
mqttc.loop_start()                                      #start loop to process callback
time.sleep(4)                                           #wait for connection setup to complete
mqttc.subscribe(topic)

#send notification message notifying that the connection was successful
mqttc.publish(topic_dummy_simulator, json.dumps("OK"))

while not end:
    time.sleep(0.01)

mqttc.loop_stop()       #Stop loop 
mqttc.disconnect()      #disconnect