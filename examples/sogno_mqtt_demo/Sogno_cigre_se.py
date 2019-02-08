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
import time
import cimpy
import network
import nv_state_estimator_cim
import measurement_generator_online
import mqtt_client


def on_message(mqttc, userdata, msg):
    message = json.loads(msg.payload)[0]
    values = message["data"]

    if message["sequence"] > 0:
        zdatameas = measurement_generator_online.Zdatameas_creation_fromPF(meas_config, zdata, values)
        Vest, Iest, Iinjest, S1est, S2est, Sinjest = nv_state_estimator_cim.DsseCall(system, zdatameas, Ymatr, Adj)
        index = np.linspace(0, len(values) - 2, int(len(values) / 2)).astype(int)

        payload["ts"]["origin"] = message["ts"]["origin"]
        payload["sequence"] = message["sequence"]
        payload["data"] = np.append(values, values)
        payload["data"][index] = Vest.mag
        payload["data"][index + 1] = Vest.phase
        payload["data"] = list(payload["data"])

        mqttc.publish("sogno-estimator", "[" + json.dumps(payload) + "]", 0)


def on_log(mqttc, userdata, level, buf):
    print("log: " + buf)


def on_connect(mqttc, userdata, flags, rc):
    if rc == 0:
        mqttc.connected_flag = False  # set flag
        print("Publisher connection OK")
    else:
        print("Bad connection, error code= ", rc)


def connect_subscribe(topic):
    mqttc = mqtt.Client("SognoSE", True)
    mqttc.username_pw_set("villas", "s3c0sim4!")
    mqttc.on_log = on_log
    mqttc.on_connect = on_connect

    mqttc.connect("137.226.248.91")

    mqttc.loop_start()  # Start loop
    time.sleep(1)  # Wait for connection setup to complete
    mqttc.subscribe(topic)


xml_files = [
    r"C:\Users\acs-opal\git\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_EQ.xml",
    r"C:\Users\acs-opal\git\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_SV.xml",
    r"C:\Users\acs-opal\git\cim-grid-data\CIGRE_MV\CIGRE_MV_no_tapchanger_With_LoadFlow_Results\Rootnet_FULL_NE_06J16h_TP.xml"]

res = cimpy.cimread(xml_files)
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

input("Press enter to stop client...")

mqttc.loop_stop()
mqttc.disconnect()
