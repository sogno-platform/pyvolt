#!/usr/bin/env python

################
# Dependancies:
# paho, mqtt library for python (sudo pip install paho-mqtt)

import paho.mqtt.client as mqtt
import time
import json

def on_log(mqttc, userdata, level, buf):
    print("log: "+buf)

def on_connect(mqttc, userdata, flags, rc):
    if rc==0:
        mqttc.connected_flag=False #set flag
        print("Publisher connection OK")
    else:
        print("Bad connection, error code= ",rc)
    
#def on_message(mqttc, userdata, msg):
#    message = json.loads(msg.payload)

def connect_subscribe(topic):
    mqttc = mqtt.Client("SognoSE", True)
    mqttc.username_pw_set("villas","s3c0sim4!")
    mqttc.on_log = on_log
    mqttc.on_connect = on_connect
    #mqttc.on_message = on_message
    
    mqttc.connect("137.226.248.91")
    
    mqttc.loop_start()  #Start loop
    mqttc.subscribe(topic)
    
    return mqttc	