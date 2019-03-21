import json
import paho.mqtt.client as mqtt
import time
import numpy as np

# The callback for when the client receives a CONNACK response from the server.
def on_connect(client, userdata, flags, rc):
    if rc==0:
        client.connected_flag=True      #set flag
        print("connected OK with returned code=",rc)
    else:
        print("Bad connection with returned code=",rc)

# The callback for when a PUBLISH message is received from the server.
def on_message(client, userdata, msg):
    global connected_clients, sequence, V_est_sogno_cigre_se, V_est_sogno_cigre_se_cim, end
    if connected_clients<2:
        #check that the two clients are connected
        if json.loads(msg.payload) == "OK":
            connected_clients +=1
            print("client" + str(connected_clients) + " connected")
            if connected_clients==2:
                mqttc.publish(topic_clients, data[sequence])
    else:
        #receive and compare the results of the state estimation
        message = json.loads(msg.payload)[0]
        from_client = message['client']
        V_est_mag = np.array(message['V_est_mag'])
        V_est_phase = np.array(message['V_est_phase'])
        if from_client == "Sogno_cigre_se_cim":
            V_est_sogno_cigre_se_cim = V_est_mag + 1j*V_est_phase
        elif from_client == "Sogno_cigre_se":
            V_est_sogno_cigre_se = V_est_mag + 1j*V_est_phase
        if V_est_sogno_cigre_se_cim is not None and V_est_sogno_cigre_se is not None:
            if (np.array_equal(np.around(V_est_sogno_cigre_se_cim,3), np.around(V_est_sogno_cigre_se,3))):
                print("V_est_sogno_cigre_se==V_est_sogno_cigre_se_cim (sequence = " + str(sequence) +")?: True ")
            else:
                print("V_est_sogno_cigre_se==V_est_sogno_cigre_se_cim (sequence = " + str(sequence) +")?: False ")
                print(V_est_sogno_cigre_se_cim - V_est_sogno_cigre_se)
            sequence += 1
            if sequence == len_data+1:
                end = True
                mqttc.publish(topic_clients, json.dumps("end"))
            else:
                mqttc.publish(topic_clients, data[sequence])
                V_est_sogno_cigre_se = None
                V_est_sogno_cigre_se_cim = None

#global variables used in callback on_message
end = False
connected_clients = 0 
sequence = 1
V_est_sogno_cigre_se = None
V_est_sogno_cigre_se_cim = None
len_data = 300          #number of sequences to use

#topic names
topic = "dpsim-powerflow_results"
topic_clients = "dpsim-powerflow"

broker_adress = "m16.cloudmqtt.com"
mqttc = mqtt.Client("SognoDemo", True)          #create new instance
mqttc.username_pw_set("ilgtdaqk", "UbNQQjmcUdqq")
mqttc.on_connect = on_connect                   #attach function to callback
mqttc.on_message = on_message                   #attach function to callback
mqttc.connect(broker_adress, 14543)             #connect to broker
mqttc.loop_start()                              #start loop to process callback
time.sleep(1)                                   #wait for connection setup to complete
mqttc.subscribe(topic)

data_file = r"C:\Users\Martin\Desktop\hiwi\git\state-estimation-client\sogno_demo\dpsim_powerflow_record_cigre.txt"
data = []
with open(data_file) as json_file:
    for line in json_file:
        data.append(line)

while not end:
    time.sleep(0.01)

mqttc.loop_stop()       #Stop loop 
mqttc.disconnect()      #disconnect