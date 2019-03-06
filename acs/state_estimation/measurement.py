from enum import Enum
import numpy as np

class ElemType(Enum):
	Node = 1							#Node Voltage
	Branch = 2							#Complex Power flow at branch

class MeasType(Enum):
	V_mag = 1							#Node Voltage
	Sinj_real= 2						#Complex Power Injection at node
	Sinj_imag= 3						#Complex Power Injection at node
	S1_real = 4							#Active Power flow at branch, measured at initial node (S1.real)
	S1_imag = 5							#Reactive Power flow at branch, measured at initial node (S1.imag)
	I_mag = 6							#Branch Current
	Vpmu_mag = 7						#Node Voltage
	Vpmu_phase = 8						#Node Voltage
	Ipmu_mag = 9						#Branch Current
	Ipmu_phase = 10						#Branch Current
	S2_real = 11						#Active Power flow at branch, measured at final node (S2.real)
	S2_imag = 12						#Reactive Power flow at branch, measured at final node (S2.imag)

class Measurement():
	def __init__(self, element, element_type, meas_type, meas_value, std_dev):
		"""
		Creates a measurement, which is used by the estimation module. Possible types of measurements are: v, p, q, i, Vpmu and Ipmu
		@element: pointer to measured element
		@element_type: Clarifies which element is measured.
		@meas_type: 
		@meas_value: measurement value.
		@std_dev: standard deviation in the same unit as the measurement.
		"""
		
		if not isinstance(element_type, ElemType):
			raise Exception("elem_type must be an object of class ElemType")
		
		if not isinstance(meas_type, MeasType):
			raise Exception("meas_type must be an object of class MeasType")
			
		self.element = element
		self.element_type = element_type
		self.meas_type = meas_type
		self.meas_value = meas_value
		self.std_dev = std_dev
		self.std_dev = std_dev
		self.mval = 0.0					#measured values (affected by uncertainty)

class Measurents_set():
	def __init__(self):
		self.measurements = []			#array with all measurements
		
	def create_measurement(self, element, element_type, meas_type, meas_value, std_dev):
		"""
		to add elements to the measurements array
		"""
		self.measurements.append(Measurement(element, element_type, meas_type, meas_value, std_dev))
	
	def meas_creation(self):
		""" 
		It calculates the measured values (affected by uncertainty) at the measurement points
		"""
		err_pu = np.random.normal(0,1,len(self.measurements))
		for index, measurement in enumerate(self.measurements):
			measurement.mval = measurement.meas_value + self.std_dev*err_pu[index]

	def meas_creation_test(self, err_pu):
		""" 
		For test purposes.
		It calculates the measured values (affected by uncertainty) at the measurement points.
		This function takes as paramenter the random gaussian distribution. 
		"""
		for index, measurement in enumerate(self.measurements):
			measurement.mval = measurement.meas_value + measurement.std_dev*err_pu[index]
	
	def getMeasurements(self, type):
		"""
		return an array with all measurements of type "type" in the array Measurents_set.measurements.
		"""
		measurements = []
		for measurement in self.measurements:
			if measurement.meas_type is type:
				measurements.append(measurement)
		
		return measurements
	
	def getNumberOfMeasurements(self, type):
		"""
		return number of measurements of type "type" in the array Measurents_set.measurements
		"""
		number = 0
		for measurement in self.measurements:
			if measurement.meas_type is type:
				number = number+1
			
		return number
		
	def getIndexOfMeasurements(self, type):
		"""
		return index of all measurements of type "type" in the array Measurents_set.measurements
		"""
		idx = np.zeros(self.getNumberOfMeasurements(type), dtype=int)
		i = 0
		for index, measurement in enumerate(self.measurements):
			if measurement.meas_type is type:
				idx[i] = index
				i = i+1
				
		return idx
		
	def getWeightsMatrix(self):
		"""
		return an array the weights (obtained as standard_deviations^-2)
		"""
		weights = np.zeros(len(self.measurements))
		for index, measurement in enumerate(self.measurements):
			#the weight is small and can bring instability during matrix inversion, so we "cut" everything below 10^-6
			if measurement.std_dev<10**(-6):
				measurement.std_dev = 10**(-6)
			weights[index] = measurement.std_dev**(-2)
		
		return weights
	
	def getMeasValues(self):
		"""
		for test purposes
		returns an array with all measured values 
		"""
		meas_val = np.zeros(len(self.measurements))
		for index, measurement in enumerate(self.measurements):
			meas_val[index] = measurement.meas_value
		
		return meas_val
		
	def getmVal(self):
		"""
		returns an array with all measured values (affected by uncertainty)
		"""
		mVal = np.zeros(len(self.measurements))
		for index, measurement in enumerate(self.measurements):
			mVal[index] = measurement.mval			
		
		""" Replace in mVal amplitude and phase of Vpmu by real and imaginary part """
		#get all measurements of type MeasType.Vpmu_mag
		Vpmu_mag_idx = self.getIndexOfMeasurements(type=MeasType.Vpmu_mag)
		#get all measurements of type MeasType.Vpmu_phase
		Vpmu_phase_idx = self.getIndexOfMeasurements(type=MeasType.Vpmu_phase)
		for vpmu_mag_index, vpmu_phase_index in zip(Vpmu_mag_idx, Vpmu_phase_idx):
			vamp = self.measurements[vpmu_mag_index].mval
			vtheta = self.measurements[vpmu_phase_index].mval
			mVal[vpmu_mag_index] = vamp*np.cos(vtheta)
			mVal[vpmu_phase_index] = vamp*np.sin(vtheta)
		
		""" Replace in z amplitude and phase of Ipmu by real and imaginary part """
		#get all measurements of type MeasType.Ipmu_mag
		Ipmu_mag_idx = self.getIndexOfMeasurements(type=MeasType.Ipmu_mag)
		#get all measurements of type MeasType.Ipmu_phase
		Ipmu_phase_idx = self.getIndexOfMeasurements(type=MeasType.Ipmu_phase)
		for ipmu_mag_index, ipmu_phase_index in zip(Ipmu_mag_idx, Ipmu_phase_idx):
			iamp = self.measurements[ipmu_mag_index].mval
			itheta = self.measurements[ipmu_phase_index].mval
			mVal[ipmu_mag_index] = iamp*np.cos(itheta)
			mVal[ipmu_phase_index] = iamp*np.sin(itheta)
		
		return mVal
	
	def getStd_Dev(self):
		"""
		for test purposes
		returns an array with all standard deviations
		"""
		std_dev = np.zeros(len(self.measurements))
		for index, measurement in enumerate(self.measurements):
			std_dev[index] = measurement.std_dev			
		
		return std_dev