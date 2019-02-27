from enum import Enum
import numpy as np

class ElemType(Enum):
	Node = 1		#Node Voltage
	Branch = 2		#Complex Power Injection at node
	
class MeasType(Enum):
	V = 1		#Node Voltage
	Sinj= 2		#Complex Power Injection at node
	S1 = 4		#Complex Power flow at branch, measured at initial node
	S2 = 5		#Complex Power flow at branch, measured at final node
	I = 6		#Branch Current
	Vpmu = 7	#Node Voltage
	Ipmu = 9	#Branch Current

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
		self.std_dev = self.meas_value*(std_dev/300)
		if self.std_dev<10**(-6)
			self.std_dev = 10**(-6)
		self.mval = 0.0		#measured values (affected by uncertainty)
		
class Measurents_set():
	def __init__(self):
		self.measurements = []
		
	def create_measurement(self, element, element_type, meas_type, meas_value, std_dev):
		self.measurements.append(Measurement(element, element_type, meas_type, meas_value, std_dev))
	
	def meas_creation(self):
		""" 
		It calculates the measured values (affected by uncertainty) at the measurement points
		"""
		err_pu = np.random.normal(0,1,len(self.measurements))
		for index, measurement in enumerate(self.measurements):
			measurement.mval = measurement.meas_value + self.std_dev*err_pu[index]
			
	def getNumberOfMeasurements(self)
		"""
		return number of measurements of each type in the array Measurents_set.measurements
		"""
		nvi, npi, nqi, npf, nqf, nii, nvpum, nipmu = 0
		for elem in self.measurements:
			if elem.meas_type is MeasType.V:
				nvi = nvi+1
			elif elem.meas_type is MeasType.Sinj:
				npi = npi+1
				nqi = nqi+1
			elif elem.meas_type is MeasType.S1 or elem.meas_type is MeasType.S2:
				npf = npf+1
				nqf = nqf+1
			elif elem.meas_type is MeasType.I:
				nii = nii+1
			elif elem.meas_type is MeasType.Vpmu:
				nvpum = nvpum+1
			elif elem.meas_type is MeasType.Ipmu:
				nipmu = nipmu+1
			
		return nvi, npi, nqi, npf, nqf, nii, nvpum, nipmu
		
	def getMeasuredActiveInjPowers(self):
		"""
		return an array with the measurements of type Sinj.real
		"""
		Pinj = np.array([])
		for elem in self.measurements:
			if elem.meas_type is MeasType.Sinj:
				Pinj = np.append(Pinj, elem.real)
		
		return Pinj
		
	def getMeasuredReactiveInjPowers(self):
		"""
		return an array with the measurements of type Sinj.imag
		"""
		Qinj = np.array([])
		for elem in self.measurements:
			if elem.meas_type is MeasType.Sinj:
				Qinj = np.append(Qinj, elem.imag)
		
		return Qinj
		
	def getMeasuredActiveBPowers(self):
		"""
		return an array with the measurements of type S1.real or S2.real
		"""
		Pbr = np.array([])
		for elem in self.measurements:
			if elem.meas_type is MeasType.S1 or elem.meas_type is MeasType.S2:
				Pbr = np.append(Pbr, elem.real)
		return Pbr
		
	def getMeasuredReactiveBPowers(self):
		"""
		return an array with the measurements of type S1.imag or S2.imag
		"""
		Qbr = np.array([])
		for elem in self.measurements:
			if elem.meas_type is MeasType.S1 or elem.meas_type is MeasType.S2:
				Qbr = np.append(Qbr, elem.imag)
		return Qbr
		
	def getWeightsMatrix(self)
		"""
		creates the weights matrix (obtained as standard_deviations^-2)
		"""
		weights = np.zeros(len(self.measurements))
		for index, measurement in enumerate(self.measurements):
			weights[index] = measurement.std_dev**(-2)
		
		return np.diag(weights)
	