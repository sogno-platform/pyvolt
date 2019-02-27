import sys
import numpy
import nv_powerflow_cim
from measurement import *

sys.path.append("../../../dataprocessing")
from villas.dataprocessing.readtools import read_timeseries_dpsim


def DsseCall(system, zdata, measurements):
	""" It identifies the type of measurements present in the measurement set and 
	calls the appropriate estimator for dealing with them."""
	# system: model of the system (nodes, lines, topology)
	# zdata: Vector of measurements in Input (voltages, currents, powers)

	# 1. Select type of Estimator.
	# If at least a PMU is present we launch the combined estimator, otherwise a simple traditional etimator 
	Vmag_meas = 0
	Vpmu_meas = 0
	for elem in measurements.measurements_set:
		if elem.meas_type=="v":
			Vmag_meas = Vmag_meas + 1
		elif elem.meas_type=="Vpmu":
			Vpmu_meas = Vpmu_meas + 1
	
	trad_code = 1 if Vmag_meas >0 else 0
	PMU_code = 2 if Vpmu_meas >0 else 0
	est_code = trad_code + PMU_code
		
	# 2. Run Estimator.
	if est_code == 1:
		Vest = DsseTrad(system, zdata)
	elif est_code == 2:
		Vest = DssePmu(system, zdata)
	else:
		Vest = DsseMixed(system, zdata)
		
	# 3. Populate arrays with data in output of estimator
	# The Estimator provides:
	# Vest: Node Voltage phasor estimated 
	# Iest: Branch Current phasor estimated   
	# Iinjest: Injection Current phasor estimated 
	# S1est:  Complex Power flow at branch, measured at initial node
	# S2est:  Complex Power flow at branch, measured at final node  
	# Sinjest: Complex Power Injection at node

	""" From here on, all the other quantities of the grid are calculated """
	results = nv_powerflow_cim.PowerflowResults(system)
	results.load_voltages(Vest)
	results.calculateI()
	results.calculateIinj()
	results.calculateSinj()
	results.calculateI()
	results.calculateS()
	
	return results

def DsseTrad(system, zdata):
	""" It performs state estimation using rectangular node voltage state variables 
	and it is customized to work without PMU measurements"""
	# Traditional state estimator
	# system: model of the system (nodes, lines, topology)
	# zdata: Vector of measurements in Input (voltages, currents, powers)	
	
	# number of nodes of the grids, identify also the number of states (2*nodes_num-1)
	nodes_num = len(system.nodes) 
			
	vidx = numpy.where(zdata.mtype==1) #voltage input measurement
	pidx = numpy.where(zdata.mtype==2) #active power input measurement
	qidx = numpy.where(zdata.mtype==3) #reactive power input measurement
	pfidx = numpy.where(zdata.mtype==4) #active power flow input measurement
	qfidx = numpy.where(zdata.mtype==5) #reactive power flow input measurement
	iidx = numpy.where(zdata.mtype==6) #current flow input measurement
	
	nvi = len(vidx[0]) # number v.measure
	npi = len(pidx[0]) # number power.act.measure
	nqi = len(qidx[0]) # number power.react.measure
	npf = len(pfidx[0]) # number power.act.flow.measure
	nqf = len(qfidx[0]) # number power.react.flow.measure
	nii = len(iidx[0]) #number c.flow.measure
	
	busvi = zdata.mfrom[vidx] #node where v.measure is taken
	buspi = zdata.mfrom[pidx] #node where power.act.measure is taken
	busqi = zdata.mfrom[qidx] #node where power.react.measure is taken
	fbuspf = zdata.mfrom[pfidx] #node from where power.act.measure starts
	tbuspf = zdata.mto[pfidx] #node to where power.act.measure goes
	fbusqf = zdata.mfrom[qfidx] #node from where power.react.measure starts
	#tbusqf = zdata.mto[qfidx] #should be added [aan] ? #node to where power.react.measure goes
	fbusiamp = zdata.mfrom[iidx] #node from where current.flow.measure starts
	tbusiamp = zdata.mto[iidx] #node to where current.flow.measure goes
	   
	z = zdata.mval #values of measurements are stored in array z
	
	#extrapolate different types of measurements
	Pinj = z[pidx] 
	Qinj = z[qidx]
	Pbr = z[pfidx]
	Qbr = z[qfidx]
	
	# zero injection measurements
	# the weight is small and can bring instability during matrix inversion, so we "cut" everything below 10^-6
	idx = numpy.where(zdata.mstddev<10**(-6)) 
	zdata.mstddev[idx] = 10**(-6)
	
	# weights matrix is obtained as stdandard_deviations^-2
	weights = zdata.mstddev**(-2)
	W = numpy.diag(weights)
	
	# Admittances of the lines of the network
	Gmatrix = system.Ymatrix.real
	Bmatrix = system.Ymatrix.imag
	Yabs_matrix = numpy.absolute(system.Ymatrix)
	Yphase_matrix = numpy.angle(system.Ymatrix)
	
	# Jacobian Matrix. Includes the derivatives of the measurements (voltages, currents, powers) with respect to the states (voltages)
	
	""" Jacobian for Power Injection Measurements (converted to equivalent 
	rectangualar current measurements) """
	# Derivative of power injection(converted to current injection) with respect to voltage: it is the admittance. 
	H2 = numpy.zeros((npi,2*nodes_num-1))
	H3 = numpy.zeros((nqi,2*nodes_num-1))
	for i in range(npi):
		m = buspi[i]-1
		m2 = m + nodes_num - 1
		H2[i][m] = - Gmatrix[m][m]
		H2[i][m2] = Bmatrix[m][m]
		H3[i][m] = - Bmatrix[m][m]
		H3[i][m2] = - Gmatrix[m][m]
		idx = numpy.subtract(system.Adj[m],1)
		H2[i][idx] = - Gmatrix[m][idx]
		H3[i][idx] = - Bmatrix[m][idx]
		if 0 in idx:
			pos = numpy.where(idx==0)
			idx = numpy.delete(idx,pos)
		idx2 = idx + nodes_num-1
		H2[i][idx2] = Bmatrix[m][idx]
		H3[i][idx2] = - Gmatrix[m][idx]
		
	""" Jacobian for branch Power Measurements (converted to equivalent 
	rectangualar current measurements)"""
	# Derivative of branch power flow(converted to current flow) with respect to voltage: it is the admittance. 
	H4 = numpy.zeros((npf,2*nodes_num-1))
	H5 = numpy.zeros((nqf,2*nodes_num-1))
	for i in range(npf):
		m = fbuspf[i]-1
		n = tbuspf[i]-1
		H4[i][m] = - Gmatrix[m][n]
		H4[i][n] = Gmatrix[m][n]
		H5[i][m] = - Bmatrix[m][n]
		H5[i][n] = Bmatrix[m][n]
		if m > 0:
			m2 = m + nodes_num-1
			H4[i][m2] = Bmatrix[m][n]
			H5[i][m2] = - Gmatrix[m][n]
		if n > 0:
			n2 = n + nodes_num-1
			H4[i][n2] = - Bmatrix[m][n]
			H5[i][n2] = Gmatrix[m][n]
		
	epsilon = 5 # treshold to stop Netwon Rapson iterations
	Vr = numpy.ones(nodes_num) #initialize voltage real part to 1 per unit
	Vx = numpy.zeros(nodes_num) #initialize voltage imaginary part to 0 per unit
	V = Vr + 1j*Vx
	num_iter = 0 # number of iterations of Newton Rapson
	
	StateVr = numpy.ones(nodes_num) #initialize voltage real part to 1 per unit
	StateVx = numpy.zeros(nodes_num-1) #initialize voltage imaginary part to 0 per unit
	State = numpy.concatenate((StateVr,StateVx),axis=0)
	
	# Iteration of Netwon Rapson method: needed to solve non-linear system of equation
	while epsilon>10**(-6):
		""" Computation of equivalent current measurements in place of the power measurements """
		# in every iteration the input power measurements are converted into currents by dividing by the voltage estimated at the previous iteration
		Irinj = (Pinj*V[buspi-1].real + Qinj*V[busqi-1].imag)/(numpy.absolute(V[buspi-1])**2)
		Ixinj = (Pinj*V[buspi-1].imag - Qinj*V[busqi-1].real)/(numpy.absolute(V[buspi-1])**2)
		z[pidx] = Irinj
		z[qidx] = Ixinj
		
		Irbr = (Pbr*V[fbuspf-1].real + Qbr*V[fbusqf-1].imag)/(numpy.absolute(V[fbuspf-1])**2)
		Ixbr = (Pbr*V[fbuspf-1].imag - Qbr*V[fbusqf-1].real)/(numpy.absolute(V[fbuspf-1])**2)
		z[pfidx] = Irbr
		z[qfidx] = Ixbr
		
		""" Voltage Magnitude Measurements """
		# at every iteration we update h(x) vector where V measure are available
		h1 = numpy.absolute(V[busvi-1])
		# the Jacobian rows where voltage measurements are presents is updated
		H1 = numpy.zeros((nvi,2*nodes_num-1))
		for i in range(nvi):
			m = busvi[i]-1
			H1[i][m] = numpy.cos(numpy.angle(V[m]))
			if m > 0:
				m2 = m + nodes_num-1
				H1[i][m2] = numpy.sin(numpy.angle(V[m]))
				
		""" Power Injection Measurements """
		# h(x) vector where power injections are present
		h2 = numpy.inner(H2,State)
		h3 = numpy.inner(H3,State)
		
		""" Power Flow Measurements """
		# h(x) vector where power flows are present
		h4 = numpy.inner(H4,State)
		h5 = numpy.inner(H5,State)
		
		""" Current Magnitude Measurements """
		# h(x) vector where current flows are present
		h6re = numpy.zeros((nii))
		h6im = numpy.zeros((nii))
		h6complex = numpy.zeros((nii),dtype=complex)
		h6 = numpy.ones((nii))
		H6 = numpy.zeros((nii,2*nodes_num-1))
		for i in range(nii):
			m = fbusiamp[i]-1
			n = tbusiamp[i]-1
			h6re[i] = Yabs_matrix[m][n]*((V[n].real-V[m].real)*numpy.cos(Yphase_matrix[m][n]) + (V[m].imag-V[n].imag)*numpy.sin(Yphase_matrix[m][n]))
			h6im[i] = Yabs_matrix[m][n]*((V[n].real-V[m].real)*numpy.sin(Yphase_matrix[m][n]) + (V[n].imag-V[m].imag)*numpy.cos(Yphase_matrix[m][n]))
			h6complex[i] = h6re[i] + 1j*h6im[i]
			if num_iter>0:
				h6[i] = numpy.absolute(h6complex[i])
			H6[i][m] = - Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6re[i] + numpy.sin(Yphase_matrix[m][n])*h6im[i])/h6[i]
			H6[i][n] = Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6re[i] + numpy.sin(Yphase_matrix[m][n])*h6im[i])/h6[i]
			if m > 0:
				m2 = m + nodes_num-1
				H6[i][m2] = - Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6im[i] - numpy.sin(Yphase_matrix[m][n])*h6re[i])/h6[i]
			if n > 0:
				n2 = n + nodes_num-1
				H6[i][n2] = Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6im[i] - numpy.sin(Yphase_matrix[m][n])*h6re[i])/h6[i]
		
		""" WLS computation """
		# all the sub-matrixes of H calcualted so far are merged in a unique matrix
		H = numpy.concatenate((H1,H2,H3,H4,H5,H6),axis=0)
		# h(x) sub-vectors are concatenated
		y = numpy.concatenate((h1,h2,h3,h4,h5,h6),axis=0)
		# "res" is the residual vector. The difference between input measurements and h(x)
		res = numpy.subtract(z,y)
		# g = transpose(H) * W * res
		g = numpy.inner(H.transpose(),numpy.inner(W,res))
		WH = numpy.inner(W,H.transpose())
		# G is the gain matrix, that will have to be inverted at each iteration
		G = numpy.inner(H.transpose(),WH.transpose())
		# inversion of G
		Ginv = numpy.linalg.inv(G)
		# Delta includes the updates of the states for the current Newton Rapson iteration
		Delta_State = numpy.inner(Ginv,g)
		# state is updated
		State = State + Delta_State
		# calculate the NR treeshold (for the next while check)
		epsilon = numpy.amax(numpy.absolute(Delta_State))
		# update the voltages 
		V.real = State[:nodes_num]
		V.imag = State[nodes_num:]
		V.imag = numpy.concatenate(([0],V.imag), axis=0)
		#V = Real_to_all(V.real, V.imag)
		
		num_iter = num_iter+1
	
	return V

def DssePmu(system, zdata):
	""" It performs state estimation using rectangular node voltage state variables 
	and it is customized to work using PMU measurements."""
	
	nodes_num = len(system.nodes)
	
	pidx = numpy.where(zdata.mtype==2)
	qidx = numpy.where(zdata.mtype==3)
	pfidx = numpy.where(zdata.mtype==4)
	qfidx = numpy.where(zdata.mtype==5)
	vmagpmuidx = numpy.where(zdata.mtype==7)
	vphasepmuidx = numpy.where(zdata.mtype==8)
	imagpmuidx = numpy.where(zdata.mtype==9)
	iphasepmuidx = numpy.where(zdata.mtype==10)
	
	npi = len(pidx[0])
	nqi = len(qidx[0])
	npf = len(pfidx[0])
	nqf = len(qfidx[0])
	nvpmu = len(vmagpmuidx[0])
	nipmu = len(imagpmuidx[0])
	
	buspi = zdata.mfrom[pidx]
	busqi = zdata.mfrom[qidx]
	fbuspf = zdata.mfrom[pfidx]
	tbuspf = zdata.mto[pfidx]
	fbusqf = zdata.mfrom[qfidx]
	busvpmu = zdata.mfrom[vmagpmuidx]
	fbusipmu = zdata.mfrom[imagpmuidx]
	tbusipmu = zdata.mto[imagpmuidx]
	   
	z = zdata.mval
	Pinj = z[pidx]
	Qinj = z[qidx]
	Pbr = z[pfidx]
	Qbr = z[qfidx]

	idx = numpy.where(zdata.mstddev<10**(-6))
	zdata.mstddev[idx] = 10**(-6)
	weights = zdata.mstddev**(-2)
	W = numpy.diag(weights)
	
	Gmatrix = system.Ymatrix.real
	Bmatrix = system.Ymatrix.imag
	
	""" Jacobian for Power Injection Measurements (converted to equivalent 
	rectangualar current measurements) """
	H2 = numpy.zeros((npi,2*nodes_num))
	H3 = numpy.zeros((nqi,2*nodes_num))
	for i in range(npi):
		m = buspi[i] - 1
		m2 = m + nodes_num
		H2[i][m] = - Gmatrix[m][m]
		H2[i][m2] = Bmatrix[m][m]
		H3[i][m] = - Bmatrix[m][m]
		H3[i][m2] = - Gmatrix[m][m]
		idx = numpy.subtract(system.Adj[m],1)
		H2[i][idx] = - Gmatrix[m][idx]
		H3[i][idx] = - Bmatrix[m][idx]
		idx2 = idx + nodes_num
		H2[i][idx2] = Bmatrix[m][idx]
		H3[i][idx2] = - Gmatrix[m][idx]
		
	""" Jacobian for branch Power Measurements (converted to equivalent 
	rectangualar current measurements)"""
	H4 = numpy.zeros((npf,2*nodes_num))
	H5 = numpy.zeros((nqf,2*nodes_num))
	for i in range(npf):
		m = fbuspf[i]-1
		n = tbuspf[i]-1
		H4[i][m] = - Gmatrix[m][n]
		H4[i][n] = Gmatrix[m][n]
		H5[i][m] = - Bmatrix[m][n]
		H5[i][n] = Bmatrix[m][n]
		m2 = m + nodes_num
		H4[i][m2] = Bmatrix[m][n]
		H5[i][m2] = - Gmatrix[m][n]
		n2 = n + nodes_num
		H4[i][n2] = - Bmatrix[m][n]
		H5[i][n2] = Gmatrix[m][n]
		
	""" Jacobian for Voltage Pmu Measurements (converted into rectangular) """
	H7 = numpy.zeros((nvpmu,2*nodes_num))
	H8 = numpy.zeros((nvpmu,2*nodes_num))
	for i in range(nvpmu):
		idx1 = vmagpmuidx[0][i]
		idx2 = vphasepmuidx[0][i]
		vamp = z[idx1]
		vtheta = z[idx2]
		z[idx1] = vamp*numpy.cos(vtheta)
		z[idx2] = vamp*numpy.sin(vtheta)
		rot_mat = numpy.array([[numpy.cos(vtheta), - vamp*numpy.sin(vtheta)], [numpy.sin(vtheta), vamp*numpy.cos(vtheta)]])
		starting_cov = numpy.array([[weights[idx1], 0], [0, weights[idx2]]])
		final_cov = numpy.inner(rot_mat,numpy.inner(starting_cov,rot_mat.transpose()))
		W[idx1][idx1] = final_cov[0][0]
		W[idx2][idx2] = final_cov[1][1]
		W[idx1][idx2] = final_cov[0][1]
		W[idx2][idx1] = final_cov[1][0]
		m = busvpmu[i]-1
		H7[i][m] = 1
		m2 = m + nodes_num
		H8[i][m2] = 1
		
	""" Jacobian for Current Pmu Measurements (converted into rectangular) """
	H9 = numpy.zeros((nipmu,2*nodes_num))
	H10 = numpy.zeros((nipmu,2*nodes_num))
	for i in range(nipmu):
		idx1 = imagpmuidx[0][i]
		idx2 = iphasepmuidx[0][i]
		iamp = z[idx1]
		itheta = z[idx2]
		z[idx1] = iamp*numpy.cos(itheta)
		z[idx2] = iamp*numpy.sin(itheta)
		rot_mat = numpy.array([[numpy.cos(itheta), - iamp*numpy.sin(itheta)], [numpy.sin(itheta), iamp*numpy.cos(itheta)]])
		starting_cov = numpy.array([[weights[idx1], 0], [0, weights[idx2]]])
		final_cov = numpy.inner(rot_mat,numpy.inner(starting_cov,rot_mat.transpose()))
		W[idx1][idx1] = final_cov[0][0]
		W[idx2][idx2] = final_cov[1][1]
		W[idx1][idx2] = final_cov[0][1]
		W[idx2][idx1] = final_cov[1][0]
		m = fbusipmu[i]-1
		n = tbusipmu[i]-1
		H9[i][m] = - Gmatrix[m][n]
		H9[i][n] = Gmatrix[m][n]
		H10[i][m] = - Bmatrix[m][n]
		H10[i][n] = Bmatrix[m][n]
		m2 = m + nodes_num
		n2 = n + nodes_num
		H9[i][m2] = Bmatrix[m][n]
		H9[i][n2] = - Bmatrix[m][n]
		H10[i][m2] = - Gmatrix[m][n]
		H10[i][n2] = Gmatrix[m][n]
		
	epsilon = 5
	Vr = numpy.ones(nodes_num)
	Vx = numpy.zeros(nodes_num)
	V = Vr + 1j*Vx
	num_iter = 0
	
	StateVr = numpy.ones(nodes_num)
	StateVx = numpy.zeros(nodes_num)
	State = numpy.concatenate((StateVr,StateVx),axis=0)
	
	H = numpy.concatenate((H2,H3,H4,H5,H7,H8,H9,H10),axis=0)
	WH = numpy.inner(W,H.transpose())
	G = numpy.inner(H.transpose(),WH.transpose())
	Ginv = numpy.linalg.inv(G)
	
	while epsilon>10**(-6):
		""" Computation of equivalent current measurements in place of the power measurements """
		Irinj = (Pinj*V[buspi-1].real + Qinj*V[busqi-1].imag)/(numpy.absolute(V[buspi-1])**2)
		Ixinj = (Pinj*V[buspi-1].imag - Qinj*V[busqi-1].real)/(numpy.absolute(V[buspi-1])**2)
		z[pidx] = Irinj
		z[qidx] = Ixinj
		
		Irbr = (Pbr*V[fbuspf-1].real + Qbr*V[fbusqf-1].imag)/(numpy.absolute(V[fbuspf-1])**2)
		Ixbr = (Pbr*V[fbuspf-1].imag - Qbr*V[fbusqf-1].real)/(numpy.absolute(V[fbuspf-1])**2)
		z[pfidx] = Irbr
		z[qfidx] = Ixbr
		
		
		""" WLS computation """
		y = numpy.inner(H,State)
		res = numpy.subtract(z,y)
		g = numpy.inner(H.transpose(),numpy.inner(W,res))

		Delta_State = numpy.inner(Ginv,g)
		
		State = State + Delta_State
		epsilon = numpy.amax(numpy.absolute(Delta_State))
		
		V.real = State[:nodes_num]
		V.imag = State[nodes_num:]
		
		num_iter = num_iter+1
		
	return V
	
def DsseMixed(system, measurements):
	""" It performs state estimation using rectangular node voltage state variables
	and it is built to work in scenarios where both conventional and PMU measurements 
	are simultaneously present."""
	
	vidx = numpy.where(zdata.mtype==1)
	pidx = numpy.where(zdata.mtype==2)
	qidx = numpy.where(zdata.mtype==3)
	pfidx = numpy.where(zdata.mtype==4)
	qfidx = numpy.where(zdata.mtype==5)
	iidx = numpy.where(zdata.mtype==6)
	vmagpmuidx = numpy.where(zdata.mtype==7)
	vphasepmuidx = numpy.where(zdata.mtype==8)
	imagpmuidx = numpy.where(zdata.mtype==9)
	iphasepmuidx = numpy.where(zdata.mtype==10)
	
	#nvi = len(vidx[0])
	#npi = len(pidx[0])
	#nqi = len(qidx[0])
	#npf = len(pfidx[0])
	#nqf = len(qfidx[0])
	#nii = len(iidx[0])
	#nvpmu = len(vmagpmuidx[0])
	#nipmu = len(imagpmuidx[0])
	
	#calculate number of measurements of each type
	nvi, npi, nqi, npf, nqf, nii, nvpum, nipmu = measurements.getNumberOfMeasurements()
	
	busvi = zdata.mfrom[vidx]
	buspi = zdata.mfrom[pidx]
	busqi = zdata.mfrom[qidx]
	fbuspf = zdata.mfrom[pfidx]
	tbuspf = zdata.mto[pfidx]
	fbusqf = zdata.mfrom[qfidx]
	fbusiamp = zdata.mfrom[iidx]
	tbusiamp = zdata.mto[iidx]
	busvpmu = zdata.mfrom[vmagpmuidx]
	fbusipmu = zdata.mfrom[imagpmuidx]
	tbusipmu = zdata.mto[imagpmuidx]
	   
	#z = zdata.mval
	#Pinj = z[pidx]
	#Qinj = z[qidx]
	#Pbr = z[pfidx]
	#Qbr = z[qfidx]
	
	Pinj = measurements.getMeasuredActiveInjPowers()
	Qinj = measurements.getMeasuredReactiveInjPowers()
	Pbr = measurements.getMeasuredActiveBPowers()
	Qbr = measurements.getMeasuredReactiveBPowers()
	
	#idx = numpy.where(zdata.mstddev<10**(-6))
	#zdata.mstddev[idx] = 10**(-6)
	#weights = zdata.mstddev**(-2)
	#W = numpy.diag(weights)
	# weights matrix is obtained as stdandard_deviations^-2
	W = measurements.getWeightsMatrix()
	
	Gmatrix = system.Ymatrix.real
	Bmatrix = system.Ymatrix.imag
	Yabs_matrix = numpy.absolute(system.Ymatrix)
	Yphase_matrix = numpy.angle(system.Ymatrix)
		
	""" Jacobian for Power Injection Measurements (converted to equivalent 
	rectangualar current measurements) """
	H2 = numpy.zeros((npi,2*len(system.nodes)))
	H3 = numpy.zeros((nqi,2*len(system.nodes)))
	for i in range(npi):
		m = buspi[i] - 1
		m2 = m + len(system.nodes)
		H2[i][m] = - Gmatrix[m][m]
		H2[i][m2] = Bmatrix[m][m]
		H3[i][m] = - Bmatrix[m][m]
		H3[i][m2] = - Gmatrix[m][m]
		idx = numpy.subtract(system.Adj[m],1)
		H2[i][idx] = - Gmatrix[m][idx]
		H3[i][idx] = - Bmatrix[m][idx]
		idx2 = idx + len(system.nodes)
		H2[i][idx2] = Bmatrix[m][idx]
		H3[i][idx2] = - Gmatrix[m][idx]
	
	""" Jacobian for branch Power Measurements (converted to equivalent 
	rectangualar current measurements)"""
	H4 = numpy.zeros((npf,2*len(system.nodes)))
	H5 = numpy.zeros((nqf,2*len(system.nodes)))
	for i in range(npf):
		m = fbuspf[i]-1
		n = tbuspf[i]-1
		H4[i][m] = - Gmatrix[m][n]
		H4[i][n] = Gmatrix[m][n]
		H5[i][m] = - Bmatrix[m][n]
		H5[i][n] = Bmatrix[m][n]
		m2 = m + len(system.nodes)
		H4[i][m2] = Bmatrix[m][n]
		H5[i][m2] = - Gmatrix[m][n]
		n2 = n + len(system.nodes)
		H4[i][n2] = - Bmatrix[m][n]
		H5[i][n2] = Gmatrix[m][n]
		
	""" Jacobian for Voltage Pmu Measurements (converted into rectangular) """
	H7 = numpy.zeros((nvpmu,2*len(system.nodes)))
	H8 = numpy.zeros((nvpmu,2*len(system.nodes)))
	for i in range(nvpmu):
		idx1 = vmagpmuidx[0][i]
		idx2 = vphasepmuidx[0][i]
		vamp = z[idx1]
		vtheta = z[idx2]
		z[idx1] = vamp*numpy.cos(vtheta)
		z[idx2] = vamp*numpy.sin(vtheta)
		rot_mat = numpy.array([[numpy.cos(vtheta), - vamp*numpy.sin(vtheta)], [numpy.sin(vtheta), vamp*numpy.cos(vtheta)]])
		starting_cov = numpy.array([[weights[idx1], 0], [0, weights[idx2]]])
		final_cov = numpy.inner(rot_mat,numpy.inner(starting_cov,rot_mat.transpose()))
		W[idx1][idx1] = final_cov[0][0]
		W[idx2][idx2] = final_cov[1][1]
		W[idx1][idx2] = final_cov[0][1]
		W[idx2][idx1] = final_cov[1][0]
		m = busvpmu[i]-1
		H7[i][m] = 1
		m2 = m + len(system.nodes)
		H8[i][m2] = 1
	
	""" Jacobian for Current Pmu Measurements (converted into rectangular) """
	H9 = numpy.zeros((nipmu,2*len(system.nodes)))
	H10 = numpy.zeros((nipmu,2*len(system.nodes)))
	for i in range(nipmu):
		idx1 = imagpmuidx[0][i]
		idx2 = iphasepmuidx[0][i]
		iamp = z[idx1]
		itheta = z[idx2]
		z[idx1] = iamp*numpy.cos(itheta)
		z[idx2] = iamp*numpy.sin(itheta)
		rot_mat = numpy.array([[numpy.cos(itheta), - iamp*numpy.sin(itheta)], [numpy.sin(itheta), iamp*numpy.cos(itheta)]])
		starting_cov = numpy.array([[weights[idx1], 0], [0, weights[idx2]]])
		final_cov = numpy.inner(rot_mat,numpy.inner(starting_cov,rot_mat.transpose()))
		W[idx1][idx1] = final_cov[0][0]
		W[idx2][idx2] = final_cov[1][1]
		W[idx1][idx2] = final_cov[0][1]
		W[idx2][idx1] = final_cov[1][0]
		m = fbusipmu[i]-1
		n = tbusipmu[i]-1
		H9[i][m] = - Gmatrix[m][n]
		H9[i][n] = Gmatrix[m][n]
		H10[i][m] = - Bmatrix[m][n]
		H10[i][n] = Bmatrix[m][n]
		m2 = m + len(system.nodes)
		n2 = n + len(system.nodes)
		H9[i][m2] = Bmatrix[m][n]
		H9[i][n2] = - Bmatrix[m][n]
		H10[i][m2] = - Gmatrix[m][n]
		H10[i][n2] = Gmatrix[m][n]
		
	epsilon = 5
	Vr = numpy.ones(len(system.nodes))
	Vx = numpy.zeros(len(system.nodes))
	V = Vr + 1j*Vx
	#V = Real_to_all(Vr, Vx)
	num_iter = 0
	
	StateVr = numpy.ones(len(system.nodes))
	StateVx = numpy.zeros(len(system.nodes))
	State = numpy.concatenate((StateVr,StateVx),axis=0)
	
	while epsilon>10**(-6):
		""" Computation of equivalent current measurements in place of the power measurements """
		Irinj = (Pinj*V[buspi-1].real + Qinj*V[busqi-1].imag)/(numpy.absolute(V[buspi-1])**2)
		Ixinj = (Pinj*V[buspi-1].imag - Qinj*V[busqi-1].real)/(numpy.absolute(V[buspi-1])**2)
		z[pidx] = Irinj
		z[qidx] = Ixinj
		
		Irbr = (Pbr*V[fbuspf-1].real + Qbr*V[fbusqf-1].imag)/(numpy.absolute(V[fbuspf-1])**2)
		Ixbr = (Pbr*V[fbuspf-1].imag - Qbr*V[fbusqf-1].real)/(numpy.absolute(V[fbuspf-1])**2)
		z[pfidx] = Irbr
		z[qfidx] = Ixbr
		
		""" Voltage Magnitude Measurements """
		h1 = numpy.absolute(V[busvi-1])
		H1 = numpy.zeros((nvi,2*len(system.nodes)))
		for i in range(nvi):
			m = busvi[i]-1
			H1[i][m] = numpy.cos(numpy.angle(V[m]))
			m2 = m + len(system.nodes)
			H1[i][m2] = numpy.sin(numpy.angle(V[m]))
				
		""" Power Injection Measurements """
		h2 = numpy.inner(H2,State)
		h3 = numpy.inner(H3,State)
		
		""" Power Flow Measurements """
		h4 = numpy.inner(H4,State)
		h5 = numpy.inner(H5,State)
		
		""" Current Magnitude Measurements """
		h6re = numpy.zeros((nii))
		h6im = numpy.zeros((nii))
		h6complex = numpy.zeros((nii),dtype=complex)
		h6 = numpy.ones((nii))
		H6 = numpy.zeros((nii,2*len(system.nodes)))
		for i in range(nii):
			m = fbusiamp[i]-1
			n = tbusiamp[i]-1
			h6re[i] = Yabs_matrix[m][n]*((V[n].real-V[m].real)*numpy.cos(Yphase_matrix[m][n]) + (V[m].imag-V[n].imag)*numpy.sin(Yphase_matrix[m][n]))
			h6im[i] = Yabs_matrix[m][n]*((V[n].real-V[m].real)*numpy.sin(Yphase_matrix[m][n]) + (V[n].imag-V[m].imag)*numpy.cos(Yphase_matrix[m][n]))
			h6complex[i] = h6re[i] + 1j*h6im[i]
			if num_iter>0:
				h6[i] = numpy.absolute(h6complex[i])
			H6[i][m] = - Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6re[i] + numpy.sin(Yphase_matrix[m][n])*h6im[i])/h6[i]
			H6[i][n] = Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6re[i] + numpy.sin(Yphase_matrix[m][n])*h6im[i])/h6[i]
			m2 = m + len(system.nodes)
			H6[i][m2] = - Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6im[i] - numpy.sin(Yphase_matrix[m][n])*h6re[i])/h6[i]
			n2 = n + len(system.nodes)
			H6[i][n2] = Yabs_matrix[m][n]*(numpy.cos(Yphase_matrix[m][n])*h6im[i] - numpy.sin(Yphase_matrix[m][n])*h6re[i])/h6[i]
		
		""" PMU Voltage Measurements """
		h7 = numpy.inner(H7,State)
		h8 = numpy.inner(H8,State)
		
		""" PMU Current Measurements """
		h9 = numpy.inner(H9,State)
		h10 = numpy.inner(H10,State)
		
		""" WLS computation """
		H = numpy.concatenate((H1,H2,H3,H4,H5,H6,H7,H8,H9,H10),axis=0)
		y = numpy.concatenate((h1,h2,h3,h4,h5,h6,h7,h8,h9,h10),axis=0)
		res = numpy.subtract(z,y)
		g = numpy.inner(H.transpose(),numpy.inner(W,res))
		WH = numpy.inner(W,H.transpose())
		G = numpy.inner(H.transpose(),WH.transpose())
		
		Ginv = numpy.linalg.inv(G)
		Delta_State = numpy.inner(Ginv,g)
		
		State = State + Delta_State
		epsilon = numpy.amax(numpy.absolute(Delta_State))
		
		V.real = State[:len(system.nodes)]
		V.imag = State[len(system.nodes):]
		V = Vr + 1j*Vx
		
		num_iter = num_iter+1
	
	return V