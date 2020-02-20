import numpy as np
from pyvolt.results import Results
from pyvolt.measurement import *


def DsseCall(system, measurements):
    """
    Performs state estimation
    It identifies the type of measurements present in the measurement set and
    calls the appropriate estimator for dealing with them.

    @param system: model of the system (nodes, lines, topology)
    @param measurements: Vector of measurements in Input (voltages, currents, powers)
    return: object of class results.Results
    """

    # select type of Estimator.
    # if at least a PMU is present we launch the combined estimator, 
    # otherwise a simple traditional etimator
    Vmag_meas = 0
    Vpmu_meas = 0
    for elem in measurements.measurements:
        if elem.meas_type == MeasType.V_mag:
            Vmag_meas += 1
        elif elem.meas_type == MeasType.Vpmu_mag:
            Vpmu_meas += 1

    trad_code = 1 if Vmag_meas > 0 else 0
    PMU_code = 2 if Vpmu_meas > 0 else 0
    est_code = trad_code + PMU_code

    # number of nodes of the grid
    nodes_num = system.get_nodes_num()

    Gmatrix = system.Ymatrix.real
    Bmatrix = system.Ymatrix.imag
    Yabs_matrix = np.absolute(system.Ymatrix)
    Yphase_matrix = np.angle(system.Ymatrix)
    Adj = system.Adjacencies

    # run Estimator.
    if est_code == 1:
        Vest = DsseTrad(nodes_num, measurements, Gmatrix, Bmatrix, Yabs_matrix, Yphase_matrix, Adj)
    elif est_code == 2:
        Vest = DssePmu(nodes_num, measurements, Gmatrix, Bmatrix, Adj)
    else:
        Vest = DsseMixed(nodes_num, measurements, Gmatrix, Bmatrix, Yabs_matrix, Yphase_matrix, Adj)

    # calculate all the other quantities of the grid
    results = Results(system)
    results.load_voltages(Vest)
    results.calculate_all()

    return results


def DsseTrad(nodes_num, measurements, Gmatrix, Bmatrix, Yabs_matrix, Yphase_matrix, Adj):
    """
    Traditional state estimator
    It performs state estimation using rectangular node voltage state variables
    and it is customized to work without PMU measurements

    @param nodes_num: number of nodes of the grid
    @param measurements: Vector of measurements in Input (voltages, currents, powers)
    @param Gmatrix: conductance matrix 
    @param Bmatrix: susceptance matrix
    @param Yabs_matrix: magnitude of the admittance matrix
    @param Yphase_matrix: phase of the admittance matrix
    @param Adj: 
    return: np.array V - estimated voltages
    """

    # calculate  weightsmatrix (obtained as stdandard_deviations^-2)
    weights = measurements.getWeightsMatrix()
    W = np.diag(weights)

    # Jacobian for Power Injection Measurements
    H2, H3 = calculateJacobiMatrixSinj(measurements, nodes_num, Gmatrix, Bmatrix, Adj, type=2)

    # Jacobian for branch Power Measurements
    H4, H5 = calculateJacobiBranchPower(measurements, nodes_num, Gmatrix, Bmatrix, type=2)

    # get array which contains the index of measurements type V_mag and I_mag
    vidx = measurements.getIndexOfMeasurements(type=MeasType.V_mag)
    iidx = measurements.getIndexOfMeasurements(type=MeasType.I_mag)
    nvi = len(vidx)
    nii = len(iidx)

    # get array which contains the index of measurements type MeasType.Sinj_real, MeasType.Sinj_imag in the array measurements.measurements
    pidx = measurements.getIndexOfMeasurements(type=MeasType.Sinj_real)
    qidx = measurements.getIndexOfMeasurements(type=MeasType.Sinj_imag)

    # get array which contains the index of measurements type MeasType.S_real, MeasType.S_imag in the array measurements.measurements
    p1br = measurements.getIndexOfMeasurements(type=MeasType.S1_real)
    p2br = measurements.getIndexOfMeasurements(type=MeasType.S2_real)
    q1br = measurements.getIndexOfMeasurements(type=MeasType.S1_imag)
    q2br = measurements.getIndexOfMeasurements(type=MeasType.S2_imag)

    # get an array with all measured values (affected by uncertainty)
    z = measurements.getMeasValues()

    V = np.ones(nodes_num) + 1j * np.zeros(nodes_num)
    State = np.concatenate((np.ones(nodes_num), np.zeros(nodes_num)), axis=0)
    epsilon = 5
    num_iter = 0

    # Iteration of Netwon Rapson method: needed to solve non-linear system of equation
    while epsilon > 10 ** (-6):
        """ Computation of equivalent current measurements in place of the power measurements """
        # in every iteration the input power measurements are converted into currents by dividing by the voltage estimated at the previous iteration
        z = convertSinjMeasIntoCurrents(measurements, V, z, pidx, qidx)
        z = convertSbranchMeasIntoCurrents(measurements, V, z, p1br, q1br, p2br, q2br)

        """ Voltage Magnitude Measurements """
        h1, H1 = update_h1_vector(measurements, V, vidx, nvi, nodes_num, type=2)

        """ Power Injection Measurements """
        # h(x) vector where power injections are present
        h2 = np.inner(H2, State)
        h3 = np.inner(H3, State)

        """ Power Flow Measurements """
        # h(x) vector where power flows are present
        h4 = np.inner(H4, State)
        h5 = np.inner(H5, State)

        """ Current Magnitude Measurements """
        h6, H6 = update_h6_vector(measurements, V, iidx, nii, Yabs_matrix, Yphase_matrix, nodes_num, num_iter, type=2)

        """ WLS computation """
        # all the sub-matrixes of H calcualted so far are merged in a unique matrix
        H = np.concatenate((H1, H2, H3, H4, H5, H6), axis=0)
        # h(x) sub-vectors are concatenated
        y = np.concatenate((h1, h2, h3, h4, h5, h6), axis=0)
        # "res" is the residual vector. The difference between input measurements and h(x)
        res = np.subtract(z, y)
        # g = transpose(H) * W * res
        g = np.inner(H.transpose(), np.inner(W, res))
        WH = np.inner(W, H.transpose())
        # G is the gain matrix, that will have to be inverted at each iteration
        G = np.inner(H.transpose(), WH.transpose())
        # inversion of G
        Ginv = np.linalg.inv(G)
        # Delta includes the updates of the states for the current Newton Rapson iteration
        Delta_State = np.inner(Ginv, g)
        # state is updated
        State = State + Delta_State
        # calculate the NR treeshold (for the next while check)
        epsilon = np.amax(np.absolute(Delta_State))
        # update the voltages
        V.real = State[:nodes_num]
        V.imag = State[nodes_num:]
        V.imag = np.concatenate(([0], V.imag), axis=0)

        num_iter = num_iter + 1

    return V


def DssePmu(nodes_num, measurements, Gmatrix, Bmatrix, Adj):
    """
    Traditional state estimator
    It performs state estimation using rectangular node voltage state variables
    and it is customized to work using PMU measurements

    @param nodes_num: number of nodes of the grid
    @param measurements: Vector of measurements in Input (voltages, currents, powers)
    @param Gmatrix
    @param Bmatrix
    @param Adj
    return: np.array V - estimated voltages
    """
    # calculate weights matrix (obtained as stdandard_deviations^-2)
    weights = measurements.getWeightsMatrix()
    W = np.diag(weights)

    # Jacobian for Power Injection Measurements
    H2, H3 = calculateJacobiMatrixSinj(measurements, nodes_num, Gmatrix, Bmatrix, Adj, type=1)

    # Jacobian for branch Power Measurements
    H4, H5 = calculateJacobiBranchPower(measurements, nodes_num, Gmatrix, Bmatrix, type=1)

    # Jacobian for branch Power Measurements
    H7, H8 = calculateJacobiVoltagePmu(measurements, nodes_num, Gmatrix, Bmatrix)
    W = update_W_matrix(measurements, weights, W, "Vpmu")

    # Jacobian for Current Pmu Measurements
    H9, H10 = calculateJacobiCurrentPmu(measurements, nodes_num, Gmatrix, Bmatrix)
    W = update_W_matrix(measurements, weights, W, "Ipmu")

    # get an array with all measured values (affected by uncertainty)
    z = measurements.getMeasValues()

    H = np.concatenate((H2, H3, H4, H5, H7, H8, H9, H10), axis=0)
    WH = np.inner(W, H.transpose())
    G = np.inner(H.transpose(), WH.transpose())
    Ginv = np.linalg.inv(G)

    # get array which contains the index of measurements type MeasType.Sinj_real, MeasType.Sinj_imag in the array measurements.measurements
    pidx = measurements.getIndexOfMeasurements(type=MeasType.Sinj_real)
    qidx = measurements.getIndexOfMeasurements(type=MeasType.Sinj_imag)

    # get array which contains the index of measurements type MeasType.S_real, MeasType.S_imag in the array measurements.measurements
    p1br = measurements.getIndexOfMeasurements(type=MeasType.S1_real)
    p2br = measurements.getIndexOfMeasurements(type=MeasType.S2_real)
    q1br = measurements.getIndexOfMeasurements(type=MeasType.S1_imag)
    q2br = measurements.getIndexOfMeasurements(type=MeasType.S2_imag)

    # get an array with all measured values (affected by uncertainty)
    z = measurements.getMeasValues()

    V = np.ones(nodes_num) + 1j * np.zeros(nodes_num)
    State = np.concatenate((np.ones(nodes_num), np.zeros(nodes_num)), axis=0)
    epsilon = 5
    num_iter = 0

    # Iteration of Netwon Rapson method: needed to solve non-linear system of equation
    while epsilon > 10 ** (-6):
        """ Computation of equivalent current measurements in place of the power measurements """
        # in every iteration the input power measurements are converted into currents by dividing by the voltage estimated at the previous iteration
        z = convertSinjMeasIntoCurrents(measurements, V, z, pidx, qidx)
        z = convertSbranchMeasIntoCurrents(measurements, V, z, p1br, q1br, p2br, q2br)

        """ WLS computation """
        y = np.inner(H, State)
        res = np.subtract(z, y)
        g = np.inner(H.transpose(), np.inner(W, res))

        Delta_State = np.inner(Ginv, g)

        State = State + Delta_State
        epsilon = np.amax(np.absolute(Delta_State))

        V.real = State[:nodes_num]
        V.imag = State[nodes_num:]

        num_iter = num_iter + 1

    return V


def DsseMixed(nodes_num, measurements, Gmatrix, Bmatrix, Yabs_matrix, Yphase_matrix, Adj):
    """
    Traditional state estimator
    It performs state estimation using rectangular node voltage state variables
    and it is built to work in scenarios where both conventional and PMU measurements
    are simultaneously present.

    @param nodes_num: number of nodes of the grid
    @param measurements: Vector of measurements in Input (voltages, currents, powers)
    @param Gmatrix
    @param Bmatrix
    @param Yabs_matrix
    @param Yphase_matrix
    @param Adj
    return: np.array V (estimated voltages)
    """

    # calculate weights matrix (obtained as stdandard_deviations^-2)
    weights = measurements.getWeightsMatrix()
    W = np.diag(weights)

    # Jacobian Matrix. Includes the derivatives of the measurements (voltages, currents, powers) with respect to the states (voltages)

    # Jacobian for Power Injection Measurements
    H2, H3 = calculateJacobiMatrixSinj(measurements, nodes_num, Gmatrix, Bmatrix, Adj, type=1)

    # Jacobian for branch Power Measurements
    H4, H5 = calculateJacobiBranchPower(measurements, nodes_num, Gmatrix, Bmatrix, type=1)

    # Jacobian for branch Power Measurements
    H7, H8 = calculateJacobiVoltagePmu(measurements, nodes_num, Gmatrix, Bmatrix)
    W = update_W_matrix(measurements, weights, W, "Vpmu")

    # Jacobian for Current Pmu Measurements
    H9, H10 = calculateJacobiCurrentPmu(measurements, nodes_num, Gmatrix, Bmatrix)
    W = update_W_matrix(measurements, weights, W, "Ipmu")

    # get array which contains the index of measurements type V_mag and I_mag
    vidx = measurements.getIndexOfMeasurements(type=MeasType.V_mag)
    iidx = measurements.getIndexOfMeasurements(type=MeasType.I_mag)
    nvi = len(vidx)
    nii = len(iidx)

    # get array which contains the index of measurements type MeasType.Sinj_real, MeasType.Sinj_imag in the array measurements.measurements
    pidx = measurements.getIndexOfMeasurements(type=MeasType.Sinj_real)
    qidx = measurements.getIndexOfMeasurements(type=MeasType.Sinj_imag)

    # get array which contains the index of measurements type MeasType.S_real, MeasType.S_imag in the array measurements.measurements
    p1br = measurements.getIndexOfMeasurements(type=MeasType.S1_real)
    p2br = measurements.getIndexOfMeasurements(type=MeasType.S2_real)
    q1br = measurements.getIndexOfMeasurements(type=MeasType.S1_imag)
    q2br = measurements.getIndexOfMeasurements(type=MeasType.S2_imag)

    # get an array with all measured values (affected by uncertainty)
    z = measurements.getMeasValues()

    V = np.ones(nodes_num) + 1j * np.zeros(nodes_num)
    State = np.concatenate((np.ones(nodes_num), np.zeros(nodes_num)), axis=0)
    epsilon = 5
    num_iter = 0

    # Iteration of Netwon Rapson method: needed to solve non-linear system of equation
    while epsilon > 10 ** (-6):
        """ Computation of equivalent current measurements in place of the power measurements """
        z = convertSinjMeasIntoCurrents(measurements, V, z, pidx, qidx)
        z = convertSbranchMeasIntoCurrents(measurements, V, z, p1br, q1br, p2br, q2br)

        """ Voltage Magnitude Measurements """
        h1, H1 = update_h1_vector(measurements, V, vidx, nvi, nodes_num, type=1)

        """ Power Injection Measurements """
        h2 = np.inner(H2, State)
        h3 = np.inner(H3, State)

        """ Power Flow Measurements """
        h4 = np.inner(H4, State)
        h5 = np.inner(H5, State)

        """ Current Magnitude Measurements """
        h6, H6 = update_h6_vector(measurements, V, iidx, nii, Yabs_matrix, Yphase_matrix, nodes_num, num_iter, type=1)

        """ PMU Voltage Measurements """
        h7 = np.inner(H7, State)
        h8 = np.inner(H8, State)

        """ PMU Current Measurements """
        h9 = np.inner(H9, State)
        h10 = np.inner(H10, State)

        """ WLS computation """
        H = np.concatenate((H1, H2, H3, H4, H5, H6, H7, H8, H9, H10), axis=0)
        y = np.concatenate((h1, h2, h3, h4, h5, h6, h7, h8, h9, h10), axis=0)
        res = np.subtract(z, y)
        g = np.inner(H.transpose(), np.inner(W, res))
        WH = np.inner(W, H.transpose())
        G = np.inner(H.transpose(), WH.transpose())

        Ginv = np.linalg.inv(G)
        Delta_State = np.inner(Ginv, g)

        State = State + Delta_State
        epsilon = np.amax(np.absolute(Delta_State))

        V.real = State[:nodes_num]
        V.imag = State[nodes_num:]

        num_iter = num_iter + 1

    return V


def calculateJacobiMatrixSinj(measurements, nodes_num, Gmatrix, Bmatrix, Adj, type):
    """
    It calculates the Jacobian for Power Injection Measurements
    (converted to equivalent rectangualar current measurements)

    @param measurements: object of class Measurement that contains all measurements (voltages, currents, powers)
    @param nodes_num: len of system.nodes
    @param Gmatrix
    @param Bmatrix
    @param Adj
    @param type: 1 for DssePmu and DssePmu, 2 for DsseTrad
    return: 1. H2: Jacobian for Pinj
            2. H3: Jacobian for Qinj
    """
    # get all measurements of type MeasType.Sinj_real
    pinj_meas = measurements.getMeasurements(type=MeasType.Sinj_real)
    # get all measurements of type MeasType.Sinj_real
    qinj_meas = measurements.getMeasurements(type=MeasType.Sinj_imag)
    if type == 1:
        H2 = np.zeros((len(pinj_meas), 2 * nodes_num))
        H3 = np.zeros((len(qinj_meas), 2 * nodes_num))
    elif type == 2:
        H2 = np.zeros((len(pinj_meas), 2 * nodes_num - 1))
        H3 = np.zeros((len(qinj_meas), 2 * nodes_num - 1))

    for index, measurement in enumerate(pinj_meas):
        m = measurement.element.index
        if type == 1:
            m2 = m + nodes_num
        elif type == 2:
            m2 = m + nodes_num - 1
        H2[index][m] = - Gmatrix[m][m]
        H2[index][m2] = Bmatrix[m][m]
        H3[index][m] = - Bmatrix[m][m]
        H3[index][m2] = - Gmatrix[m][m]
        idx = np.subtract(Adj[m], 1)
        H2[index][idx] = - Gmatrix[m][idx]
        H3[index][idx] = - Bmatrix[m][idx]
        if type == 1:
            idx2 = idx + nodes_num
        elif type == 2:
            if 0 in idx:
                pos = np.where(idx == 0)
                idx = np.delete(idx, pos)
            idx2 = idx + nodes_num - 1
        H2[index][idx2] = Bmatrix[m][idx]
        H3[index][idx2] = - Gmatrix[m][idx]

    return H2, H3


def calculateJacobiBranchPower(measurements, nodes_num, Gmatrix, Bmatrix, type):
    """
    It calculates the Jacobian for branch Power Measurements
    (converted to equivalent rectangualar current measurements)

    @param measurements: object of class Measurement that contains all measurements (voltages, currents, powers)
    @param nodes_num: len of system.nodes
    @param Gmatrix
    @param Bmatrix
    @param type: 1 for DssePmu and DssePmu, 2 for DsseTrad
    return:	1. H4: Jacobian for S_real
            2. H5: Jacobian for S_imag
    """

    # get all measurements of type MeasType.S_real and MeasType.S_imag
    p1_meas = measurements.getMeasurements(type=MeasType.S1_real)
    p2_meas = measurements.getMeasurements(type=MeasType.S2_real)
    q1_meas = measurements.getMeasurements(type=MeasType.S1_imag)
    q2_meas = measurements.getMeasurements(type=MeasType.S2_imag)

    if type == 1:
        H4 = np.zeros((len(p1_meas) + len(p2_meas), 2 * nodes_num))
        H5 = np.zeros((len(q1_meas) + len(q2_meas), 2 * nodes_num))
    elif type == 2:
        H4 = np.zeros((len(p1_meas) + len(p2_meas), 2 * nodes_num - 1))
        H5 = np.zeros((len(q1_meas) + len(q2_meas), 2 * nodes_num - 1))

    for i, measurement in enumerate(p1_meas):
        m = measurement.element.start_node.index
        n = measurement.element.end_node.index
        H4[i][m] = - Gmatrix[m][n]
        H4[i][n] = Gmatrix[m][n]
        H5[i][m] = - Bmatrix[m][n]
        H5[i][n] = Bmatrix[m][n]
        if type == 1:
            m2 = m + nodes_num
            H4[i][m2] = Bmatrix[m][n]
            H5[i][m2] = - Gmatrix[m][n]
            n2 = n + nodes_num
            H4[i][n2] = - Bmatrix[m][n]
            H5[i][n2] = Gmatrix[m][n]
        elif type == 2:
            if m > 0:
                m2 = m + nodes_num - 1
                H4[i][m2] = Bmatrix[m][n]
                H5[i][m2] = - Gmatrix[m][n]
            if n > 0:
                n2 = n + nodes_num - 1
                H4[i][n2] = - Bmatrix[m][n]
                H5[i][n2] = Gmatrix[m][n]

    for i, measurement in enumerate(iterable=p2_meas, start=len(p1_meas)):
        n = measurement.element.start_node.index
        m = measurement.element.end_node.index
        H4[i][m] = - Gmatrix[m][n]
        H4[i][n] = Gmatrix[m][n]
        H5[i][m] = - Bmatrix[m][n]
        H5[i][n] = Bmatrix[m][n]
        if type == 1:
            m2 = m + nodes_num
            H4[i][m2] = Bmatrix[m][n]
            H5[i][m2] = - Gmatrix[m][n]
            n2 = n + nodes_num
            H4[i][n2] = - Bmatrix[m][n]
            H5[i][n2] = Gmatrix[m][n]
        elif type == 2:
            if m > 0:
                m2 = m + nodes_num - 1
                H4[i][m2] = Bmatrix[m][n]
                H5[i][m2] = - Gmatrix[m][n]
            if n > 0:
                n2 = n + nodes_num - 1
                H4[i][n2] = - Bmatrix[m][n]
                H5[i][n2] = Gmatrix[m][n]

    return H4, H5


def calculateJacobiVoltagePmu(measurements, nodes_num, Gmatrix, Bmatrix):
    """
    It calculates the Jacobian for Voltage Pmu Measurements
    (converted to equivalent rectangualar current measurements)

    @param measurements: object of class Measurement that contains all measurements (voltages, currents, powers)
    @param nodes_num: len of system.nodes
    @param Gmatrix
    @param Bmatrix
    return: 1. H7: Jacobian for S_real
            2. H8: Jacobian for S_imag
    """

    # get all measurements of type MeasType.Vpmu_mag
    Vpmu_mag_meas = measurements.getMeasurements(type=MeasType.Vpmu_mag)
    # get all measurements of type MeasType.Vpmu_phase
    Vpmu_phase_meas = measurements.getMeasurements(type=MeasType.Vpmu_phase)
    H7 = np.zeros((len(Vpmu_mag_meas), 2 * nodes_num))
    H8 = np.zeros((len(Vpmu_mag_meas), 2 * nodes_num))

    # TODO: index of Vmag = index of Vphase???
    for index, measurement in enumerate(Vpmu_mag_meas):
        vamp = measurement.meas_value_ideal
        vtheta = Vpmu_phase_meas[index].meas_value_ideal
        m = measurement.element.index
        H7[index][m] = 1
        m2 = m + nodes_num
        H8[index][m2] = 1

    return H7, H8


def calculateJacobiCurrentPmu(measurements, nodes_num, Gmatrix, Bmatrix):
    """
    It calculates the Jacobian for Current Pmu Measurements
    (converted to equivalent rectangualar current measurements)

    @param measurements: object of class Measurement that contains all measurements (voltages, currents, powers)
    @param nodes_num: len of system.nodes
    @param Gmatrix
    @param Bmatrix
    return: 1. H9: Jacobian for S_real
            2. H10: Jacobian for S_imag
    """

    # get all measurements of type MeasType.Ipmu_mag
    Ipmu_mag_meas = measurements.getMeasurements(type=MeasType.Ipmu_mag)
    # get all measurements of type MeasType.Vpmu_phase
    Ipmu_phase_meas = measurements.getMeasurements(type=MeasType.Vpmu_phase)
    H9 = np.zeros((len(Ipmu_mag_meas), 2 * nodes_num))
    H10 = np.zeros((len(Ipmu_mag_meas), 2 * nodes_num))

    for index, measurement in enumerate(Ipmu_mag_meas):
        iamp = measurement.meas_value_ideal
        itheta = Ipmu_phase_meas[index].meas_value_ideal
        m = measurement.element.start_node.index
        n = measurement.element.end_node.index
        H9[index][m] = - Gmatrix[m][n]
        H9[index][n] = Gmatrix[m][n]
        H10[index][m] = - Bmatrix[m][n]
        H10[index][n] = Bmatrix[m][n]
        m2 = m + nodes_num
        n2 = n + nodes_num
        H9[index][m2] = Bmatrix[m][n]
        H9[index][n2] = - Bmatrix[m][n]
        H10[index][m2] = - Gmatrix[m][n]
        H10[index][n2] = Gmatrix[m][n]

    return H9, H10


def update_W_matrix(measurements, weights, W, type):
    """
    adds to the matrix W the values related to pmu measurements (Vpmu or Ipmu)

    @param measurements: object of class Measurement that contains all measurements (voltages, currents, powers)
    @param weights: weights matrix
    @param W: np.diag(weights)
    @param type: "Vpmu" or "Ipmu"
    return: updated W matrix
    """

    if type == "Vpmu":
        # get index of all measurements of type "MeasType.Vpmu_mag" in the array MeasurementSet.measurements
        index_mag = measurements.getIndexOfMeasurements(MeasType.Vpmu_mag)
        # get index of all measurements of type "MeasType.Vpmu_phase" in the array MeasurementSet.measurements
        index_phase = measurements.getIndexOfMeasurements(MeasType.Vpmu_phase)
    elif type == "Ipmu":
        # get index of all measurements of type "MeasType.Ipmu_mag" in the array MeasurementSet.measurements
        index_mag = measurements.getIndexOfMeasurements(MeasType.Ipmu_mag)
        # get index of all measurements of type "MeasType.Ipmu_phase" in the array MeasurementSet.measurements
        index_phase = measurements.getIndexOfMeasurements(MeasType.Ipmu_phase)

    for index, (idx_mag, idx_theta) in enumerate(zip(index_mag, index_phase)):
        value_amp = measurements.measurements[idx_mag].meas_value
        value_theta = measurements.measurements[idx_theta].meas_value
        rot_mat = np.array([[np.cos(value_theta), - value_amp * np.sin(value_theta)],
                            [np.sin(value_theta), value_amp * np.cos(value_theta)]])
        starting_cov = np.array([[weights[idx_mag], 0], [0, weights[idx_theta]]])
        final_cov = np.inner(rot_mat, np.inner(starting_cov, rot_mat.transpose()))
        W[idx_mag][idx_mag] = final_cov[0][0]
        W[idx_theta][idx_theta] = final_cov[1][1]
        W[idx_mag][idx_theta] = final_cov[0][1]
        W[idx_theta][idx_mag] = final_cov[1][0]

    return W


def update_h1_vector(measurements, V, vidx, nvi, nodes_num, type):
    """
    update h1 and H1 vectors

    @param measurements: Vector of measurements in Input (voltages, currents, powers)
    @param V: vector of the estimated voltages
    @param vidx: array which contains the index of measurements type V_mag in measurements.measurements
    @param nvi: len of vidx
    @param nodes_num: number of nodes of the grid - len(system.nodes)
    @param type: 1 for DssePmu and DssePmu, 2 for DsseTrad
    return: vector h1 and H1
    """

    # at every iteration we update h(x) vector where V measure are available
    h1 = np.zeros(nvi)
    # the Jacobian rows where voltage measurements are presents is updated
    if type == 1:
        H1 = np.zeros((nvi, 2 * nodes_num))
    elif type == 2:
        H1 = np.zeros((nvi, 2 * nodes_num - 1))
    for i, index_vmag in enumerate(vidx):
        # get index of the node
        node_index = measurements.measurements[index_vmag].element.index
        h1[i] = np.absolute(V[node_index])
        H1[i][node_index] = np.cos(np.angle(V[node_index]))
        if type == 1:
            m2 = node_index + nodes_num
            H1[i][m2] = np.sin(np.angle(V[node_index]))
        elif type == 2:
            if m > 0:
                m2 = node_index + nodes_num - 1
                H1[i][m2] = np.sin(V.phase[m])

    return h1, H1


def update_h6_vector(measurements, V, iidx, nii, Yabs_matrix, Yphase_matrix, nodes_num, num_iter, type):
    """
    update h6 and H6 vectors where current flows are present

    @param measurements: Vector of measurements in Input (voltages, currents, powers)
    @param V: vector of the estimated voltages
    @param iidx: array which contains the index of measurements type I_mag in measurements.measurements
    @param nii: len of iidx
    @param Yabs_matrix:
    @param Yphase_matrix:
    @param nodes_num: number of nodes of the grid - len(system.nodes)
    @param num_iter: number of current iteration
    @param type: 1 for DssePmu and DssePmu, 2 for DsseTrad
    return: vector h6 and H6
    """

    # Current Magnitude Measurements
    h6re = np.zeros((nii))
    h6im = np.zeros((nii))
    h6complex = np.zeros((nii), dtype=complex)
    h6 = np.ones((nii))
    if type == 1:
        H6 = np.zeros((nii, 2 * nodes_num))
    elif type == 2:
        H6 = np.zeros((nii, 2 * nodes_num - 1))

    for i, index_imag in enumerate(iidx):
        # get index of the start node
        m = measurements.measurements[index_imag].element.start_node.index
        # get index of the end node
        n = measurements.measurements[index_imag].element.end_node.index
        h6re[i] = Yabs_matrix[m][n] * (
                (V[n].real - V[m].real) * np.cos(Yphase_matrix[m][n]) + (V[m].imag - V[n].imag) * np.sin(
            Yphase_matrix[m][n]))
        h6im[i] = Yabs_matrix[m][n] * (
                (V[n].real - V[m].real) * np.sin(Yphase_matrix[m][n]) + (V[n].imag - V[m].imag) * np.cos(
            Yphase_matrix[m][n]))
        h6complex[i] = h6re[i] + 1j * h6im[i]
        if num_iter > 0:
            h6[i] = np.absolute(h6complex[i])
        H6[i][m] = - Yabs_matrix[m][n] * (
                np.cos(Yphase_matrix[m][n]) * h6re[i] + np.sin(Yphase_matrix[m][n]) * h6im[i]) / h6[i]
        H6[i][n] = Yabs_matrix[m][n] * (np.cos(Yphase_matrix[m][n]) * h6re[i] + np.sin(Yphase_matrix[m][n]) * h6im[i]) / \
                   h6[i]
        if type == 1:
            m2 = m + nodes_num
            H6[i][m2] = - Yabs_matrix[m][n] * (
                    np.cos(Yphase_matrix[m][n]) * h6im[i] - np.sin(Yphase_matrix[m][n]) * h6re[i]) / h6[i]
            n2 = n + nodes_num
            H6[i][n2] = Yabs_matrix[m][n] * (
                    np.cos(Yphase_matrix[m][n]) * h6im[i] - np.sin(Yphase_matrix[m][n]) * h6re[i]) / h6[i]
        if type == 2:
            if m > 0:
                m2 = m + nodes_num - 1
                H6[i][m2] = - Yabs_matrix[m][n] * (
                        np.cos(Yphase_matrix[m][n]) * h6im[i] - np.sin(Yphase_matrix[m][n]) * h6re[i]) / h6[i]
            if n > 0:
                n2 = n + nodes_num - 1
                H6[i][n2] = Yabs_matrix[m][n] * (
                        np.cos(Yphase_matrix[m][n]) * h6im[i] - np.sin(Yphase_matrix[m][n]) * h6re[i]) / h6[i]

    return h6, H6


def convertSinjMeasIntoCurrents(measurements, V, z, pidx, qidx):
    """
    In every iteration the input power measurements are converted into currents
    by dividing by the voltage estimated at the previous iteration and this values
    are replaced in the array z

    @param measurements: Vector of measurements in Input (voltages, currents, powers)
    @param V: vector of the estimated voltages
    @param z: array with all measured values (affected by uncertainty) --> MeasurementSet.getMeasValues
    @param pidx: array which contains the index of measurements type Sinj_real in measurements.measurements
    @param qidx: array which contains the index of measurements type Sinj_imag in measurements.measurements
    returns: updated z array
    """

    for p_index, q_index in zip(pidx, qidx):
        # get values of the measurements p_inj and q_inj   (affected by uncertainty-->meas_value)
        p_inj = measurements.measurements[p_index].meas_value
        q_inj = measurements.measurements[q_index].meas_value
        # get index of the node
        node_index = measurements.measurements[
            p_index].element.index  # == measurements.measurements[q_index].element.index
        z[p_index] = (p_inj * V[node_index].real + q_inj * V[node_index].imag) / (np.absolute(V[node_index]) ** 2)
        z[q_index] = (p_inj * V[node_index].imag - q_inj * V[node_index].real) / (np.absolute(V[node_index]) ** 2)

    return z


def convertSbranchMeasIntoCurrents(measurements, V, z, p1br, q1br, p2br, q2br):
    """
    In every iteration the input power measurements are converted into currents
    by dividing by the voltage estimated at the previous iteration and this values
    are replaced in the array z

    @param measurements: Vector of measurements in Input (voltages, currents, powers)
    @param V: vector of the estimated voltages
    @param z: array with all measured values (affected by uncertainty) --> MeasurementSet.getMeasValues
    @param p1br: array which contains the index of measurements type S1_real in measurements.measurements
    @param q1br: array which contains the index of measurements type S1_imag in measurements.measurements
    @param p2br: array which contains the index of measurements type S2_real in measurements.measurements
    @param q2br: array which contains the index of measurements type S2_imag in measurements.measurements
    returns: updated z array
    """
    for pbr_index, qbr_index in zip(p1br, q1br):
        # get values of the measurements pbr_inj and qbr_inj   (affected by uncertainty-->meas_value)
        p_br = measurements.measurements[pbr_index].meas_value
        q_br = measurements.measurements[qbr_index].meas_value
        # get index of the start node
        node_index = measurements.measurements[
            pbr_index].element.start_node.index  # == measurements.measurements[qbr_index].element.start_node.index
        z[pbr_index] = (p_br * V[node_index].real + q_br * V[node_index].imag) / (np.absolute(V[node_index]) ** 2)
        z[qbr_index] = (p_br * V[node_index].imag - q_br * V[node_index].real) / (np.absolute(V[node_index]) ** 2)

    for pbr_index, qbr_index in zip(p2br, q2br):
        # get values of the measurements pbr_inj and qbr_inj   (affected by uncertainty-->meas_value)
        p_br = measurements.measurements[pbr_index].meas_value
        q_br = measurements.measurements[qbr_index].meas_value
        # get index of the start node
        node_index = measurements.measurements[
            pbr_index].element.end_node.index  # == measurements.measurements[qbr_index].element.start_node.index
        z[pbr_index] = (p_br * V[node_index].real + q_br * V[node_index].imag) / (np.absolute(V[node_index]) ** 2)
        z[qbr_index] = (p_br * V[node_index].imag - q_br * V[node_index].real) / (np.absolute(V[node_index]) ** 2)

    return z
