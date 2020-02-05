import numpy as np
from .network import BusType
from .results import Results


def solve(system):
    """It performs powerflow by using rectangular node voltage state variables and considering the current mismatch function.
    
    Solve the non-linear powerflow problem stated by

    r = z-h(state) = 0

    following the Newton-Raphson approach

    delta_state = H^-1 * r
    new_state = old_state + delta_state

    r: residual function (current mismatch) 
    z: expected currents
    state: rectangular voltages (i.e. [V0_re, V1_re, ..., VN_re, V0_im, V1_im, ... , VN_im])
    h: currents calculated from state
    H: Jacobian matrix
    V: same as state but with complex numbers (i.e. [V0_re+j*V0_im, V1_re+j*V1_im, ...])
    """

    nodes_num = system.get_nodes_num()
    z = np.zeros(2 * nodes_num)
    h = np.zeros(2 * nodes_num)
    H = np.zeros((2 * nodes_num, 2 * nodes_num))

    for node in system.nodes:
        if node.ideal_connected_with =='':
            i = node.index
            m = 2 * i
            i2 = i + nodes_num
            node_type = node.type
            if node_type == BusType.SLACK:
                z[m] = np.real(node.voltage_pu)
                z[m + 1] = np.imag(node.voltage_pu)
                H[m][i] = 1
                H[m + 1][i2] = 1
            elif node_type is BusType.PQ:
                H[m][i] = - np.real(system.Ymatrix[i][i])
                H[m][i2] = np.imag(system.Ymatrix[i][i])
                H[m + 1][i] = - np.imag(system.Ymatrix[i][i])
                H[m + 1][i2] = - np.real(system.Ymatrix[i][i])
                idx1 = np.subtract(system.Adjacencies[i], 1)
                idx2 = idx1 + nodes_num
                H[m][idx1] = - np.real(system.Ymatrix[i][idx1])
                H[m][idx2] = np.imag(system.Ymatrix[i][idx1])
                H[m + 1][idx1] = - np.imag(system.Ymatrix[i][idx1])
                H[m + 1][idx2] = - np.real(system.Ymatrix[i][idx1])
            elif node_type is BusType.PV:
                z[m + 1] = np.real(node.power)
                H[m][i] = - np.real(system.Ymatrix[i][i])
                H[m][i2] = np.imag(system.Ymatrix[i][i])
                idx1 = np.subtract(system.Adjacencies[i], 1)
                idx2 = idx1 + nodes_num
                H[m][idx1] = - np.real(system.Ymatrix[i][idx1])
                H[m][idx2] = np.imag(system.Ymatrix[i][idx1])

    epsilon = 10 ** (-10)
    #epsilon = 0.01
    diff = 5
    V = np.ones(nodes_num) + 1j * np.zeros(nodes_num)
    num_iter = 0

    state = np.concatenate((np.ones(nodes_num), np.zeros(nodes_num)), axis=0)

    while diff > epsilon:
        for node in system.nodes:
            if node.ideal_connected_with =='':
                i = node.index
                m = 2 * i
                i2 = i + nodes_num
                node_type = node.type
                if node_type is BusType.SLACK:
                    h[m] = np.inner(H[m], state)
                    h[m + 1] = np.inner(H[m + 1], state)
                elif node_type is BusType.PQ:
                    z[m] = (np.real(node.power_pu) * np.real(V[i]) + 
                           np.imag(node.power_pu) * np.imag(V[i])) / (np.abs(V[i]) ** 2)
                    z[m + 1] = (np.real(node.power_pu) * np.imag(V[i]) - 
                           np.imag(node.power_pu) * np.real(V[i])) / (np.abs(V[i]) ** 2)
                    h[m] = np.inner(H[m], state)
                    h[m + 1] = np.inner(H[m + 1], state)
                elif node_type is BusType.PV:
                    z[m] = (np.real(node.power_pu) * np.real(V[i]) + 
                            np.imag(node.power_pu) * np.imag(V[i]))(np.abs(V[i]) ** 2)
                    h[m] = np.inner(H[m], state)
                    h[m + 1] = np.abs(V[i])
                    H[m + 1][i] = np.cos(np.angle(V[i]))
                    H[m + 1][i2] = np.sin(np.angle(V[i]))

        r = np.subtract(z, h)
        Hinv = np.linalg.inv(H)
        delta_state = np.inner(Hinv, r)
        state = state + delta_state
        diff = np.amax(np.absolute(delta_state))

        V = state[:nodes_num] + 1j * state[nodes_num:]
        num_iter = num_iter + 1

    # calculate all the other quantities of the grid
    powerflow_results = Results(system)
    powerflow_results.load_voltages(V)
    powerflow_results.calculate_all()

    return powerflow_results, num_iter
