import numpy as np

def find_interpolant(nodes,data):
    '''
    :param  NODES: ARRAY of the nodes, e.g. the nodes along the z-axis and x_axis
    :param DATA: ARRAY of the aerodynamic loading's magnitude. It can be seen as the value of each node.

    The LENGTH of both INPUTS have to be the same

    :return: MATRIX (N x 4) with the INTERPOLANTS, it return the 1st column the a_i, the 2nd column the b_i, the 3rd column the c_i, the 4th column the d_i
    '''
    #Number of points
    N = len(nodes)

    # Computes the function difference (f_{i+1}-f_{i})
    diff_magnitude = np.diff(data)

    # Computes the step-size h_{i}=x_{i+1}-x_{i}
    diff_nodes = np.diff(nodes)

    # Computes the h_{i-1}/6*M_{i-1}+(h_{i-1}+h_{i})/3*M_{i}+h_{i}/6*M_{i+1}=(f_{i+1}-f_{i})\h_{i}(f_{i}-f_{i-1})/h_{i-1}
    a_matrix = np.zeros((N - 2, N - 2))  # left hand side
    b_matrix = np.zeros(N - 2)  # right hand side

    for row in range(0, N - 2):
        a_matrix[row, row] = (diff_nodes[row] + diff_nodes[row + 1]) / 3
        if row != N - 3:
            a_matrix[row, row + 1] = diff_nodes[row + 1] / 6
        if row != 0:
            a_matrix[row, row - 1] = diff_nodes[row] / 6

        b_matrix[row] = (data[row + 2] - data[row + 1]) / diff_nodes[row + 1] - (
                    data[row + 1] - data[row]) / diff_nodes[row]

    # Calculate the Constraint Matrix (M_0,M_1,M_2) --> A*x=B
    m_matrix = np.linalg.tensorsolve(a_matrix, b_matrix)

    # Adds M_{0} and M_{n} equal to 0 (just adding two elements, one at the beginning and one at the end)
    m_matrix = np.insert(np.append(m_matrix, [0]), 0, 0)

    # Computes coefficients a_{i}, b_{i}, c_{i} and d_{i}
    coeff_a_matrix = np.diff(m_matrix) / (6 * diff_nodes)

    coeff_b_matrix = m_matrix / 2
    coeff_b_matrix = np.delete(coeff_b_matrix, -1)

    coeff_c_matrix = np.zeros(N - 1)
    for element in range(0, N - 1):
        coeff_c_matrix[element] = diff_magnitude[element] / diff_nodes[element] - diff_nodes[element] / 3 * m_matrix[element] - \
                                  diff_nodes[element] / 6 * m_matrix[element + 1]

    coeff_d_matrix = np.delete(np.array(data),-1)

    # Coeff_matrix is a matrix containing all the coefficients
    coeff_matrix = np.transpose(np.array([coeff_a_matrix, coeff_b_matrix, coeff_c_matrix, coeff_d_matrix]))
    return coeff_matrix

def compute_interpolant(nodes,interpolant,resolution):
    '''
    :param NODES: ARRAY of the nodes, e.g. the nodes along the z-axis and x_axis
    :param INTERPOLANTS: MATRIX (N x 4) with the INTERPOLANTS, it return the 1st column the a_i, the 2nd column the b_i, the 3rd column the c_i, the 4th column the d_i
    :param  RESOLUTION: Minimum spacing needed in millimiters [mm]
    :return: 2 ARRAYS: -NEW_NODES
                       -NEW_LOADING
    '''

    # Computes the step-size h_{i}=x_{i+1}-x_{i}
    diff_nodes = np.diff(nodes)

    #Finds the maximum spacing between nodes
    max_diff_nodes = np.amax(np.abs(diff_nodes))

    #Number of steps needed between two MAIN NODES
    N = np.ceil(max_diff_nodes/(resolution*10**(-3)))+1

    new_nodes = []
    new_loading = []

    for element in range(0, interpolant.shape[0]):
        node_i_spline = nodes[element]
        diff_node_spline = diff_nodes[element]
        step_spline = np.linspace(0, diff_node_spline, N)  # Taking N steps for each spline
        for step in range(len(step_spline) - 1):
            spline_function = interpolant[element, 0] * (step_spline[step]) ** 3 + interpolant[element, 1] * (
            step_spline[step]) ** 2 + interpolant[element, 2] * (step_spline[step]) + interpolant[element, 3]
            new_loading.append(spline_function)
            new_nodes.append(node_i_spline + step_spline[step])

    new_nodes = np.array(new_nodes)
    new_loading = np.array(new_loading)

    return new_nodes,new_loading