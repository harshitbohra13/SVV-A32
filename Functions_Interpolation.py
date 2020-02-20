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
    The length of NODES and the ROWS of the INTERPOLANTS have to be equal
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
            if element==interpolant.shape[0]-1:
                step+=1
                new_nodes.append(nodes[-1])
                spline_function = interpolant[element, 0] * (step_spline[step]) ** 3 + interpolant[element, 1] * (
                step_spline[step]) ** 2 + interpolant[element, 2] * (step_spline[step]) + interpolant[element, 3]
                new_loading.append(spline_function)

    new_nodes = np.array(new_nodes)
    new_loading = np.array(new_loading)

    return new_nodes,new_loading

def compute_value_interpolation(nodes,interpolant,position):
    '''
    The length of NODES and the ROWS of the INTERPOLANTS have to be equal
    :param NODES: array of nodes, e.g. the nodes along the z-axis or x-axis
    :param INTERPOLANTS: MATRIX (N x 4) with the INTERPOLANTS, it return the 1st column the a_i, the 2nd column the b_i, the 3rd column the c_i, the 4th column the d_i
    :param POSITION: the exact x-axis or z-axis location that you would like to get the value
    :return:-the VALUE (INTIGER) of the interpolation at a such location
            -the ARRAY containing the coefficients at that specific position
    '''

    #Number of nodes
    N = len(nodes)

    location_spline = 0

    #Finding in which spline corresponds the location of the POSITION
    for element in range(0,N):
        if nodes[element]<position<nodes[element+1] or nodes[element]>position>nodes[element+1] or nodes[element]==position:
            location_spline = element
            break;

    #Computing RESULT using the cubic spline formula
    diff_position=position-nodes[element]
    a_coeff = interpolant[location_spline,0]
    b_coeff = interpolant[location_spline, 1]
    c_coeff = interpolant[location_spline, 2]
    d_coeff = interpolant[location_spline, 3]
    coeff_array = np.array([a_coeff,b_coeff,c_coeff,d_coeff])
    result = a_coeff*diff_position**3+b_coeff*diff_position**2+c_coeff*diff_position+d_coeff
    return result,coeff_array

def integrate_polynomial(coefficients_array):
    '''
    :param coefficients_array: [a_i,b_i,c_i,d_i] = a_i*x**3+b_i*x**2+c_i*x+d
    :return: array integrated [a_i/4,b_i/3,c_i/2,d_i,0 'it should be constant'] = a_i/4*x**4+b_i/3*x**3+c_i/2*x**2+d*x+Constant
    '''

    #Number of coefficients
    n = len(coefficients_array)

    #Matrix containing 1,1/2,1/3,1/4...
    a_matrix = np.zeros((n,n))
    for element in range(0,n):
        a_matrix[n-1-element][n-1-element] = 1/(element+1)

    x_matrix = np.transpose(np.array([coefficients_array]))

    #Computing the result array and adding a 0 at the end
    b_matrix = np.transpose(np.dot(a_matrix,x_matrix))
    b_matrix = np.append(b_matrix,[0])
    return b_matrix