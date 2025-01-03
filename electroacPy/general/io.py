"""
Functions to import / export datasets (mostly import for now)
"""

import numpy as np
from copy import copy
from os import listdir
from os.path import isfile, join
import re
import pyuff
from scipy.signal.windows import tukey


def load_impulse(folderPath, isConvertedFromMAT=False):
    """
    Import impulse response data from CSV files.

    Parameters
    ----------
    folderPath : str
        Path to the folder where IR.csv files are stored.

    isConvertedFromMAT : bool, optional
        If True, indicates that the data has been converted from .mat to .csv format,
        which may affect the number of rows to skip when loading data. Default is False.

    Returns
    -------
    theta : numpy array
        Array of angles in degrees corresponding to the imported impulse responses.

    IR_MATRIX : numpy array
        Impulse response matrix of shape (len(theta), len(IR)).

    Notes
    -----
    This function imports impulse response data from CSV files stored in the specified
    'folderPath'. The function returns a tuple containing two numpy arrays: 'theta' and
    'IR_MATRIX'. The 'theta' array holds the angles in degrees corresponding to the
    imported impulse responses. The 'IR_MATRIX' array holds the actual impulse response
    data in a matrix format, where each row represents an angle and each column represents
    a time step in the impulse response.

    If the 'isConvertedFromMAT' parameter is set to True, the function adjusts the number of
    rows to skip when loading data to accommodate the conversion from .mat to .csv format.
    """
    files = [f for f in listdir(folderPath) if isfile(join(folderPath, f))]
    files = [f for f in files if "IR" in f]

    # load file to get length of time data
    if isConvertedFromMAT is False:
        dummyData = np.loadtxt(join(folderPath, files[0]), skiprows=4, delimiter=',')
        lenData = len(dummyData[:, 0])
        skiprow = 4
    elif isConvertedFromMAT is True:
        dummyData = np.loadtxt(join(folderPath, files[0]), skiprows=1, delimiter=',')
        lenData = len(dummyData[:, 0])
        skiprow = 1

    # build angle array
    pattern = r'[\d\.]+'
    numerical_values = re.findall(pattern, files[0][:-4])

    where_theta = files[0].index(numerical_values[0])

    theta = np.zeros(len(files))
    for t in range(len(theta)):
        # print("where is theta: ", where_theta)
        # print("theta: ", theta[t])
        theta[t] = float(files[t][where_theta:where_theta + 3])
    theta = np.sort(theta)

    # build IR matrix
    IR_MATRIX = np.zeros([len(theta), lenData])
    for t in range(len(theta)):
        filepath = join(folderPath, files[0][:where_theta] + str(theta[t]) + '_' + 'IR.csv')
        data = np.loadtxt(filepath, skiprows=skiprow, delimiter=',')
        time = data[:, 0]
        signal = data[:, 1]
        IR_MATRIX[t, :] = signal
    return theta, time, IR_MATRIX


def processIR(IR_matrix, Fs, minTime, maxTime, alpha=0.02):
    indMin = int(minTime * Fs)
    indMax = int(maxTime * Fs)
    indTot = indMax - indMin
    window = tukey(int(indTot), alpha)

    IR_out = np.zeros([len(IR_matrix), indTot])

    for i in range(len(IR_matrix)):
        IR_out[i, :] = IR_matrix[i, indMin:indMax] * window
    return IR_out

def processFRF(IR_matrix, Fs, Nfft=None):
    if Nfft == None:
        Nfft = np.shape(IR_matrix)[1]
    FRF_out = np.zeros([len(IR_matrix), Nfft], dtype=complex)
    for i in range(len(IR_matrix)):
        FRF_out[i, :] = np.fft.fft(IR_matrix[i, :], Nfft)
    freq = np.arange(Nfft) * Fs / Nfft
    return freq, FRF_out

# load acceleration data from uff file
def loadUFF(uffFile):
    """
    Load acceleration / velocity data from a UFF file and return frequency axis, acceleration data,
    and point positions in space.

    Parameters
    ----------
    uffFile : str
        Path to the *.uff file.

    Returns
    -------
    uff_data : tuple
        A tuple containing the following elements:
        - freq : numpy array
            Frequency axis corresponding to the loaded acceleration data.
        - a_matrix : numpy array, shape (n_acc_points, nfft)
            Matrix of complex-valued acceleration data.
        - X : numpy array, shape (nPoints, 3)
            Matrix containing point positions in 3D space (x, y, z).

    Notes
    -----
    This function loads acceleration data from a UFF (Universal File Format) file and extracts
    relevant information. The function returns a tuple 'uff_data' containing three elements:

    - 'freq': A numpy array representing the frequency axis corresponding to the loaded
      acceleration data.
    - 'a_matrix': A numpy array with shape (n_acc_points, nfft), containing complex-valued
      acceleration data. Each row corresponds to a different point in space.
    - 'X': A numpy array with shape (nPoints, 3) representing the point positions in 3D space.
      Each row represents a point, and the columns correspond to the coordinates (x, y, z).
    """
    file = pyuff.UFF(uffFile)
    data = file.read_sets()

    # create acceleration matrix
    nData = len(data)

    point_count_acc = 0
    point_count_vel = 0
    index_acc = []
    index_vel = []
    index_58 = []
    for i in range(nData):
        if data[i]['type'] == 58:
            if data[i]['ordinate_axis_lab'] == 'Acceleration' and data[i]['orddenom_axis_lab'] == 'Voltage':
                # print('here')
                point_count_acc += 1
                index_acc.append(i)
            if data[i]['ordinate_axis_lab'] == 'Velocity':
                point_count_vel += 1
                index_vel.append(i)
            index_58.append(i)

    nfft = len(data[index_58[0]]['data'])
    a_matrix = np.zeros([point_count_acc, nfft], dtype=complex)
    v_matrix = np.zeros([point_count_vel, nfft], dtype=complex)
    # node_id = np.zeros(len(index_58))
    for i, val in enumerate(index_acc):
        a_matrix[i, :] = data[val]['data']
    for i, val in enumerate(index_vel):
        v_matrix[i, :] = data[val]['data']
        # node_id[i] = data[val]['rsp_node']

    freq = data[index_58[0]]['x']

    # space data
    point_count = 0
    index_2411 = []
    for i in range(len(data)):
        if data[i]['type'] == 2411:
            point_count += 1
            index_2411.append(i)

    x, y, z = data[index_2411[0]]['x'], \
        data[index_2411[0]]['y'], \
        data[index_2411[0]]['z']
    X = np.zeros([len(x), 3])
    X[:, 0] = x
    X[:, 1] = y
    X[:, 2] = z

    # output data
    uff_data = (freq, a_matrix, v_matrix, X)
    return uff_data


def loadREW_FRF(file, skiprows=14):
    data = np.loadtxt(file, skiprows=skiprows)
    F = data[:, 0]
    Hf = data[:, 1]
    Ht = data[:, 2]
    return F, Hf, Ht