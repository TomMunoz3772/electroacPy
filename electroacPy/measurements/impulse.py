#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:47:02 2023

@author: tom.munoz
"""
import numpy as np
from os import listdir
from os.path import isfile, join
import re

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
    return theta, IR_MATRIX
