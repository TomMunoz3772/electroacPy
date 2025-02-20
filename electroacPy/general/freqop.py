"""
Functions to manipulate frequency related problems (decimate, octave smoothing)
"""

import numpy as np
from tqdm import tqdm
import pyuff

## manipulate frequency array (Hz)
def laplace(f):
    """
    Convert frequency values from Hertz to complex frequency in the Laplace domain.

    Parameters
    ----------
    f : numpy array
        Frequency axis in Hertz.

    Returns
    -------
    s : numpy array
        Complex frequency values corresponding to the Laplace domain (2*pi*freq*1j).

    Notes
    -----
    This function takes an array of frequency values in Hertz and converts them to
    complex frequency values in the Laplace domain using the formula 2*pi*freq*1j,
    where 'freq' is the input frequency values. The resulting 's' values can be used
    for Laplace domain analysis of systems.
    """
    return 2 * np.pi * f * 1j


def octave_smoothing(X_FFT, freq_data, nthOctBand):
    """
    Apply octave band smoothing to frequency response data.

    Parameters
    ----------
    X_FFT : numpy array
        Dataset to be smoothed (Frequency Response in SPL).

    freq_data : numpy array
        Frequency axis corresponding to the data.

    nthOctBand : int
        Order of the octave band to use for filtering.

    Returns
    -------
    smoothed_data : numpy array
        Frequency response data after applying octave band smoothing.

    Notes
    -----
    This function applies octave band smoothing to the input frequency response data
    'X_FFT'. The smoothing is performed using the specified 'nthOctBand' order of
    octave band filtering. The 'freq_data' parameter is the frequency axis corresponding
    to the data. The function calculates and applies the smoothing coefficients based on
    the octave band filtering approach, resulting in smoothed frequency response data.
    """
    f_min = np.argwhere(freq_data > 1)[0][0]  # first value superior to zero
    x_fft = np.zeros(len(X_FFT))  # data to be smoothed
    if nthOctBand > 0:
        print('Filtering Data, 1/{} octave band'.format(nthOctBand))
        for i in tqdm(np.arange(f_min, len(freq_data))):
            fv = freq_data[i]
            s = (fv / nthOctBand) / np.pi
            g = np.exp(-(((freq_data - fv) ** 2) / (2 * (s ** 2))))
            g_s = g / np.sum(g)
            x_fft[i] = np.sum(g_s * X_FFT)  # smoothed coefficients
        if X_FFT.all() >= 0:
            x_fft[x_fft < 0] = 0  # set all negative data to zero
    return x_fft

def smooth_directivity(X_FFT, freq_data, nthOctBand):
    """
    Apply octave band smoothing to frequency response data.

    Parameters
    ----------
    X_FFT : numpy array
        Dataset to be smoothed (Frequency Response in SPL).

    freq_data : numpy array
        Frequency axis corresponding to the data.

    nthOctBand : int
        Order of the octave band to use for filtering.

    Returns
    -------
    smoothed_data : numpy array
        Frequency response data after applying octave band smoothing.

    Notes
    -----
    This function applies octave band smoothing to the input frequency response data
    'X_FFT'. The smoothing is performed using the specified 'nthOctBand' order of
    octave band filtering. The 'freq_data' parameter is the frequency axis corresponding
    to the data. The function calculates and applies the smoothing coefficients based on
    the octave band filtering approach, resulting in smoothed frequency response data.
    """
    f_min = np.argwhere(freq_data > 1)[0][0]  # first value superior to zero
    x_fft = np.zeros(np.shape(X_FFT))  # data to be smoothed
    if nthOctBand > 0:
        print('Filtering Data, 1/{} octave band'.format(nthOctBand))
        for i in tqdm(np.arange(f_min, len(freq_data))):
            fv = freq_data[i]
            s = (fv / nthOctBand) / np.pi
            g = np.exp(-(((freq_data - fv) ** 2) / (2 * (s ** 2))))
            g_s = g / np.sum(g)
            for j in range(len(X_FFT)):
                x_fft[j, i] = np.sum(g_s * X_FFT[j, :])  # smoothed coefficients
        if X_FFT.all() >= 0:
            x_fft[x_fft < 0] = 0  # set all negative data to zero
    return x_fft

def decimate_frequency_axis(old_axis, new_axis):
    """
    Decimate an existing frequency axis into a new one by finding the closest values.

    Parameters:
    old_axis (numpy.ndarray): The existing frequency axis from which values will be selected.
    new_axis (numpy.ndarray): The new frequency axis to decimate onto.

    Returns:
    numpy.ndarray: The decimated frequency axis containing values 
    from the old_axis that are closest to the values in new_axis.

    Example:
        >>> old_axis = np.linspace(0, 100, 1000)
        >>> new_axis = np.logspace(1, 2, 50)
        >>> decimated_result = decimate_frequency_axis(old_axis, new_axis)
        >>> print(decimated_result)
    """
    decimated_axis = np.empty_like(new_axis)
    closest_idx = np.zeros(len(new_axis), dtype=int)

    for i, new_val in enumerate(new_axis):
        closest_idx_tmp = np.argmin(np.abs(old_axis - new_val))
        decimated_axis[i] = old_axis[closest_idx_tmp]
        closest_idx[i] = int(closest_idx_tmp)

    return decimated_axis, closest_idx


## create frequency arrays
def freq_log10(start, stop, numpoints, Nd=2):
    """
    Generate a logarithmically spaced frequency array.

    Parameters
    ----------
    start : float
        Starting frequency (lowest frequency) in Hertz.

    stop : float
        Stopping frequency (highest frequency) in Hertz.

    numpoints : int
        Number of points in the generated frequency array.
    
    Nd : int
        Number of decimals

    Returns
    -------
    freq_array : numpy array
        Logarithmically spaced frequency array spanning from 'start' to 'stop'.

    Notes
    -----
    This function generates a numpy array of 'numpoints' logarithmically spaced
    frequency values between 'start' and 'stop'. The frequency values are evenly
    distributed on a logarithmic scale, making them suitable for logarithmic plots
    or analyses where the frequency range spans multiple orders of magnitude.
    """
    freq_array = np.logspace(np.log10(start), np.log10(stop), numpoints)
    return np.round(freq_array, Nd)


def freq_uff(uff_file, freq_axis_user):
    """
    Create frequency axis from uff dataset. Will return closest frequencies to freq_axis_user
    :param uff_file:
    :param freq_axis_wanted:
    :return:
    """
    uff_data          = pyuff.UFF(uff_file)
    data              = uff_data.read_sets()
    freq_axis_meas    = data[-1]['x']
    freq_out, ind_out = decimate_frequency_axis(freq_axis_meas, freq_axis_user)
    return np.unique(freq_out)
