from . import gain
from . import acoustics
from . import freqop
from . import geometry
from . import io
from . import signals
from . import plot
from .mesh import CAD as meshCAD
import numpy as np
from re import findall

## General functions
def findInArray(array, value):
    """
    Find the index and associated value of the closest data point in the array.

    Parameters
    ----------
    array : numpy array
        Array to be scanned through.

    value : float or complex or other data type
        Value to find in the array.

    Returns
    -------
    ind_value : int
        Index of the closest data point to the given value.

    value_out : float or complex or other data type
        Associated value in the input array that is closest to the given value.

    Notes
    -----
    This function calculates the index and the associated value of the closest data point
    in the input 'array' to the specified 'value'. It finds the index corresponding to the
    data point with the smallest absolute difference from the given value and returns both
    the index and the associated value from the array.
    """

    ind_value = int((np.abs(array - value)).argmin())
    value_out = array[ind_value]
    return ind_value, value_out


def normMinMax(dataset, minValue, maxValue):
    """
    Linearly normalize a dataset between specified minimum and maximum boundaries.

    Parameters
    ----------
    dataset : numpy array
        Dataset to be normalized.
    minValue : float
        Minimum boundary for normalization.
    maxValue : float
        Maximum boundary for normalization.

    Returns
    -------
    normDataset : numpy array
        Normalized dataset with values between minValue and maxValue.

    Notes
    -----
    This function linearly normalizes the input 'dataset' between the specified 'minValue'
    and 'maxValue' boundaries. It rescales the values of the dataset so that they span the
    range from 'minValue' to 'maxValue', while maintaining the relative differences between
    data points. The resulting 'normDataset' array contains the normalized dataset.
    """
    normDataset = (maxValue-minValue) * (dataset-np.min(dataset))/(np.max(dataset)-np.min(dataset))+minValue
    return normDataset

def wrap(phase):
    return np.arctan2(np.sin(phase), np.cos(phase))


def parallel(*elements):
    Y = 0
    for i in range(len(elements)):
        Y += 1/elements[i]
    Z = 1/Y
    return Z

def extract_numbers_to_list(input_string):
    # Use regular expression to find all numbers in the input string
    numbers = findall(r'\d+', input_string)

    # Convert the matched strings to actual integers and store them in a list
    numbers_list = [int(number) for number in numbers]

    return numbers_list

def slice_array_into_parts(arr, num_parts):
    if len(arr) % num_parts != 0:
        raise ValueError("Array length is not divisible by the number of parts")

    part_size = len(arr) // num_parts
    sliced_parts = [arr[i:i + part_size] for i in range(0, len(arr), part_size)]
    return sliced_parts