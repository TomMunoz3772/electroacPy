import numpy as np

## general gain calculation
def dB(x):
    """
    Convert a value to decibels (dB).

    Parameters
    ----------
    x : float or numpy array
        Value(s) to be converted to decibels.

    Returns
    -------
    dB_value : float or numpy array
        Value(s) converted to decibels (dB).

    Notes
    -----
    This function takes a numeric value or an array of values 'x' and converts
    them to their equivalent decibel (dB) representation using the formula
    20 * log10(abs(x)). The resulting dB value(s) provide a logarithmic
    representation of the input value(s) suitable for expressing ratios
    and differences in magnitude.
    """
    return 20*np.log10(np.abs(x))

def dB_10log(x):
    """
    Convert a value to decibels using the base 10 logarithm.

    Parameters
    ----------
    x : float or numpy array
        Value(s) to be converted to decibels using the base 10 logarithm.

    Returns
    -------
    dB_value : float or numpy array
        Value(s) converted to decibels using the base 10 logarithm.

    Notes
    -----
    This function takes a numeric value or an array of values 'x' and converts
    them to their equivalent decibel (dB) representation using the base 10
    logarithm, i.e., 10 * log10(x). The resulting dB value(s) provide a
    logarithmic representation of the input value(s) suitable for expressing
    ratios and differences in magnitude.
    """
    return 10 * np.log10(x)

def dB_ref(x, freq, ref):
    """
    Estimate the gain in decibels relative to a reference frequency.

    Parameters
    ----------
    x : numpy array
        Array of values to be normalized and converted to decibels.

    freq : numpy array
        Frequency axis corresponding to the values in 'x'.

    ref : float
        Reference frequency (Hz) to which the gain will be normalized.

    Returns
    -------
    gain_dB : numpy array
        Gain values in decibels relative to the reference frequency.

    Notes
    -----
    This function estimates the gain in decibels of the input 'x' array with respect to
    a reference frequency 'ref'. The input 'freq' should be the frequency axis corresponding
    to the values in 'x'. The gain is calculated by normalizing the input values to their
    maximum value at the reference frequency and then converting them to decibels.
    The resulting gain_dB values represent the gain of 'x' relative to the reference frequency.
    """
    indref = int((np.abs(freq - ref)).argmin())
    return 20*np.log10(np.abs(x)/np.max(np.abs(x[indref])))

def dB_zero(x):
    """
    Convert values to decibels relative to the maximum value.

    Parameters
    ----------
    x : numpy array
        Array of values to be normalized and converted to decibels.

    Returns
    -------
    gain_dB : numpy array
        Gain values in decibels relative to the maximum value.

    Notes
    -----
    This function converts the input 'x' array of values to decibels relative to its
    maximum value. The gain_dB values represent the relative magnitude of the input values
    expressed in decibels. The function normalizes the input values to their maximum value
    and then converts them to decibels using the formula 20 * log10(abs(x)/max(abs(x))).
    """
    return 20*np.log10(np.abs(x)/np.max(np.abs(x)))

## SPL calculation
def dBSPL(x):
    """
    Convert RMS signal values to dBSPL (decibels Sound Pressure Level).

    Parameters
    ----------
    x : numpy array
        Array of RMS signal values to be converted to dBSPL.

    Returns
    -------
    dBSPL_value : numpy array
        dBSPL values corresponding to the input RMS signal values.

    Notes
    -----
    This function converts the input 'x' array of RMS signal values to their equivalent
    values in dBSPL (decibels Sound Pressure Level). The conversion is performed using the
    formula 20 * log10(abs(x) / 2e-5), where 2e-5 Pa represents the reference sound pressure.
    The resulting dBSPL values provide a logarithmic representation of the RMS signal magnitude
    in decibels relative to the reference sound pressure level.
    """
    return 20*np.log10(np.abs(x)/2e-5)

def SPL(x):
    """
    Convert peak signal values to dBSPL (decibels Sound Pressure Level), RMS.

    Parameters
    ----------
    x : numpy array
        Array of peak signal values to be converted to dBSPL.

    Returns
    -------
    dBSPL_value : numpy array
        dBSPL values corresponding to the input peak signal values.

    Notes
    -----
    This function converts the input 'x' array of peak signal values to their equivalent
    values in dBSPL (decibels Sound Pressure Level). The conversion is performed using the
    formula 20 * log10(abs(x) / (2e-5 * sqrt(2))), where 2e-5 Pa represents the reference
    sound pressure and sqrt(2) is used to convert the peak values to RMS values.
    The resulting dBSPL values provide a logarithmic representation of the peak signal magnitude
    in decibels relative to the reference sound pressure level.
    """
    return 20*np.log10(np.abs(x)/2e-5/np.sqrt(2))

