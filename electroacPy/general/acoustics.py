import numpy as np
import scipy as sp
from scipy.interpolate import interp1d

c = 343
rho = 1.22


## Thiele/Small param estimator
def getThieleSmallParam(impedance_file, Bl, peak_max_bound=3000, skiprows=1, delimiter=",", freqCol=0, HCol=1, PCol=2, usecols=None):
    """
    Get the Thiele/Small param from impedance curve and force factor Bl.

    Parameters:
    - impedance_file: .txt file
        containing the impedance as H*exp(1j*Phi) -> H column 2 and Phi column 3 (freq being 1)
    - Bl: int, float
        force factor of the motor

    Return:
    - Optimized Thiele/Small Parameters

    Note:
    Me just following matlab code; though the translation Matlab -> Python didn't work
    so I used curve_fit instead of fmin (I let Scipy do the error reduction because me dum dum)
    """
    from electroacPy.general import findInArray, parallel
    
    # 1. Load impedance data
    data = np.loadtxt(impedance_file, skiprows=skiprows, delimiter=delimiter, usecols=usecols)

    freq = data[:, freqCol]
    try:
        Ze_meas = data[:, HCol] * np.exp(1j * data[:, PCol]) # H*e^(j*Phi)
    except:
        Ze_meas = data[:, HCol]

    # 2. Get some visible values
    ind_bound, _ = findInArray(freq, peak_max_bound)
    Zmax = np.max(np.abs(Ze_meas[:ind_bound]))
    ind_fc, _ = findInArray(np.abs(Ze_meas), Zmax)
    fc = freq[ind_fc]

    # 3. Reduce the size of frequency range if it goes above 20 kHz
    ind_fmax, fmax = findInArray(freq, 20e3)
    freq = freq[:ind_fmax + 1]

    # 4. Start Values
    Re = np.min(np.abs(Ze_meas))
    Le = 1e-5
    Res = Zmax - Re
    R2, L2 = 0.5, 2e-4
    R3, L3 = 0.25, 1e-4
    Wc, Qmc = 2 * np.pi * fc, getQMC(freq, Ze_meas)
    x0 = [Re, Le, Res, R2, L2, R3, L3, Wc, Qmc] # initial parameters

    # 5. Merci Scipy
    popt, pcov = sp.optimize.curve_fit(sealed_impedance_calc,
                                       freq, np.abs(Ze_meas), p0=x0, maxfev=16000)
    Re, Le = popt[0], popt[1]
    Res = popt[2]
    R2, L2 = popt[3], popt[4]
    R3, L3 = popt[5], popt[6]
    Wc, Qmc = popt[7], popt[8]
    fc = Wc / 2 / np.pi
    sc = 1j * Wc

    # -----------------
    # 6. Let's finalize
    Mmc = Qmc / Wc * Bl**2/Res
    Cmc = 1 / Wc**2 / Mmc

    # impedance at resonance
    Zfc = np.abs(Re + sc*Le + parallel(sc*L2, R2) + parallel(sc*L3, R3))

    # Q factors
    Qec = np.sqrt(Mmc/Cmc) * (Zfc / Bl**2)
    Qtc = Qmc*Qec/(Qmc+Qec)

    # out param
    params = {}
    params['Re']  = Re
    params['Le']  = Le
    params['Res'] = Res
    params['R2']  = R2
    params['L2']  = L2
    params['R3']  = R3
    params['L3']  = L3
    params['Wc']  = Wc
    params['fc']  = fc
    params['Qmc'] = Qmc
    params['Cmc'] = Cmc
    params['Mmc'] = Mmc
    params['Zfc'] = Zfc
    params['Qec'] = Qec
    params['Qtc'] = Qtc
    return params

def getQMC(freq, Z, stepf=1e-3):
    """
    Find the quality factor of the impedance peak of a driver in a sealed box.
    :param freq:
    :param Z:
    :param stepf:
    :return:
    """
    from electroacPy.general import findInArray

    # increase a bit the resolution (doesn't change curve's shape)
    finterp = np.arange(freq[0], freq[-1] + stepf, stepf)
    Zinterp = np.interp(finterp, freq, np.abs(Z))

    # get index of Zmax (index of f0) and Zmax
    ind_f0, Zmax = findInArray(np.abs(Zinterp), np.max(np.abs(Zinterp)))
    Ef0 = finterp[ind_f0]

    ind_fmin, Zfmin = findInArray(np.abs(Zinterp[:ind_f0]),
                                      np.max(np.abs(Zinterp)) * 0.707)  # get -3 dB value before resonance
    ind_fmax, Zfmax = findInArray(np.abs(Zinterp[ind_f0:]), np.max(np.abs(Zinterp)) * 0.707)
    range = [finterp[ind_fmin], finterp[ind_f0 + ind_fmax]]
    Qs = finterp[ind_f0] / (range[1] - range[0])
    return Qs


def sealed_impedance_calc(freq, Re, Le, Res, R2, L2, R3, L3, Wc, Qmc):
    """
    sealed impedance calculator, should be use within the curve_fit() algorithm
    :param freq:
    :param Re:
    :param Le:
    :param Res:
    :param R2:
    :param L2:
    :param R3:
    :param L3:
    :param Wc:
    :param Qmc:
    :return:
    """
    from electroacPy.general import parallel

    s = 1j * 2 * np.pi * freq

    # mechanical impedance
    fncZm = Res / (1 + Qmc * (s / Wc + Wc / s))

    # electrical impedance
    fncZe = Re + s * Le + parallel(s * L2, R2) + parallel(s * L3, R3)

    # total impedance
    fncZ = fncZe + fncZm
    return np.abs(fncZ)

def Zes_impedance(freq, Bl, Rms, Mms, Cms):
    s = 2j * np.pi * freq
    Zems = Bl ** 2 / (Rms + s*Mms + 1/s/Cms)
    return np.abs(Zems)

def estimate_pmec(Bl, Qms, Res, fs):
    Mms = Qms * Bl**2 / 2 / np.pi / fs / Res
    Rms = Bl**2 / Res
    Cms = Res**2 / Bl**4 * Mms / Qms**2
    return Mms, Rms, Cms

## Room Acoustics
def art60(rt60, surface, volume, f=1000, c=343):
    """
    Compute the absorption coefficient based on the given reverberation time (RT60),
    total surface area of absorbing materials in the room (surface), and the volume
    of the room (volume). Compute for f=1000 Hz.

    Parameters:
    - rt60 (float): Reverberation time in seconds.
    - surface (float): Total surface area of absorbing materials in square meters.
    - volume (float): Volume of the room in cubic meters.

    Returns:
    float: Absorption coefficient calculated based on the provided inputs.

    Example:
    >>> art60(1.5, 50, 100)
    """

    om = 2*np.pi*f
    k = om*c
    a = k * volume / rt60 / surface
    return a


def a2z(absorptionCoeff, c=343, rho=1.22):
    Z_air = rho * c
    Zn = Z_air * (1-np.sqrt(1-absorptionCoeff)) / (1+np.sqrt(1-absorptionCoeff))
    return Zn


def rectangularRoomModes(Lx, Ly, Lz, N):
    """
    Compute the Nth room modes of a rectangular room given its dimensions.

    Parameters:
    - Lx (float): Length of the room in the x-direction.
    - Ly (float): Width of the room in the y-direction.
    - Lz (float): Height of the room in the z-direction.
    - N (int): Amount of modes to compute (N^3).

    Returns:
    - modes (list): A list containing the Nth room modes in the format (fx, fy, fz),
                   where fx, fy, and fz are the frequencies along the x, y, and z directions, respectively.

    Note:
    - Room modes represent the resonant frequencies at which standing waves can exist in the room.
    - The mode numbers (N) correspond to the number of half-wavelengths along each dimension.
    - The formula for the Nth room mode frequencies is given by:
        f[x, y, z] = 343/2 * ((nx/Lx)^2 + (ny/Ly)^2 + (nz/Lz)^2)^0.5
    """

    f = np.zeros([N, N, N])
    for nx in range(N):
        for ny in range(N):
            for nz in range(N):
                f[nx, ny, nz] = 343/2 * np.sqrt((nx/Lx)**2+(ny/Ly)**2+(nz/Lz)**2)
    return f


def interpolate_Impedance(impedanceData, freq_in, freq_out):
    """
    Interpolate impedance data from freq_in array to freq_out array.

    Parameters:
    - impedanceData: numpy array, impedance values corresponding to freq_in.
    - freq_in: numpy array, input frequency array.
    - freq_out: numpy array, output frequency array for interpolation.

    Returns:
    - interpolatedImpedance: numpy array, interpolated impedance values corresponding to freq_out.
    """

    # Create an interpolation function
    interpolation_function = interp1d(freq_in, impedanceData, kind='linear', fill_value="extrapolate")

    # Use the interpolation function to compute impedance values at freq_out
    interpolatedImpedance = interpolation_function(freq_out)

    return interpolatedImpedance


#%% Export tools
def export_directivity(folder_name, file_prefix,
                       frequency_array, angle_array, pmic_array, frd=False):
    """
    Export a directivity evaluation extracted from a loudspeakerSystem object.
    Can be imported in other software, such as VituixCAD or VACS.

    Parameters
    ----------
    folder_name : str
        export path.
    frequency_array : numpy array
        frequency range of length Nfft.
    angle_array : numpy array
        angular position of microphones
    pmic_array : numpy array
        acoustic pressure (complex) of shape [Nfft, N_angles].

    Returns
    -------
    folder with *.txt file containing frequency, SPL and phase of extracted 
    directivity.

    """
    import os
    from electroacPy.general.gain import SPL
    
    # Create the folder if it doesn't exist
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    
    if frd is False:
        extension = "txt"
    elif frd is True:
        extension = "frd"
    
    # Loop through each microphone (angle)
    for i, angle in enumerate(angle_array):
        # Compute SPL values for the current microphone
        spl_values = SPL(pmic_array[:, i])
        
        # Compute phase (in degrees) using np.angle and np.rad2deg
        phase_values = np.rad2deg(np.angle(pmic_array[:, i]))
        
        # Prepare the data to save: frequency, SPL, and phase
        data_to_save = np.column_stack((frequency_array, spl_values, phase_values))
        
        # File name format: folder_name_angle.txt
        file_name = os.path.join(folder_name, f"{file_prefix}_{angle}.{extension}")
        
        # Save the data to the text file using np.savetxt with the correct header
        np.savetxt(file_name, data_to_save, header="frequency  SPL  phase", fmt="%.3e")
        
        # print(f"Saved data to {file_name}")