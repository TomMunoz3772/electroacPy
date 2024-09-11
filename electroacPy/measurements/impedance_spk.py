import generalToolbox as gtb
import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, sin, tan, exp, pi
import scipy as sp
from scipy.fft import fft
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit



class impedance:
    def __init__(self, impedance_file, Bl, peak_max_bound=3000, skiprows=1, delimiter=",",
                 freqCol=0, HCol=1, PCol=2, usecols=None, decimate=1):
        self.impedance_file = impedance_file
        self.peak_max_bound = peak_max_bound


        # 1. load data
        data = np.loadtxt(impedance_file, skiprows=skiprows, delimiter=delimiter, usecols=usecols)
        self.freq = data[::decimate, freqCol]
        self.finterp = np.arange(self.freq[0], self.freq[-1] + 0.1, 0.1)
        try:
            self.Ze_meas = data[::decimate, HCol] * np.exp(1j * data[::decimate, PCol])  # H*e^(j*Phi)
        except:
            self.Ze_meas = data[::decimate, HCol]

        # Run a first optimization to fit electrical parameters
        # and some mechanical parameters
        self.params = self.get_first_opt()
        self.Zin = get_Zin(self.finterp, self.params['Re'], self.params['Le'],
                           self.params['Res'], self.params['R2'], self.params['L2'],
                           self.params['R3'], self.params['L3'], self.params['Ws'],
                           self.params['Qms'])

        # Run a second optimization to fit remaining mechanical parameters
        s = 2j * np.pi * self.finterp
        self.Zm = self.params['Res'] / (1 + self.params['Qms'] * (s / self.params['Ws'] + self.params['Ws'] / s))
        Mms, Rms, Cms = estimate_pmec(Bl, self.params['Qms'], self.params['Res'], self.params['fs'])
        self.params['Mms'] = Mms
        self.params['Cms'] = Cms
        self.params['Rms'] = Rms
        self.params['Bl'] = Bl
        # self.get_mec_opt()



    def get_first_opt(self):
        # 2. Get some visible values
        ind_bound, _ = gtb.findInArray(self.freq, self.peak_max_bound)
        Zmax = np.max(np.abs(self.Ze_meas[:ind_bound]))
        ind_fc, _ = gtb.findInArray(np.abs(self.Ze_meas), Zmax)
        fc = self.freq[ind_fc]

        # 3. Reduce the size of frequency range if it goes above 20 kHz
        ind_fmax, fmax = gtb.findInArray(self.freq, 20e3)
        freq = self.freq[:ind_fmax + 1]
        om = 2 * np.pi * freq
        s = 1j * om
        Ze_meas = self.Ze_meas[:ind_fmax + 1]

        # 4. Start Values
        Re = np.min(np.abs(self.Ze_meas))
        Le = 1e-5
        Res = Zmax - Re
        R2, L2 = 0.5, 2e-4
        R3, L3 = 0.25, 1e-4
        Wc, Qmc = 2 * np.pi * fc, getQMC(freq, Ze_meas)
        x0 = [Re, Le, Res, R2, L2, R3, L3, Wc, Qmc]  # initial parameters

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

        # impedance at resonance
        Zfc = np.abs(Re + sc * Le + gtb.parallel(sc * L2, R2) + gtb.parallel(sc * L3, R3))

        # # Q factors
        # Qec = np.sqrt(Mmc / Cmc) * (Zfc / Bl ** 2)
        # Qtc = Qmc * Qec / (Qmc + Qec)

        # out param
        params = {}
        params['Re'] = Re
        params['Le'] = Le
        params['Res'] = Res
        params['R2'] = R2
        params['L2'] = L2
        params['R3'] = R3
        params['L3'] = L3
        params['Ws'] = Wc
        params['fs'] = fc
        params['Qms'] = Qmc
        params['Zfs'] = Zfc
        # params['Cmc'] = Cmc
        # params['Mmc'] = Mmc
        # params['Zfc'] = Zfc
        # params['Qec'] = Qec
        # params['Qtc'] = Qtc
        return params

    # def get_mec_opt(self):
    #     Bl_est = estimate_Bl(self.finterp, self.Zm, self.Zin,
    #                          self.params['Res'], self.peak_max_bound)
    #     # print('estimated Bl: ', Bl_est)
    #     Mms_est, Rms_est, Cms_est = estimate_pmec(Bl_est, self.params['Qms'],
    #                                               self.params['Res'],
    #                                               self.params['fs'])
    #     # print('estimated Mms: ', Mms_est)
    #     # print('estimated Rms: ', Rms_est)
    #     # print('estimated Cms: ', Cms_est)
    #
    #     x1 = [Bl_est, Rms_est, Mms_est, Cms_est]
    #     Zems = self.Zin - (2j*np.pi*self.finterp * self.params['Le'] + self.params['Re'])
    #     popt2, pcov2 = sp.optimize.curve_fit(Zes_impedance, self.finterp, np.abs(self.Zm), p0=x1, maxfev=100000)
    #     self.params['Bl'] = popt2[0]
    #     self.params['Rms'] = popt2[1]
    #     self.params['Mms'] = popt2[2]
    #     self.params['Cms'] = popt2[3]
    #     return None

    def get_Zin(self, frequency):
        s = 2j*np.pi*frequency
        Zms = self.params['Rms'] + s*self.params['Mms'] + 1/s/self.params['Cms']
        Zm = self.params['Bl']**2 / Zms
        Ze = (self.params['Re'] + s*self.params['Le'] +
              gtb.parallel(self.params['R2'], s*self.params['L2']) +
              gtb.parallel(self.params['R3'], s*self.params['L3']))
        Zin = Ze+Zm
        return Zin

    def plot_Zin(self, measured=False, save=False):
        if measured is True:
            Z = self.Ze_meas
            freq = self.freq
            title='measured_impedance'
        elif measured is False:
            Z = self.Zin
            freq = self.finterp
            title = 'fitted_impedance'
        elif measured == "both":
            Z1 = self.Ze_meas
            freq_1 = self.freq
            Z2 = self.Zin
            freq_2 = self.finterp
            title = 'measured_fitted_impedance'

        if measured is True or measured is False:
            fig, ax = plt.subplots()
            ax.semilogx(freq, np.abs(Z), 'b',
                        label="modulus")
            ax.legend(loc='upper left')
            ax.grid(linestyle='dotted', which='both')
            ax.set(ylabel="|Z| [Ohm]", xlabel="Frequency [Hz]")
            ax2 = ax.twinx()
            ax2.semilogx(freq, np.angle(Z), ':r',
                         label="phase")
            ax2.legend(loc='upper right')
            ax2.set(ylabel="phase [rads]", ylim=[-pi / 2, pi / 2])
            ax2.grid()
            plt.tight_layout()
            if save is True:
                plt.savefig(title, dpi=256)

        elif measured == 'both':
            fig, ax = plt.subplots()
            ax.semilogx(freq_1, np.abs(Z1), 'b',
                        label="modulus - measured")
            ax.semilogx(freq_2, np.abs(Z2), 'g',
                        label="modulus - fitted")
            ax.legend(loc='upper left')
            ax.grid(linestyle='dotted', which='both')
            ax.set(ylabel="|Z| [Ohm]", xlabel="Frequency [Hz]")
            ax2 = ax.twinx()
            ax2.semilogx(freq_1, np.angle(Z1), ':r',
                         label="phase - measured", alpha=0.45)
            ax2.semilogx(freq_2, np.angle(Z2), ':', color='black',
                         label="phase - fitted")
            ax2.legend(loc='upper right')
            ax2.set(ylabel="phase [rads]", ylim=[-pi / 2, pi / 2])
            ax2.grid()
            plt.tight_layout()
            if save is True:
                plt.savefig(title, dpi=256)
        return None

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
    s = 1j * 2 * np.pi * freq

    # mechanical impedance
    fncZm = Res / (1 + Qmc * (s / Wc + Wc / s))

    # electrical impedance
    fncZe = Re + s * Le + gtb.parallel(s * L2, R2) + gtb.parallel(s * L3, R3)

    # total impedance
    fncZ = fncZe + fncZm
    return np.abs(fncZ)

def get_Zin(freq, Re, Le, Res, R2, L2, R3, L3, Wc, Qmc):
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
    s = 1j * 2 * np.pi * freq

    # mechanical impedance
    fncZm = Res / (1 + Qmc * (s / Wc + Wc / s))

    # electrical impedance
    fncZe = Re + s * Le + gtb.parallel(s * L2, R2) + gtb.parallel(s * L3, R3)

    # total impedance
    fncZ = fncZe + fncZm
    return fncZ

def getQMC(freq, Z, stepf=1e-3):
    """
    Find the quality factor of the impedance peak of a driver in a sealed box.
    :param freq:
    :param Z:
    :param stepf:
    :return:
    """
    # increase a bit the resolution (doesn't change curve's shape)
    finterp = np.arange(freq[0], freq[-1] + stepf, stepf)
    Zinterp = np.interp(finterp, freq, np.abs(Z))

    # get index of Zmax (index of f0) and Zmax
    ind_f0, Zmax = gtb.findInArray(np.abs(Zinterp), np.max(np.abs(Zinterp)))
    Ef0 = finterp[ind_f0]

    ind_fmin, Zfmin = gtb.findInArray(np.abs(Zinterp[:ind_f0]),
                                      np.max(np.abs(Zinterp)) * 0.707)  # get -3 dB value before resonance
    ind_fmax, Zfmax = gtb.findInArray(np.abs(Zinterp[ind_f0:]), np.max(np.abs(Zinterp)) * 0.707)
    range = [finterp[ind_fmin], finterp[ind_f0 + ind_fmax]]
    Qs = finterp[ind_f0] / (range[1] - range[0])
    return Qs

def Zes_impedance(freq, Bl, Rms, Mms, Cms):
    s = 2j * np.pi * freq
    Zems = Bl ** 2 / (Rms + s*Mms + 1/s/Cms)
    return np.abs(Zems)

def estimate_pmec(Bl, Qms, Res, fs):
    Mms = Qms * Bl**2 / 2 / np.pi / fs / Res
    Rms = Bl**2 / Res
    Cms = Res**2 / Bl**4 * Mms / Qms**2
    return Mms, Rms, Cms

def estimate_Bl(freq, Zm, Zin, Res, max_bound):
    ind_max, _ = gtb.findInArray(freq, max_bound)
    freq = freq[:ind_max]
    Zm = Zm[:ind_max]
    Zin = Zin[:ind_max]
    max_Zm = np.max(np.abs(Zm))
    max_Zin = np.max(np.abs(Zin))

    Bl_est = 0.01
    Zfc = max_Zin - max_Zm / Bl_est**2
    while Zfc < Res:
        Bl_est += 0.01
        Zfc = max_Zin - max_Zm / Bl_est
    return Bl_est