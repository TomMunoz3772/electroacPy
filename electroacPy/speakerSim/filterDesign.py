#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:55:45 2023

@author: tom.munoz
"""

import numpy as np
import matplotlib.pyplot as plt
import electroacPy.general as gtb
from scipy.signal import butter, freqs, freqz
from numpy import sin, cos, sqrt, pi
import scipy as sp

class xovers:
    """
    A class for designing analog crossover filters using biquad filter stages.

    This class provides methods to add various types of analog crossover filters to a design,
    including low-pass, high-pass, band-pass, peak EQ, low-shelf, high-shelf, gain, and delay filters.
    It also supports deleting filters, updating the frequency response, and plotting the overall response.

    Attributes:
        f_axis (array-like): The frequency axis for which the filters are designed.
        h (array-like, dtype=complex): The current total frequency response of the crossover.
        filterList (list): List of filter names in the crossover.
        frf_List (list): List of individual filter frequency responses.

    Methods:
        addLowPass(filterName, order, fc): Add an analog low-pass filter.
        addHighPass(filterName, order, fc): Add an analog high-pass filter.
        addBandPass(filterName, order, fc1, fc2): Add an analog band-pass filter.
        addLowPassBQ(filterName, fc, Q, dBGain=0, coeff=False): Add a biquad low-pass filter.
        addHighPassBQ(filterName, fc, Q, dBGain=0, coeff=False): Add a biquad high-pass filter.
        addBandPassBQ(filterName, fc, Q, dBGain=0, coeff=False): Add a biquad band-pass filter.
        addPeakEQ(filterName, fc, Q, dBGain=0, coeff=False): Add a biquad peak equalizer filter.
        addLowShelf(filterName, fc, Q, dBGain=0, coeff=False): Add a biquad low-shelf filter.
        addHighShelf(filterName, fc, Q, dBGain=0, coeff=False): Add a biquad high-shelf filter.
        addGain(filterName, dBGain): Add a gain filter.
        addDelay(filterName, dt, timeConvention='-jwt'): Add a time delay filter.
        deleteFilter(filterName): Delete a filter from the crossover.
        plot(divide_stages=True, amplitude=40): Plot the overall frequency response.
    """
    def __init__(self, f_axis):
        self.identifier = "XO"
        self.f_axis = f_axis
        self.h = np.ones(len(f_axis), dtype=complex)
        self.filterList = []
        self.frf_List = []

        # ref 2 study
        self.referenceStudy = None
        self.ref2bem = None

    def addLowPass(self, filterName, order, fc):
        """
        Add an analog low-pass filter to the crossover design.

        Parameters:
            filterName (str): Name of the filter to be added.
            order (int): Order of the low-pass filter (1 or higher).
            fc (float): Cutoff frequency of the filter in Hertz.

        Returns:
            None

        Notes:
            This method adds an analog low-pass filter to the crossover design using the given
            order and cutoff frequency. The filter coefficients are computed using the `butter`
            function from the SciPy library. The frequency response is calculated using the `freqs`
            function.

            After adding the filter, the crossover's `filterList` and `frf_List` attributes are updated,
            and the total frequency response `h` is recalculated using the `updateCrossovers` function.
        """
        self.checkExistingFilters(filterName)

        b, a = butter(order, fc, btype='low', analog=True)
        _, h = freqs(b, a, self.f_axis)

        # update crossover
        self.filterList.append(filterName)
        self.frf_List.append(h)
        self.h = updateCrossovers(self.f_axis, self.filterList, self.frf_List)
        return None

    def addHighPass(self, filterName, order, fc):
        """
        Add an analog high-pass filter to the crossover design.

        Parameters:
            filterName (str): Name of the filter to be added.
            order (int): Order of the high-pass filter (1 or higher).
            fc (float): Cutoff frequency of the filter in Hertz.

        Returns:
            None

        Notes:
            This method adds an analog high-pass filter to the crossover design using the given
            order and cutoff frequency. The filter coefficients are computed using the `butter`
            function from the SciPy library. The frequency response is calculated using the `freqs`
            function.

            After adding the filter, the crossover's `filterList` and `frf_List` attributes are updated,
            and the total frequency response `h` is recalculated using the `updateCrossovers` function.
        """
        self.checkExistingFilters(filterName)

        b, a = butter(order, fc, btype='high', analog=True)
        _, h = freqs(b, a, self.f_axis)

        # update crossover
        self.filterList.append(filterName)
        self.frf_List.append(h)
        self.h = updateCrossovers(self.f_axis, self.filterList, self.frf_List)
        return None

    def addBandPass(self, filterName, order, fc1, fc2):
        """
        Add an analog band-pass filter to the crossover design.

        Parameters:
            filterName (str): Name of the filter to be added.
            order (int): Order of the band-pass filter (1 or higher).
            fc1 (float): Lower cutoff frequency of the filter in Hertz.
            fc2 (float): Upper cutoff frequency of the filter in Hertz.

        Returns:
            None

        Notes:
            This method adds an analog band-pass filter to the crossover design using the given
            order, lower cutoff frequency (fc1), and upper cutoff frequency (fc2). The filter
            coefficients are computed using the `butter` function from the SciPy library. The
            frequency response is calculated using the `freqs` function.

            After adding the filter, the crossover's `filterList` and `frf_List` attributes are updated,
            and the total frequency response `h` is recalculated using the `updateCrossovers` function.

        """
        # analog band-pass filter
        self.checkExistingFilters(filterName)

        b, a = butter(order, [fc1, fc2], btype='bandpass', analog=True)
        _, h = freqs(b, a, self.f_axis)

        # update crossover
        self.filterList.append(filterName)
        self.frf_List.append(h)
        self.h = updateCrossovers(self.f_axis, self.filterList, self.frf_List)
        return None

    def addLowPassBQ(self, filterName, fc, Q, dBGain=0, coeff=False):
        """
        Add a biquad (BQ) low-pass filter to the crossover design.

        Parameters:
            filterName (str): Name of the filter to be added.
            fc (float): Cutoff frequency of the filter in Hertz.
            Q (float): Quality factor of the filter.
            dBGain (float): Gain in decibels applied to the filter (default: 0).
            coeff (bool): If True, return filter coefficients (default: False).

        Returns:
            list or None: If coeff is True, returns a list containing filter coefficients [b, a].
                          Otherwise, returns None.

        Notes:
            This method adds a biquad low-pass filter to the crossover design using the given
            cutoff frequency (fc), quality factor (Q), and optional gain (dBGain).

        """
        self.checkExistingFilters(filterName)

        Fs = self.f_axis[-1] * 4 #2
        wc = 2 * pi * fc / Fs
        alpha = sin(wc) / (2 * Q)

        # thks T
        A = 10 ** (dBGain / 20)
        b = np.zeros(3)
        a = np.zeros(3)
        b[0] = A * (1 - cos(wc)) / 2
        b[1] = A * (1 - cos(wc))
        b[2] = A * (1 - cos(wc)) / 2
        a[0] = 1 + alpha
        a[1] = -2 * cos(wc)
        a[2] = 1 - alpha
        _, h = freqz(b/a[0], a/a[0], worN=self.f_axis, fs=Fs)

        # update crossover
        self.filterList.append(filterName)
        self.frf_List.append(h)
        self.h = updateCrossovers(self.f_axis, self.filterList, self.frf_List)
        # print(b)
        # print(a)

        out = []
        if coeff is True:
            out.append(b)
            out.append(a)
        elif coeff is False:
            out = None
        return out

    def addHighPassBQ(self, filterName, fc, Q, dBGain=0, coeff=False):
        """
        Add a biquad (BQ) high-pass filter to the crossover design.

        Parameters:
            filterName (str): Name of the filter to be added.
            fc (float): Cutoff frequency of the filter in Hertz.
            Q (float): Quality factor of the filter.
            dBGain (float): Gain in decibels applied to the filter (default: 0).
            coeff (bool): If True, return filter coefficients (default: False).

        Returns:
            list or None: If coeff is True, returns a list containing filter coefficients [b, a].
                          Otherwise, returns None.

        Notes:
            This method adds a biquad high-pass filter to the crossover design using the given
            cutoff frequency (fc), quality factor (Q), and optional gain (dBGain).

        """
        self.checkExistingFilters(filterName)

        Fs = self.f_axis[-1] * 4 #2
        wc = 2 * pi * fc / Fs
        alpha = sin(wc) / (2 * Q)

        # thks T
        A = 10 ** (dBGain / 20)
        b = np.zeros(3)
        a = np.zeros(3)
        b[0] = A * (1 + cos(wc)) / 2
        b[1] = -A * (1 + cos(wc))
        b[2] = A * (1 + cos(wc)) / 2
        a[0] = 1 + alpha
        a[1] = -2 * cos(wc)
        a[2] = 1 - alpha
        _, h = freqz(b/a[0], a/a[0], worN=self.f_axis, fs=Fs)

        # update crossovers
        self.filterList.append(filterName)
        self.frf_List.append(h)
        self.h = updateCrossovers(self.f_axis, self.filterList, self.frf_List)

        out = []
        if coeff is True:
            out.append(b)
            out.append(a)
        elif coeff is False:
            out = None
        return out

    def addBandPassBQ(self, filterName, fc, Q, dBGain=0, coeff=False):
        """
        Add a biquad (BQ) band-pass filter to the crossover design.

        Parameters:
            filterName (str): Name of the filter to be added.
            fc (float): Cutoff frequency of the filter in Hertz.
            Q (float): Quality factor of the filter.
            dBGain (float): Gain in decibels applied to the filter (default: 0).
            coeff (bool): If True, return filter coefficients (default: False).

        Returns:
            list or None: If coeff is True, returns a list containing filter coefficients [b, a].
                          Otherwise, returns None.

        Notes:
            This method adds a biquad band-pass filter to the crossover design using the given
            cutoff frequency (fc), quality factor (Q), and optional gain (dBGain).

        """
        self.checkExistingFilters(filterName)

        Fs = self.f_axis[-1] * 4 #2
        wc = 2 * pi * fc / Fs
        alpha = sin(wc) / (2 * Q)

        # thks T
        A = 10 ** (dBGain / 20)
        b = np.zeros(3)
        a = np.zeros(3)
        b[0] = A*alpha
        b[1] = 0
        b[2] = -A*alpha
        a[0] = 1 + alpha
        a[1] = -2*cos(wc)
        a[2] = 1 - alpha
        _, h = freqz(b/a[0], a/a[0], worN=self.f_axis, fs=Fs)

        # update crossovers
        self.filterList.append(filterName)
        self.frf_List.append(h)
        self.h = updateCrossovers(self.f_axis, self.filterList, self.frf_List)

        out = []
        if coeff is True:
            out.append(b)
            out.append(a)
        elif coeff is False:
            out = None
        return out

    def addPeakEQ(self, filterName, fc, Q, dBGain=0, coeff=False):
        """
        Add a biquad (BQ) peak-eq filter to the crossover design.

        Parameters:
            filterName (str): Name of the filter to be added.
            fc (float): Cutoff frequency of the filter in Hertz.
            Q (float): Quality factor of the filter.
            dBGain (float): Gain in decibels applied to the filter (default: 0).
            coeff (bool): If True, return filter coefficients (default: False).

        Returns:
            list or None: If coeff is True, returns a list containing filter coefficients [b, a].
                          Otherwise, returns None.

        Notes:
            This method adds a biquad peak-eq filter to the crossover design using the given
            cutoff frequency (fc), quality factor (Q), and optional gain (dBGain).

        """
        self.checkExistingFilters(filterName)

        Fs = self.f_axis[-1] * 4 #2
        wc = 2 * pi * fc / Fs
        alpha = sin(wc) / (2 * Q)

        # thks T
        A = 10 ** (dBGain / 40)
        b = np.zeros(3)
        a = np.zeros(3)
        b[0] = 1 + alpha * A
        b[1] = -2 * cos(wc)
        b[2] = 1 - alpha * A
        a[0] = 1 + alpha / A
        a[1] = -2 * cos(wc)
        a[2] = 1 - alpha / A
        _, h = freqz(b/a[0], a/a[0], worN=self.f_axis, fs=Fs)

        # update crossovers
        self.filterList.append(filterName)
        self.frf_List.append(h)
        self.h = updateCrossovers(self.f_axis, self.filterList, self.frf_List)

        out = []
        if coeff is True:
            out.append(b)
            out.append(a)
        elif coeff is False:
            out = None
        return out

    def addLowShelf(self, filterName, fc, Q, dBGain=0, coeff=False):
        """
        Add a biquad (BQ) low-shelf filter to the crossover design.

        Parameters:
            filterName (str): Name of the filter to be added.
            fc (float): Cutoff frequency of the filter in Hertz.
            Q (float): Quality factor of the filter.
            dBGain (float): Gain in decibels applied to the filter (default: 0).
            coeff (bool): If True, return filter coefficients (default: False).

        Returns:
            list or None: If coeff is True, returns a list containing filter coefficients [b, a].
                          Otherwise, returns None.

        Notes:
            This method adds a biquad low-shelf filter to the crossover design using the given
            cutoff frequency (fc), quality factor (Q), and optional gain (dBGain).

        """
        self.checkExistingFilters(filterName)

        Fs = self.f_axis[-1] * 4 #2
        wc = 2 * pi * fc / Fs
        alpha = sin(wc) / (2 * Q)

        # thks T
        A = 10 ** (dBGain / 40)
        b = np.zeros(3)
        a = np.zeros(3)
        b[0] = A*((A+1) - (A-1)*cos(wc) + 2*sqrt(A)*alpha)
        b[1] = 2*A*((A-1) - (A+1)*cos(wc))
        b[2] = A*((A+1) - (A-1)*cos(wc) - 2*sqrt(A)*alpha)
        a[0] = (A+1) + (A-1)*cos(wc) + 2*sqrt(A)*alpha
        a[1] = -2*((A-1) + (A+1)*cos(wc))
        a[2] = (A+1) + (A-1)*cos(wc) - 2*sqrt(A)*alpha
        _, h = freqz(b/a[0], a/a[0], worN=self.f_axis, fs=Fs)

        # update crossovers
        self.filterList.append(filterName)
        self.frf_List.append(h)
        self.h = updateCrossovers(self.f_axis, self.filterList, self.frf_List)

        out = []
        if coeff is True:
            out.append(b)
            out.append(a)
        elif coeff is False:
            out = None
        return out

    def addHighShelf(self, filterName, fc, Q, dBGain=0, coeff=False):
        """
        Add a biquad (BQ) high-shelf filter to the crossover design.

        Parameters:
            filterName (str): Name of the filter to be added.
            fc (float): Cutoff frequency of the filter in Hertz.
            Q (float): Quality factor of the filter.
            dBGain (float): Gain in decibels applied to the filter (default: 0).
            coeff (bool): If True, return filter coefficients (default: False).

        Returns:
            list or None: If coeff is True, returns a list containing filter coefficients [b, a].
                          Otherwise, returns None.

        Notes:
            This method adds a biquad high-shelf filter to the crossover design using the given
            cutoff frequency (fc), quality factor (Q), and optional gain (dBGain).

        """
        self.checkExistingFilters(filterName)

        Fs = self.f_axis[-1] * 4 #2
        wc = 2 * pi * fc / Fs
        alpha = sin(wc) / (2 * Q)

        # thks T
        A = 10 ** (dBGain / 40)
        b = np.zeros(3)
        a = np.zeros(3)
        b[0] = A*( (A+1) + (A-1)*cos(wc) + 2*sqrt(A)*alpha )
        b[1] = -2*A*( (A-1) + (A+1)*cos(wc) )
        b[2] = A*( (A+1) + (A-1)*cos(wc) - 2*sqrt(A)*alpha )
        a[0] = (A+1) - (A-1)*cos(wc) + 2*sqrt(A)*alpha
        a[1] = 2*( (A-1) - (A+1)*cos(wc) )
        a[2] = (A+1) - (A-1)*cos(wc) - 2*sqrt(A)*alpha
        _, h = freqz(b/a[0], a/a[0], worN=self.f_axis, fs=Fs)

        # update crossovers
        self.filterList.append(filterName)
        self.frf_List.append(h)
        self.h = updateCrossovers(self.f_axis, self.filterList, self.frf_List)

        out = []
        if coeff is True:
            out.append(b)
            out.append(a)
        elif coeff is False:
            out = None
        return out

    def addGain(self, filterName, dBGain):
        """
        Add a dB gain to the crossover design.

        Parameters:
            filterName (str): Name of the filter to be added.
            dBGain (float): Gain in decibels applied to the filter (default: 0).

        Returns:
            list or None: If coeff is True, returns a list containing filter coefficients [b, a].
                          Otherwise, returns None.

        Notes:
            This method adds a general gain to the crossover design using the given
            the gain (dBGain).

        """
        self.checkExistingFilters(filterName)

        h = np.ones(len(self.f_axis), dtype=complex) * 10**(dBGain/20)

        # update crossovers
        self.filterList.append(filterName)
        self.frf_List.append(h)
        self.h = updateCrossovers(self.f_axis, self.filterList, self.frf_List)
        return None

    def addTransferFunction(self, filterName, tf):
        self.checkExistingFilters(filterName)

        self.filterList.append(filterName)
        self.frf_List.append(tf)
        self.h = updateCrossovers(self.f_axis, self.filterList, self.frf_List)
        return None

    def addDelay(self, filterName, dt, timeConvention='-jwt'):
        """
        Add a time delay filter to the crossover design.

        Parameters:
            filterName (str): Name of the filter to be added.
            dt (float): time delay in seconds.
            timeConvention (str): time convention to be used (-jwt; +jwt)

        Returns:
            list or None: If coeff is True, returns a list containing filter coefficients [b, a].
                          Otherwise, returns None.

        Notes:
            This method adds a time delay to the crossover design using the given
            delay time (dt) and its time convention (timeConvention).

        """
        self.checkExistingFilters(filterName)

        if timeConvention == '-jwt' or timeConvention == '-':
            tc = -1
        elif timeConvention == '+jwt' or timeConvention == 'jwt' or timeConvention == '+':
            tc = 1
        om = 2*pi*self.f_axis
        h = np.exp(tc * 1j*om*dt)

        # update crossovers
        self.filterList.append(filterName)
        self.frf_List.append(h)
        self.h = updateCrossovers(self.f_axis, self.filterList, self.frf_List)
        return None

    def addPhaseFlip(self, filterName):
        self.checkExistingFilters(filterName)
        h = np.ones(len(self.f_axis), dtype=complex) * -1

        self.filterList.append(filterName)
        self.frf_List.append(h)
        self.h = updateCrossovers(self.f_axis, self.filterList, self.frf_List)
        return None

    def deleteFilter(self, filterName):
        # delete corresponding filter stage from crossover
        if filterName in self.filterList:
            index = self.filterList.index(filterName)
            self.filterList.pop(index)
            self.frf_List.pop(index)

            # update crossovers
        self.h = updateCrossovers(self.f_axis, self.filterList, self.frf_List)
        return None

    def checkExistingFilters(self, filterName):
        if filterName in self.filterList:
            self.deleteFilter(filterName)
        return None

    def plot(self, divide_stages=True, amplitude=64):
        # plot total response

        # get max value for ylim
        maxList = np.zeros(len(self.filterList) + 1)
        for i in range(len(self.filterList)):
            maxList[i] = np.max(gtb.gain.dB(self.frf_List[i])) + 3
        maxList[-1] = np.max(gtb.gain.dB(self.h)) + 3
        max_dB = +amplitude / 2#np.max(maxList)
        min_dB = -amplitude / 2 #max_dB - amplitude

        fig, ax = plt.subplots(2, 1)
        if divide_stages is True:
            for i in range(len(self.filterList)):
                ax[0].semilogx(self.f_axis, gtb.gain.dB(self.frf_List[i]), label=self.filterList[i])
                ax[1].semilogx(self.f_axis, np.angle(self.frf_List[i]), '--', label=self.filterList[i])
            ax[0].semilogx(self.f_axis, gtb.gain.dB(self.h), 'k', label='H(f)')
            ax[1].semilogx(self.f_axis, np.angle(self.h), '--k', label='H(f)')


        elif divide_stages is False:
            ax[0].semilogx(self.f_axis, gtb.gain.dB(self.h), 'k', label='H(f)')
            ax[1].semilogx(self.f_axis, np.angle(self.h), '--k', label='H(f)')

        else:
            Exception('divide_stages must be either True or False.')

        for i in range(2):
            ax[i].legend()
            ax[i].grid(which='both')
        ax[0].set(ylabel='Gain [dB]', ylim=[min_dB, max_dB])
        ax[1].set(xlabel='Frequency [Hz]', ylabel='Phase [rad]')
        plt.tight_layout()
        plt.show()
        return None

    def save(self, fileName):
        np.savez(fileName,
                 identifier = self.identifier,
                 f_axis     = self.f_axis,
                 h          = self.h,
                 filterList = self.filterList,
                 frf_list   = self.frf_List,
                 )
        return None


# update crossover class
def updateCrossovers(f_axis, filterList, frfList):
    # update frequency response of crossover
    nFilters = len(filterList)
    h = np.ones(len(f_axis), dtype=complex)
    for i in range(nFilters):
        h *= frfList[i]
    return h


# ACOUSTIC RESPONSE EQUALIZER
class response_eq:
    def __init__(self, freq, H, bounds):
        self.freq = freq
        self.H = np.abs(H)
        self.bounds = bounds
        self.fmin = bounds[0]
        self.fmax = bounds[1]
        self.ifmin, _ = gtb.findInArray(freq, self.fmin)
        self.ifmax, _ = gtb.findInArray(freq, self.fmax)
        self.target = np.ones(len(freq))
        self.xover = xovers(freq)

        # initialize some values
        self.intersec = None
        self.peaks = None


        # run analyser
        self.analyse_response()

    def analyse_response(self):
        self.intersec = self.find_intersections()
        self.fLogMean = self.find_logMean()
        self.peaks = self.find_peaks()
        for i in range(len(self.intersec)):
            try:
                p_opt, _ = sp.optimize.curve_fit(filter_optimize, self.freq[self.intersec[i]:self.intersec[i+1]],
                                                 - gtb.gain.dB(self.H[self.intersec[i]:self.intersec[i+1]]),
                                                 p0=[self.fLogMean[i], 1.5, gtb.gain.dB(1 + (1-self.peaks[i]))],
                                                 maxfev=10000)
                fc = p_opt[0]
                q = p_opt[1]
                g = p_opt[2]
                self.xover.addPeakEQ("peq_"+str(i), fc, q, g)
                self.H *= self.xover.frf_List[i]
            except:
                None
        return None

    def find_intersections(self):
        target_function = self.target
        transfer_function = self.H
        intersections = []
        for i in np.arange(self.ifmin, self.ifmax, 1): #range(len(target_function) - 1):
            if (target_function[i] - transfer_function[i]) * (target_function[i+1] - transfer_function[i+1]) < 0:
                intersections.append(i)
        return intersections

    def find_peaks(self):
        peaks = np.zeros(len(self.intersec) - 1)
        for i in range(len(peaks)):
            if np.mean(self.H[self.intersec[i]:self.intersec[i+1]]) < 1:
                peaks[i] = np.min(np.abs(self.H[self.intersec[i]:self.intersec[i+1]]))
            elif np.mean(self.H[self.intersec[i]:self.intersec[i+1]]) > 1:
                peaks[i] = np.max(np.abs(self.H[self.intersec[i]:self.intersec[i+1]]))
        return peaks

    def find_logMean(self):
        fLogMean = np.zeros(len(self.intersec) - 1)
        for i in range(len(fLogMean)):
            fLogMean[i] = np.sqrt(self.freq[self.intersec[i]] * self.freq[self.intersec[i+1]])
        return fLogMean


def filter_optimize(freq, fc, q, g):
    """
    Build a simple peaking EQ to fit data using scipy.optimize.curve_fit()
    """
    xo = xovers(freq)
    xo.addPeakEQ('peq', fc, q, g)
    return gtb.gain.dB(xo.h)
