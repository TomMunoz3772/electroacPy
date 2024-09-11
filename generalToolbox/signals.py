# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 16:07:10 2022

@author: tom
"""
import numpy as np

def SyncSineSweep(f1, f2, T, Fs, weight=0.8):
    """

    Parameters
    ----------
    f1 : start frequency.
    f2 : stop frequency.
    T : time length of sine sweep.
    Fs : Frequency sampling

    Returns
    -------
    t : time vector
    sss : synchronized sine sweep.

    """
    Ts = 1 / Fs
    t = np.arange(0, T, Ts)
    L = 1 / f1 * np.round(T * f1 / np.log(f2 / f1))
    syncSine = weight * np.sin(2 * np.pi * f1 * L * np.exp(t / L))
    return t, syncSine
