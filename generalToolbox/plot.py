import numpy as np
import generalToolbox as gtb
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, Slider, TextBox, Button
from matplotlib.backend_bases import MouseEvent, MouseButton
from matplotlib.contour import QuadContourSet
from matplotlib import cm, colors
import vtk
import pyvista
import os

def plotTotalDirectivity(theta, freq, pMic, xscale='log', fmin=20, fmax=20e3,
                         dBmin=-20, dBmax=False, title="title"):
    """
    Plot the directivity from pMic obtained with getMicPressure for a circular microphone array

    Parameters
    ----------
    theta : narray
        observation angles.
    freq : narray
        frequency array.
    pMic : narray
        Pressur obtained with getMicPressure().
    xscale : str, optional
        scaling of the x axis (logarithmic, linear). The default is 'log'.
    fmin : float, optional
        min xlim. The default is 20.
    fmax : float, optional
        max xlim. The default is 20e3.
    dBmin : float, optional
        min ylim. The default is -20.
    dBmax : float, optional
        max ylim. The default is False.
    norm : bool, optional
        set the normalization. The default is True.

    Returns
    -------
    None.

    """
    # compute directivity
    maxAngle = np.array([np.max(abs(pMic), 1)])
    maxMatrix = np.repeat(maxAngle, len(theta), 0).T
    directivity = gtb.gain.dB(np.abs(pMic) / maxMatrix)
    dBmax = int(np.max(gtb.gain.SPL(pMic))) + 3
    dBmin = dBmax - 40

    fig = plt.figure(figsize=(11, 8))
    fig.suptitle(title)
    fig.subplots_adjust(top=0.88)
    ax1 = fig.add_subplot(221)
    gca1 = ax1.contourf(freq, np.rad2deg(theta), gtb.gain.SPL(pMic).T,
                        np.arange(dBmin, dBmax + 3, 3), cmap='turbo')
    ax1.set_xscale('log')
    ax1.set(xlabel="Frequency [Hz]", ylabel="Angle [deg]", xlim=[fmin, fmax],
            title='Raw directivity (SPL)')
    plt.colorbar(gca1)

    # Total directivity NORMALIZED
    ax2 = fig.add_subplot(223)
    gca2 = ax2.contourf(freq, np.rad2deg(theta), directivity.T, np.arange(-21, 3, 3), cmap='turbo')
    ax2.set_xscale('log')
    ax2.set(xlabel="Frequency [Hz]", ylabel="Angle [deg]", xlim=[fmin, fmax],
            title='Absolute directivity')
    plt.colorbar(gca2)

    # pressure response
    obsAngle = 0
    ind_angle, value_angle = gtb.findInArray(theta, np.deg2rad(obsAngle))
    ax3 = fig.add_subplot(222)
    ax3.semilogx(freq, gtb.gain.SPL(pMic[:, ind_angle]), label=int(np.rad2deg(value_angle)))
    ax3.set(xlabel='Frequency [Hz]', ylabel='SPL [dB]', ylim=[dBmin, dBmax],
            title='Pressure response')
    ax3.grid(which='both')
    ax3.legend(loc='best')

    # polar
    freqObs = 1000
    ind_freq, value_freq = gtb.findInArray(freq, freqObs)
    ax4 = fig.add_subplot(224, polar=True)
    ax4.plot(theta, gtb.gain.SPL(pMic[ind_freq, :]), label=int(value_freq))
    ax4.set(ylim=[dBmin, dBmax], title='Polar directivity')
    ax4.legend(loc='best')

    # defining the cursor
    cursor = Cursor(ax1, horizOn=True, vertOn=True, color='black', linewidth=1.2,
                    useblit=True)

    def onclick(event):
        # global coord
        # coord.append((event.xdata, event.ydata))
        if event.button is MouseButton.LEFT:
            x = event.xdata
            y = event.ydata

            # update subplots -- FRF
            obsAngle = y
            ind_angle, value_angle = gtb.findInArray(theta, np.deg2rad(obsAngle))
            ax3.semilogx(freq, gtb.gain.SPL(pMic[:, ind_angle]), label=int(np.rad2deg(value_angle)))
            ax3.legend(loc='best')

            # update subplots -- polar
            freqObs = x
            ind_freq, value_freq = gtb.findInArray(freq, freqObs)
            ax4.plot(theta, gtb.gain.SPL(pMic[ind_freq, :]), label=int(value_freq))
            ax4.legend(loc='best')

            fig.canvas.draw()  # redraw the figure

        if event.button is MouseButton.RIGHT:
            if len(ax3.lines) > 1:
                ax3.lines[-1].remove()
                ax4.lines[-1].remove()
                ax3.legend(loc='best')
                ax4.legend(loc='best')
                fig.canvas.draw()

    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.tight_layout()
    plt.show()
    return cursor