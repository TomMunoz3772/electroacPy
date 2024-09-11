from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

import numpy as np
import generalToolbox as gtb


class MplDirectivityCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=8.1, height=6, title=""):
        self.fig = Figure(figsize=(width, height))
        self.fig.suptitle(title)
        self.ax_SPL_directivity = self.fig.add_subplot(221)
        self.ax_SPL_response = self.fig.add_subplot(222)
        self.ax_dB_directivity = self.fig.add_subplot(223)
        self.ax_polar_response = self.fig.add_subplot(224, polar=True)
        super(MplDirectivityCanvas, self).__init__(self.fig)
        self.setParent(parent)

    def initialize_plot(self, frequency, angle, pMic2Plot):
        # SPL DIRECTIVITY
        self.ax_SPL_directivity.contourf(frequency, angle,
                                         gtb.gain.SPL(pMic2Plot).T, cmap="turbo")

        # Normalized directivity
        maxAngle = np.array([np.max(abs(pMic2Plot), 1)])
        maxMatrix = np.repeat(maxAngle, len(angle), 0).T
        directivity = gtb.gain.dB(np.abs(pMic2Plot) / maxMatrix)
        self.ax_dB_directivity.contourf(frequency, angle,
                                        directivity.T, levels=np.arange(-21, 3, 3), cmap="turbo")


        # SPL PLOT
        idx_0, val_0 = gtb.findInArray(angle, 0)
        self.ax_SPL_response.semilogx(frequency, gtb.gain.SPL(pMic2Plot[:, idx_0]))

        # polar
        ind_freq, val_f = gtb.findInArray(frequency, 1000)
        self.ax_polar_response.plot(np.deg2rad(angle), gtb.gain.SPL(pMic2Plot[ind_freq, :]),
                                    label=int(val_f))
        self.ax_polar_response.legend(loc='best')

        # rescale properly
        self.fig.tight_layout()
        self.draw()
        return None

    def update_plot(self, frequency, angle, pMic2Plot, parameters):
        self.ax_dB_directivity.clear()
        self.ax_SPL_directivity.clear()
        self.ax_SPL_response.clear()
        self.ax_polar_response.clear()

        # LOAD PARAMETERS
        if "dBmin" in parameters:
            dBmin = parameters["dBmin"]
            dBmax = parameters["dBmax"]
            dBstep = parameters["dBstep"]
            LEVELS = np.arange(dBmin, dBmax + dBstep, dBstep)

            colormap = parameters["cmap"]
            xscale = parameters["xscale"]
            fmin = parameters["fmin"]
            fmax = parameters["fmax"]
            freq2plot = parameters["freq2plot"]
            angle2plot = parameters["angle2plot"]
        else:
            dBmax = np.max(gtb.gain.SPL(pMic2Plot)) + 6
            dBmin = dBmax - 46
            dBstep = 6
            LEVELS = np.arange(dBmin, dBmax + dBstep, dBstep)
            colormap = "turbo"
            xscale = "log"
            fmin = frequency[0]
            fmax = frequency[-1]
            angle2plot = []
            freq2plot = []

        # PLOT
        # CONTOUR SPL
        self.ax_SPL_directivity.contourf(frequency, angle,
                                         gtb.gain.SPL(pMic2Plot).T, cmap=colormap, levels=LEVELS)
        self.ax_SPL_directivity.set(xlabel="Frequency [Hz]", ylabel="Angle", xlim=[fmin, fmax])
        self.ax_SPL_directivity.set_xscale(xscale)

        # TODO -> CHANGE TO A NORMALIZED DIRECTIVITY
        maxAngle = np.array([np.max(abs(pMic2Plot), 1)])
        maxMatrix = np.repeat(maxAngle, len(angle), 0).T
        directivity = gtb.gain.dB(np.abs(pMic2Plot) / maxMatrix)
        self.ax_dB_directivity.contourf(frequency, angle,
                                        directivity.T, levels=np.arange(-21, 3, 3), cmap=colormap)
        self.ax_dB_directivity.set(xlabel="Frequency [Hz]", ylabel="Angle", xlim=[fmin, fmax])
        self.ax_dB_directivity.set_xscale(xscale)

        # FREQUENCY RESPONSE (1D)
        if bool(angle2plot) is True:
            for i in range(len(angle2plot)):
                idx, val = gtb.findInArray(angle, angle2plot[i])
                self.ax_SPL_response.plot(frequency, gtb.gain.SPL(pMic2Plot[:, idx]).T, label="{}".format(val))
        else:
            idx0, _ = gtb.findInArray(angle, 0)
            self.ax_SPL_response.plot(frequency, gtb.gain.SPL(pMic2Plot[:, idx0]).T, label="on-axis")
        self.ax_SPL_response.set(xlabel="Frequency [Hz]", ylabel="SPL [dB]", xlim=[fmin, fmax],
                                 ylim=[dBmin, dBmax])
        self.ax_SPL_response.set_xscale(xscale)
        self.ax_SPL_response.grid(linestyle='dotted', which="both")
        self.ax_SPL_response.legend(loc="best")


        # POLAR
        if bool(freq2plot) is True:
            for i in range(len(freq2plot)):
                ind_freq, val_f = gtb.findInArray(frequency, freq2plot[i])
                self.ax_polar_response.plot(np.deg2rad(angle),
                                            gtb.gain.SPL(pMic2Plot[ind_freq, :]), label=int(val_f))
        else:
            ind_freq, val_f = gtb.findInArray(frequency, 1000)
            self.ax_polar_response.plot(np.deg2rad(angle),
                                        gtb.gain.SPL(pMic2Plot[ind_freq, :]), label=int(val_f))
        self.ax_polar_response.legend(loc='best')

        # RESCALE PROPERLY
        self.fig.tight_layout()
        self.draw()
        return None
