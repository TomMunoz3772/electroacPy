#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:50:43 2023

@author: tom.munoz
"""
import numpy as np
from electroacPy.globalVariables import air
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
import generalToolbox as gtb
import generalToolbox.lp_loaders as lpl
import os
from copy import copy

pi = np.pi

class electroAcousticDriver:
    def __init__(self, U, Le, Re, Cms, Mms, Rms, Bl, Sd, f_array,
                 c=air.c, rho=air.rho):
        """
        Create an electro-acoustic driver from its Thiele-Small parameters.

        Parameters
        ----------
        U : float
            Input voltage (usually 1).
        Le : float
            Coil inductance (H).
        Re : float
            Coil resistance (Ohms).
        Cms : float
            Suspension's compliance (m/N).
        Mms : float
            Moving mass of the driver (kg).
        Rms : float
            Mechanical resistance (N.s/m).
        Bl : float
            Force factor (N.A / T.m).
        Sd : float
            Diaphragm surface area (m^2).
        f_array : array-like
            Frequencies at which the solutions are computed.
        c : float, optional
            Speed of sound. The default is air.c.
        rho : float, optional
            Density of air. The default is air.rho.

        Returns
        -------
        None

        Notes
        -----
        This class represents an electro-acoustic driver and calculates its electrical and mechanical impedance,
        equivalent acoustical impedance, input pressure, radiation impedance, quality factors, and other parameters.

        """
        # identifier
        self.identifier = "EAC"

        # medium properties
        w = 2*pi*f_array
        Zc = rho * c
        k = w / c
        
        # self param
        self.U = U
        self.Le = Le
        self.Re = Re
        self.Cms = Cms
        self.Mms = Mms
        self.Rms = Rms
        self.Bl = Bl
        self.Sd = Sd
        self.f_array = f_array

        # medium returns
        self.c = c
        self.k = k
        self.f = f_array
        self.w = w
        self.rho = rho

        # speaker radius
        r = np.sqrt(Sd / pi)
        self.r = r
        self.Sd = Sd

        # impedance
        s = 1j*w
        self.Ze = Re + s*Le
        Zms = Rms + s*Mms + 1/(s*Cms)

        # equivalent acoustical impedance
        self.Zac = (1/Sd**2)*(Bl**2/self.Ze)
        self.Zas = Zms / Sd**2
        self.Zms = Zms
        self.ZeTot = self.Ze + Bl**2/Zms
        self.Bl = Bl
        # equivalent input pressure
        self.Ps = U * Bl / (self.Ze * Sd)

        # Radiation impedance (speaker front impedance)
        Mrad = 8 * rho * r / 3 / pi / Sd
        Rrad = Zc / Sd * (k*r)**2 / 2
        self.Zaf = Rrad + 1j*w*Mrad
        self.Zs = self.Zac + self.Zas + self.Zaf # total acoustical impedance coming from the driver

        # Quality factors and others
        self.Qes = Re / (Bl)**2 * np.sqrt(Mms/Cms)
        self.Qms = 1/Rms * np.sqrt(Mms/Cms)
        self.Qts = self.Qms*self.Qes / (self.Qms + self.Qes)
        self.Vas = rho*c**2*Sd**2*Cms
        self.Fs = 1 / (2*pi*np.sqrt(Cms*Mms))
        self.EBP = self.Fs/self.Qes


        # Ref Signals #
        # Velocity
        # av = [Mms*Le, Le*Rms + Mms*Re, Rms*Re + Le / Cms + Bl**2, Re/Cms]
        # bv = [Bl, 0]
        # _, self.Hv = freqs(bv, av, worN=f_array)
        self.Hv = Bl/self.Ze  / (self.Zms + Bl**2 / self.Ze)

        # Displacement
        # ax = av
        # bx = [Bl]
        # _, self.Hx = freqs(bx, ax, worN=f_array)
        self.Hx = self.Hv / s

        # Acceleration
        # aa = av
        # ba = [Bl, 0, 0]
        # _, self.Ha = freqs(ba, aa, worN=f_array)
        self.Ha = self.Hv * s

        # Acoustic simulation reference
        self.ref2bem = None

        # in box velocity and impedance
        self.inBox    = False
        self.isPorted = False  # easier to manage if radiator in study_ is not a speakerBox object
        self.v        = self.Hv

        # references
        self.ref2bem   = False
        self.poly_data = False  # is class from polytech?

        
    def plotZe(self):
        """
        Plot the electrical impedance ZeTot in both modulus and phase.

        Returns
        -------
        None

        """
        fig, ax = plt.subplots()
        ax.semilogx(self.f_array, np.abs(self.ZeTot), 
                    label="modulus")
        ax.legend(loc='upper left')
        ax.grid(linestyle='dotted', which='both')
        ax.set(ylabel="|Z| [Ohm]", xlabel="Frequency [Hz]")
        ax2 = ax.twinx()
        ax2.semilogx(self.f_array, np.angle(self.ZeTot), '--', 
                     label="phase")
        ax2.legend(loc='upper right')
        ax2.set(ylabel="phase [rads]", ylim=[-pi, pi])
        ax2.grid()
        plt.tight_layout()
        return plt.show()
    
    def plotXVA(self):
        """
        Plot the displacement, velocity, and acceleration frequency responses.

        Returns
        -------
        None

        """
        fig, ax = plt.subplots(3, 1)
        ax[0].semilogx(self.f_array, np.abs(self.Hx*1e3), label='Displacement')
        ax[1].semilogx(self.f_array, np.abs(self.Hv), label='Velocity')
        ax[2].semilogx(self.f_array, np.abs(self.Ha), label='Acceleration')
        ax[2].set(xlabel="Frequency [Hz]")
        ax[0].set(ylabel="mm")
        ax[1].set(ylabel="m/s")
        ax[2].set(ylabel="m/s^2")
        for i in range(3):
            ax[i].grid(which='both')
            ax[i].legend(loc='best')
        plt.tight_layout()
        return plt.show()

    def getThieleSmallParam(self):
        """
        Print out the Thiele/Small parameters of the electro-acoustic driver.

        Returns
        -------
        None

        """
        greetingStr = "Thiele/Small parameters"
        print(greetingStr)
        print("-"*len(greetingStr))
        print("--- Electrical ---")
        print("Re = ", self.Re, " Ohm")
        print("Le = ", self.Le*1e3, " mH")
        print("Bl = ", self.Bl, " N/A")
        
        print("--- Mechanical ---")
        print("Rms = ", round(self.Rms, 2), " N.s/m")
        print("Mms = ", round(self.Mms*1e3, 2), " g")
        print("Cms = ", round(self.Cms*1e3, 2), " mm/N")
        print("Kms = ", round(1/self.Cms), "N/m")
        print("Sd = ", round(self.Sd*1e4, 2), " cm^2")
    
        print("--- Quality Factors ---")
        print("Qes = ", round(self.Qes, 2))
        print("Qms = ", round(self.Qms, 2))
        print("Qts = ", round(self.Qts, 2))
        
        print("--- Others ---")
        print("Fs = ", round(self.Fs, 2), " Hz")
        print("Vas = ", self.Vas, " m^3")
        return None

    def sealedAlignment(self):
        """
        Compute Volume from Qtc value

        Parameters
        ----------
        driver : class
            electro_acoustic_driver object.
        Qtc : total quality factor (mechanical, electrical, acoustical)
        c : speed of sound. The default is air.c.
        rho : air density. The default is air.rho.

        Returns
        -------
        Vb : sealed enclosure volume.
        fc : resonance frequency of the driver inside the enclosure (without radiation mass)
        """
        driver = self
        c = self.c
        rho = self.rho

        ## box parameters
        Vb = driver.Vas
        fc = driver.Fs * np.sqrt(driver.Vas / Vb + 1)
        Qtc = fc / driver.Fs * driver.Qts
        Cab = Vb / rho / c**2

        ## radiated pressure at 1 m
        f_axis = driver.f_array
        omega = 2*np.pi*f_axis
        k = omega/c

        Zac = driver.Zac
        Zas = driver.Zas
        Zaf = driver.Zaf
        Zab = 1/1j/omega/Cab
        Ps = driver.Ps
        Qs = Ps / (Zac + Zas + Zab) # removed Zaf

        p = 1j*k*rho*c*Qs * np.exp(-1j*k*1) / (2 * np.pi * 1)
        Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * (Zab))# + driver.Zaf))

        ## data plot
        fig, ax = plt.subplots(2, 1)
        ax[0].semilogx(f_axis, gtb.gain.SPL(p), label='SPL')
        ax[1].semilogx(f_axis, np.abs(Ze), label='Magnitude')

        ax[1].set(xlabel='Frequency [Hz]', ylabel='impedance')
        ax[0].set(ylabel='SPL [dB] at 1 meter')
        ax[0].legend(loc='best')
        ax[1].legend(loc='best')
        for i in range(2):
            ax[i].grid(which='both')

        plt.subplots_adjust(bottom=0.25)

        # creation of text boxes
        defaultQtc = round(Qtc, 4)
        graphBox = fig.add_axes([0.2, 0.05, 0.1, 0.075])
        bQtc = TextBox(graphBox, "Qtc: ")

        volumeBox = fig.add_axes([0.4, 0.05, 0.1, 0.075])
        bVb = TextBox(volumeBox, "Vb (L): ")

        def update(expr):
            # update Cab value and [p, Ze]
            Qtc = float(expr)
            fc = Qtc / driver.Qts * driver.Fs
            Vb = driver.Vas / ((fc/driver.Fs)**2 - 1)
            Cab = Vb / rho / c ** 2
            Zab = 1 / 1j / omega / Cab

            Qs = Ps / (Zac + Zas + Zab)  # remove Zaf
            p = 1j * k * rho * c * Qs * np.exp(1j * k * 1) / (2 * np.pi * 1)
            Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * (Zab)) #+ driver.Zaf))

            ax[0].clear()
            ax[1].clear()

            ax[0].semilogx(f_axis, gtb.gain.SPL(p), label='SPL')
            ax[1].semilogx(f_axis, np.abs(Ze), label='Magnitude')
            ax[1].set(xlabel='Frequency [Hz]', ylabel='Impedance [Ohm]')
            ax[0].set(ylabel='SPL [dB] at 1 meter')
            ax[0].legend(loc='best')
            ax[1].legend(loc='upper left')
            for i in range(2):
                ax[i].grid(which='both')
            plt.draw()
            bVb.set_val(str(round(Vb*1e3, 4)))

        def updateVolume(expr):
            Vbexpr = float(expr)*1e-3
            fc = driver.Fs * np.sqrt(driver.Vas / Vbexpr + 1)
            Qtc = fc/driver.Fs * driver.Qts
            bQtc.set_val(str(round(Qtc, 4)))

        # Qtc box
        bQtc.on_submit(update)
        bQtc.set_val(str(defaultQtc))  # set default value

        # volume box (input in L)
        bVb.on_submit(updateVolume)
        bVb.set_val(str(round(Vb*1e3, 4)))
        return bQtc, bVb, plt.show()

    def portedAlignment(self):
        """
        Adjust ported box alignment and plot results.

        Returns:
        -------
        None
        """

        driver = copy(self)
        c = self.c
        rho = self.rho

        ## load data
        f_axis = driver.f_array
        omega = 2 * np.pi * f_axis
        k = omega / c
        s = 1j * omega

        ## default box parameters
        self.Vb = copy(driver.Vas)

        fc = driver.Fs * np.sqrt(driver.Vas / self.Vb + 1)
        Qtc = fc / driver.Fs * driver.Qts

        eta = 1e-5

        # ports dimensions
        self.Lp = 343 / driver.Fs / 100  # default and arbitrary length
        self.rp = self.Lp / 2  # default and arbitrary radius
        self.Sp = np.pi * self.rp ** 2

        # box impedance
        self.Cab = self.Vb / rho / c ** 2  # compliance of the enclosure
        self.Rab = rho * c / eta / self.Vb
        self.Zbox = gtb.parallel(1 / s / self.Cab, self.Rab)

        # port impedance
        self.Map = rho * self.Lp / self.Sp  # acoustical mass of the port
        self.Mal = 0.85 * 2 * self.rp  # length correction
        self.Map += self.Mal

        # radiation impedance (port)
        Pp = 2 * np.pi * self.rp  # perimeter of the port
        alpha = np.sqrt(f_axis) * (0.95e-5 + 2.03e-5) * Pp / 2 / self.Sp
        kl = k * (1 + alpha * (1 - 1j))
        d0 = (0.6133 + 0.85) * self.rp  # length correction
        self.Zrad = rho * c / self.Sp * (1 / 4 * (kl * self.rp) ** 2 + 1j * kl * d0)
        self.Zp = s * self.Map + self.Zrad  # total impedance of the port (tube + radiation)

        # enclosure impedance
        self.Zab = gtb.parallel(1 / s / self.Cab, self.Rab, s * self.Map + self.Zrad)

        # total system impedance calculation
        ZaTot = driver.Zs
        Ps = driver.Ps

        # compute volume velocities
        Qs = Ps / (ZaTot + self.Zab)
        Qp = Qs * self.Zbox / (self.Zbox + self.Zp)

        p = 1j * k * rho * c * (Qs - Qp) * np.exp(-1j * k * 1) / 2 / np.pi / 1

        # total electrical impedance
        Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * (self.Zab + driver.Zaf))

        ## data plot
        fig, ax = plt.subplots(2, 1)
        ax[0].semilogx(f_axis, gtb.gain.SPL(p), label='SPL')
        ax[1].semilogx(f_axis, np.abs(Ze), label='Magnitude')

        ax[1].set(xlabel='Frequency [Hz]', ylabel='Impedance [Ohm]')
        ax[0].set(ylabel='SPL [dB] at 1 meter')
        ax[0].legend(loc='best')
        ax[1].legend(loc='best')
        for i in range(2):
            ax[i].grid(which='both')
        plt.subplots_adjust(bottom=0.25)

        # creation of text boxes
        volumeBox = fig.add_axes([0.2, 0.05, 0.075, 0.075])
        bVb = TextBox(volumeBox, "Vb (L): ")

        lengthBox = fig.add_axes([0.4, 0.05, 0.075, 0.075])
        bLp = TextBox(lengthBox, "Lp (cm): ")

        radiusBox = fig.add_axes([0.6, 0.05, 0.075, 0.075])
        bRp = TextBox(radiusBox, "rp (cm): ")

        sectionBox = fig.add_axes([0.8, 0.05, 0.075, 0.075])
        bSp = TextBox(sectionBox, r"Sp (cm$^2$): ")

        def updateVolume(expr):
            self.Vb = float(expr) * 1e-3  # back in m^3

            # box impedance
            self.Cab = self.Vb / rho / c ** 2  # compliance of the enclosure
            self.Rab = rho * c / eta / self.Vb
            self.Zbox = gtb.parallel(1 / s / self.Cab, self.Rab)

            # enclosure impedance
            self.Zab = gtb.parallel(1 / s / self.Cab, self.Rab, s * self.Map + self.Zrad)

            # update data
            updateZab(self.Zab, self.Zbox, self.Zp)
            return None

        def updatePortLength(expr):
            self.Lp = float(expr) * 1e-2  # back in m

            # box impedance
            # Cab = Vb / rho / c ** 2  # compliance of the enclosure
            # Rab = rho * c / eta / Vb
            # Zbox = parallel(1 / s / Cab, Rab)

            # port impedance
            self.rp = np.sqrt(self.Sp / pi)  # radius of the port
            self.Map = rho * self.Lp / self.Sp  # acoustical mass of the port
            self.Mal = 0.85 * 2 * self.rp  # length correction
            self.Map += self.Mal

            # radiation impedance (port)
            Pp = 2 * np.pi * self.rp  # perimeter of the port
            alpha = np.sqrt(f_axis) * (0.95e-5 + 2.03e-5) * Pp / 2 / self.Sp
            kl = k * (1 + alpha * (1 - 1j))
            d0 = (0.6133 + 0.85) * self.rp  # length correction
            self.Zrad = rho * c / self.Sp * (1 / 4 * (kl * self.rp) ** 2 + 1j * kl * d0)
            self.Zp = s * self.Map + self.Zrad  # total impedance of the port (tube + radiation)

            # enclosure impedance
            self.Zab = gtb.parallel(1 / s / self.Cab, self.Rab, s * self.Map + self.Zrad)

            # update data
            updateZab(self.Zab, self.Zbox, self.Zp)
            return None

        def updatePortRadius(expr):
            self.rp = float(expr) * 1e-2  # back in m
            self.Sp = np.pi * self.rp ** 2

            # # box impedance
            # Cab = Vb / rho / c ** 2  # compliance of the enclosure
            # Rab = rho * c / eta / Vb
            self.Zbox = gtb.parallel(1 / s / self.Cab, self.Rab)

            # port impedance
            self.Map = rho * self.Lp / self.Sp  # acoustical mass of the port
            self.Mal = 0.85 * 2 * self.rp  # length correction
            self.Map += self.Mal

            # radiation impedance (port)
            Pp = 2 * np.pi * self.rp  # perimeter of the port
            alpha = np.sqrt(f_axis) * (0.95e-5 + 2.03e-5) * Pp / 2 / self.Sp
            kl = k * (1 + alpha * (1 - 1j))
            d0 = (0.6133 + 0.85) * self.rp  # length correction
            self.Zrad = rho * c / self.Sp * (1 / 4 * (kl * self.rp) ** 2 + 1j * kl * d0)
            self.Zp = s * self.Map + self.Zrad  # total impedance of the port (tube + radiation)

            # enclosure impedance
            self.Zab = gtb.parallel(1 / s / self.Cab, self.Rab, s * self.Map + self.Zrad)

            # update data
            bSp.set_val(str(round(self.Sp * 1e4, 2)))
            updateZab(self.Zab, self.Zbox, self.Zp)
            return None

        def updatePortSection(expr):
            self.Sp = float(expr) * 1e-4  # back in m^2
            self.rp = np.sqrt(self.Sp / pi)  # radius of the port

            # box impedance
            # Cab = Vb / rho / c ** 2  # compliance of the enclosure
            # Rab = rho * c / eta / Vb
            self.Zbox = gtb.parallel(1 / s / self.Cab, self.Rab)

            # port impedance
            self.Map = rho * self.Lp / self.Sp  # acoustical mass of the port
            self.Mal = 0.85 * 2 * self.rp  # length correction
            self.Map += self.Mal

            # radiation impedance (port)
            Pp = 2 * np.pi * self.rp  # perimeter of the port
            alpha = np.sqrt(f_axis) * (0.95e-5 + 2.03e-5) * Pp / 2 / self.Sp
            kl = k * (1 + alpha * (1 - 1j))
            d0 = (0.6133 + 0.85) * self.rp  # length correction
            self.Zrad = rho * c / self.Sp * (1 / 4 * (kl * self.rp) ** 2 + 1j * kl * d0)
            self.Zp = s * self.Map + self.Zrad  # total impedance of the port (tube + radiation)

            # enclosure impedance
            self.Zab = gtb.parallel(1 / s / self.Cab, self.Rab, s * self.Map + self.Zrad)

            # update data
            bRp.set_val(str(round(self.rp * 1e2, 2)))
            updateZab(self.Zab, self.Zbox, self.Zp)
            return None

        def updateZab(Zab, Zbox, Zp):
            # compute volume velocities
            Qs = Ps / (ZaTot + Zab)
            Qp = Qs * Zbox / (Zbox + Zp)

            p = 1j * k * rho * c * (Qs - Qp) * np.exp(-1j * k * 1) / (2 * np.pi * 1)

            # total electrical impedance
            Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * (Zab + driver.Zaf))

            ## data plot
            ax[0].clear()
            ax[1].clear()
            ax[0].semilogx(f_axis, gtb.gain.SPL(p), label='SPL')
            ax[1].semilogx(f_axis, np.abs(Ze), label='Magnitude')

            ax[1].set(xlabel='Frequency [Hz]', ylabel='Impedance [Ohm]')
            ax[0].set(ylabel='SPL [dB] at 1 meter')
            ax[0].legend(loc='best')
            ax[1].legend(loc='best')
            for i in range(2):
                ax[i].grid(which='both')
            plt.subplots_adjust(bottom=0.25)
            return None

        # volume box (input in L)
        bVb.on_submit(updateVolume)
        bVb.set_val(str(round(self.Vb * 1e3, 2)))  # 1e3 to set in L

        # length box (input in cm)
        bLp.on_submit(updatePortLength)
        bLp.set_val(str(round(self.Lp * 1e2, 2)))  # 1e2 to set in cm

        # radius box (input in cm)
        bRp.on_submit(updatePortRadius)
        bRp.set_val(str(round(self.rp * 1e2, 2)))  # 1e2 to set in cm

        # section box (input in cm^2)
        bSp.on_submit(updatePortSection)
        bSp.set_val(str(round(self.Sp * 1e4, 2)))  # 1e4 to set in cm^2
        return bVb, bLp, bRp, bSp, plt.show()

# def loadLPM(lpmfile, freq_array, U=1, LeZero=False,
#             number_of_drivers=1,
#             wiring='parallel',
#             c=air.c,
#             rho=air.rho):
#     """
#     Return electro_acoustic_driver object from LPM file (Klippel measurements)
#     :param lpmfile:
#     :param freq_array:
#     :param U:
#     :param c:
#     :param rho:
#     :return:
#     """
#     data = pd.read_csv(lpmfile, sep="\t", header=0,
#                        names=["Parameters", "value", "unit", "description"],
#                        encoding='unicode_escape')
#     idxRe = data.index[data['Parameters'] == 'Re'].to_list()[0]
#     idxLe = data.index[data['Parameters'] == 'Le'].to_list()[0]
#     idxCms = data.index[data['Parameters'] == 'Cms'].to_list()[0]
#     idxMms = data.index[data['Parameters'] == 'Mms'].to_list()[0]
#     idxRms = data.index[data['Parameters'] == 'Rms'].to_list()[0]
#     idxSd = data.index[data['Parameters'] == 'Sd'].to_list()[0]
#     idxBl = data.index[data['Parameters'] == 'Bl'].to_list()[0]

#     # print('idxLe: ', idxLe)
#     # print(type(data.iat[idxLe, 1]))

#     Re = float(data.iat[idxRe, 1])
#     Le = float(data.iat[idxLe, 1]) * 1e-3
#     Mms = float(data.iat[idxMms, 1]) * 1e-3
#     Rms = float(data.iat[idxRms, 1])
#     Cms = float(data.iat[idxCms, 1]) * 1e-3
#     Bl = float(data.iat[idxBl, 1])
#     Sd = float(data.iat[idxSd, 1]) * 1e-4

#     if LeZero is True:
#         Le = 0

#     if number_of_drivers > 1:
#         if wiring == 'parallel':
#             n = number_of_drivers
#             drv = electroAcousticDriver(U, Le/n, Re/n, Cms/n, Mms*n, Rms*n, Bl, Sd*n, freq_array, c, rho)
#         elif wiring == 'series':
#             n = number_of_drivers
#             drv = electroAcousticDriver(U, Le*n, Re*n, Cms/n, Mms*n, Rms*n, Bl*n, Sd*n, freq_array, c, rho)
#         else:
#             ValueError("'wiring' must be either 'parallel' or 'series'.")
#     else:
#         drv = electroAcousticDriver(U, Le, Re, Cms, Mms, Rms, Bl, Sd, freq_array, c, rho)
#     return drv

def loadLPM(lpmfile, freq_array, U=1, LeZero=False,
            number_of_drivers=1,
            wiring='parallel',
            c=air.c,
            rho=air.rho):
    
    # define loader based on extension
    _, extension = os.path.splitext(lpmfile)
    if extension == ".qsp":
        loader    = lpl.qspeaker_lp_loader
        weight_Le  = 1e-3
        weight_Sd  = 1
        weight_Mms = 1
        weight_Cms = 1
    elif extension == ".sdrv":
        loader = lpl.speakerSim_lp_loader
        weight_Le  = 1
        weight_Sd  = 1
        weight_Mms = 1
        weight_Cms = 1
    elif extension == ".wdr":
        loader = lpl.winSd_lp_loader
        weight_Le  = 1
        weight_Sd  = 1
        weight_Mms = 1
        weight_Cms = 1
    elif extension == ".bastaelement":
        loader = lpl.basta_lp_loader
        weight_Le  = 1
        weight_Sd  = 1
        weight_Mms = 1
        weight_Cms = 1
    elif extension == ".txt":
        with open(lpmfile, 'r') as file:
            first_line = file.readline().strip()
        if first_line == 'Electrical Parameters':
            loader = lpl.klippel_lp_loader
            weight_Le  = 1e-3
            weight_Sd  = 1e-4
            weight_Mms = 1e-3
            weight_Cms = 1e-3
        else:
            loader = lpl.hornResp_lp_loader
            weight_Le  = 1e-3
            weight_Sd  = 1e-4
            weight_Mms = 1e-3
            weight_Cms = 1
    
    # create driver object
    data = loader(lpmfile)
    Le = data["Le"] * weight_Le
    Re = data["Re"]
    Cms = data["Cms"] * weight_Cms
    Mms = data["Mms"] * weight_Mms
    Rms = data["Rms"]
    Bl = data["Bl"]
    Sd = data["Sd"] * weight_Sd
    
    if LeZero is True:
        Le = 0
    
    
    if number_of_drivers > 1:
        if wiring == 'parallel':
            n = number_of_drivers
            drv = electroAcousticDriver(U, Le/n, Re/n, Cms/n, Mms*n, 
                                        Rms*n, Bl, Sd*n, freq_array, c, rho)
        elif wiring == 'series':
            n = number_of_drivers
            drv = electroAcousticDriver(U, Le*n, Re*n, Cms/n, 
                                        Mms*n, Rms*n, Bl*n, Sd*n, 
                                        freq_array, c, rho)
        else:
            ValueError("'wiring' must be either 'parallel' or 'series'.")
    else:
        drv = electroAcousticDriver(U, Le, Re, Cms, Mms, Rms, 
                                    Bl, Sd, freq_array, c, rho)
    return drv
    
    
    
