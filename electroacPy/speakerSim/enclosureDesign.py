#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:56:03 2023

@author: tom.munoz
"""

from electroacPy.globalVariables import air
from electroacPy.speakerSim.electroAcousticDriver import electroAcousticDriver as ead
from generalToolbox.freqop import freq_log10 as f10
from generalToolbox.freqop import laplace
from generalToolbox import parallel
from generalToolbox.geometry import compute_circle_point_cloud_surface_area as ccpcsa
import matplotlib.pyplot as plt
import numpy as np

pi = np.pi

class speakerBox:
    def __init__(self, Vb, Lp=False, Sp=False, rp=False, flanged=False,
                 frequencyRange=f10(20, 2500, 50),
                 eta=1e-5, c=air.c, rho=air.rho):
        """
        Create an enclosure with the possibility to set electroAcousticDriver objects in it.

        :param Vb: float,
            volume
        :param Lp: float,
            port's length, if False enclosure is sealed
        :param Sp: float,
            port's area
        :param rp: float,
            port's radius
        :param frequencyRange: float,
            study range
        """
        # default enclosure parameters
        self.Vb = Vb
        self.Lp = Lp  # if not given, enclosure is sealed
        self.Sp = Sp  # must be given (or rp) if Lp is defined
        self.rp = rp
        self.frequencyRange = frequencyRange
        self.eta = eta # losses coefficient -> box leakage
        self.c = c
        self.rho = rho
        w = 2*pi*frequencyRange
        s = 1j*w
        k = w / self.c
        # Check if enclosure is ported or not
        if Lp != False:  # PORTED ENCLOSURE - From Beranek
            self.isPorted = True
            if self.Sp is False and self.rp is False:
                raise ValueError('Cross sectional area Sp or radius rp must be defined.')
            elif self.Sp is False:
                self.Sp = np.pi * self.rp ** 2
                self.rp = rp
            elif self.rp is False:
                self.rp = np.sqrt(Sp / np.pi)
                self.Sp = Sp

            # box impedance
            Cab = Vb / self.rho / self.c ** 2  # compliance of the enclosure
            Rab = self.rho * self.c / eta / Vb
            Zbox = parallel(1 / s / Cab, Rab)

            # port impedance
            if flanged is True:
                ll = 8*self.rp/3/np.pi
            else:
                ll = 2*self.rp/np.pi
            Map = (self.Lp + 0.64*self.rp) * self.rho / self.Sp
            Rap = np.sqrt(2*w*self.rho*1.86e-5) / self.Sp * (ll / self.rp + 0.7)
            # Zap = s*Map + Rap

            # port radiation
            # RAR2 = np.pi * self.frequencyRange**2 * self.rho / self.c
            # MA2 = (self.Lp + 0.64*self.rp)*self.rho / (np.pi*self.rp**2)
            RAR2 = 0.159 * w ** 2 * self.rho / self.c
            MAR2 = 0.270 * self.rho / self.rp
            Zap =  s*MAR2 + RAR2 + s*Map + Rap

            # total enclosure impedance
            self.Zab = Zbox # box impedance
            self.Zap = Zap
            self.Za = parallel(Zap, Zbox) # total acoustic impedance

            # port impedance
            # Map = self.rho * Lp / self.Sp  # acoustical mass of the port
            # Mal = 0.85 * 2 * self.rp  # length correction
            # # Mal = 0.54*rho/2/rp
            # Map += Mal

            # radiation impedance (port)
            # Pp = 2 * np.pi * self.rp  # perimeter of the port
            # alpha = np.sqrt(frequencyRange) * (0.95e-5 + 2.03e-5) * Pp / 2 / self.Sp
            # kl = k * (1 + alpha * (1 - 1j))
            # d0 = (0.6133 + 0.85) * self.rp  # length correction
            # Zrad = self.rho * self.c / self.Sp * (1 / 4 * (kl * self.rp) ** 2 + 1j * kl * d0)
            # Zrad = (0.270*self.rho/self.rp)*s + 0.159*w**2*self.rho / self.c
            # Zp = s * Map + Zrad  # total impedance of the port (tube + radiation)
            # enclosure impedance
            # self.Zab = parallel(1/s/Cab, Rab)
            # self.Zap = s*Map + Zrad
            # self.Za = parallel(self.Zab, self.Zap)


        else: # SEALED ENCLOSURE
            if type(Vb) == float or type(Vb) == int or type(Vb) == np.float64:
                Cab = Vb / self.rho / self.c ** 2  # compliance of the enclosure
                if eta == 0:
                    self.Zab = 1 / (s * Cab)
                    self.Za = self.Zab  # total acoustic impedance
                else:
                    Rab = self.rho * self.c / eta / Vb
                    self.Zab =  1 / ((s * Cab) + 1 / Rab)
                    self.Za = self.Zab   # total acoustic impedance
            elif type(Vb) == np.ndarray: # if the impedance if defined by the user (could be ported?)
                self.Zab = Vb
                self.Za = self.Zab       # total acoustic impedance
            self.isPorted = False


        # to store results
        self.v   = np.zeros(len(frequencyRange), dtype=complex)
        self.vp  = np.zeros(len(frequencyRange), dtype=complex)
        self.Ze  = np.zeros(len(frequencyRange), dtype=complex) # total electrical impedance of driver in enclosure

        # reference to driver
        self.whichDriver = None

        # acoustic simulation reference
        self.isFEM     = False
        self.isBEM     = False
        self.ref2bem   = None   # is it referenced to bem mesh ?
        self.poly_data = False  # is class from polytech?


    def getDriverResponse(self, driver, Nd=1, wiring='parallel'):
        """
        Compute driver response in enclosure. Can setup multiple drivers inside common enclosure.

        Parameters
        ----------
        driver : electroAcousticDriver object
        Nd : int,
            number of drivers inside enclosure. Default is 1
        wirint : str,
            how drivers are wired inside enclosure. Default is 'parallel'

        Return
        ------
        None
        """

        if driver.identifier == 'EAC':
            d = driver
            n = Nd
            if wiring == 'parallel':
                drvTmp = ead(d.U, d.Le/n, d.Re/n, d.Cms/n, d.Mms*n, d.Rms*n, d.Bl, d.Sd*n, d.f_array, d.c, d.rho)
            elif wiring == 'series':
                drvTmp = ead(d.U, d.Le*n, d.Re*n, d.Cms/n, d.Mms*n, d.Rms*n, d.Bl*n, d.Sd*n, d.f_array, d.c, d.rho)
            else:
                raise ValueError("Speaker wiring not understood. Accepted values are: 'parallel', 'series'.")

            if self.isPorted is False:
                self.v, self.Ze = sealedBoxQ_V2(drvTmp, self.Za)
            else:
                self.v, self.vp, self.Ze = portedBoxQ_V2(drvTmp, self.Zab, self.Zap, self.Sp)
            driver.inBox = True
            
        elif driver.identifier == 'PLV':
            l = driver
            n = Nd
            if self.isPorted is False:
                raise TypeError('No need to compute LEM spkbox. Measurement should be taken in enclosure')
                #self.v, self.Ze = sealedBoxQ(l, self.Vb, velocity=True, impedance=True)
            else:
                self.v, self.Ze = portedBoxQ(l, self.Vb, self.Lp, self.Sp,
                                             velocity=True, impedance=True, eta=self.eta, c=self.c, rho=self.rho)
        return None

    def plotZe(self):
        """
        Plot the electrical impedance ZeTot in both modulus and phase of driver in enclosure.

        Returns
        -------
        Matplotlib figure

        """
        try:
            fig, ax = plt.subplots()
            ax.semilogx(self.frequencyRange, np.abs(self.Ze),
                        label="modulus")
            ax.legend(loc='upper left')
            ax.grid(linestyle='dotted', which='both')
            ax.set(ylabel="|Z| [Ohm]", xlabel="Frequency [Hz]")
            ax2 = ax.twinx()
            ax2.semilogx(self.frequencyRange, np.angle(self.Ze), '--',
                         label="phase")
            ax2.legend(loc='upper right')
            ax2.set(ylabel="phase [rads]", ylim=[-pi / 2, pi / 2])
            ax2.grid()
            plt.tight_layout()
        except:
            print('No driver set in specified enclosure.')
        return plt.show()

    def plotXVA(self):
        """
        Plot displacement, velocity and acceleration response of driver in enclosure

        Returns
        -------
        Matplotlib figure
        """
        w = self.frequencyRange * 2 * pi
        try:
            fig, ax = plt.subplots(3, 1)
            ax[0].semilogx(self.frequencyRange, np.abs(self.v/1j/w * 1e3), label='Displacement')
            ax[1].semilogx(self.frequencyRange, np.abs(self.v), label='Velocity')
            ax[2].semilogx(self.frequencyRange, np.abs(self.v*1j*w), label='Acceleration')
            ax[2].set(xlabel="Frequency [Hz]")
            ax[0].set(ylabel="mm")
            ax[1].set(ylabel="m/s")
            ax[2].set(ylabel="m/s^2")
            for i in range(3):
                ax[i].grid(which='both')
                ax[i].legend(loc='best')
            plt.tight_layout()
        except:
            print('No driver set in specified enclosure.')
        return None


## New alignment function
def sealedBoxQ_V2(driver, Zab):
    """
    Compute the volume velocity at the driver output (sealed enclosure). Takes Vb as the enclosure volume if it is float
     or int. Takes Vb as the acoustical impedance if it is passed as a ndarray (e.g. Zab from FEM calculations)

    Parameters
    ----------
    driver : electro_acoustic_driver class.
    Vb : float, int, complex
        if float or int: Volume of the enclosure. m**3
        if complex: acoustical impedance seen by the driver
    eta: float, optional
        acoustical loss factor in enclosure. Default is 0
    c : TYPE, optional
        celerity of sound in air. The default is self.c.
    rho : TYPE, optional
        density of air. The default is self.c.

    Returns
    -------
    Qs : Volume velocity at the driver.
    if impedance if True:
        Ze : Total electrical impedance of the speaker  ()
    """
    if driver.identifier == 'EAC':
        Qs = driver.Ps / (driver.Zs + Zab)
        Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * (Zab + driver.Zaf))
        v = Qs / driver.Sd
    elif driver.identifier == 'PLV':
        raise TypeError('If enclosure is sealed, no need to get LEM data (in-enclosure acceleration data).')
    return v, Ze

def portedBoxQ_V2(driver, Zab, Zap, Sp):
    """
    Compute the volume velocities at the driver and port of a bass-reflex
    enclosure. Will also output the total impedance expressed in the electrical
    domain.

    Parameters
    ----------
    driver : electro_acoustic_driver class.
    Zab : Impedance of the "box" part of enclosure
    Zap : Impedance of port
    Sp : Cross-section area of the port. m**2

    Returns
    -------
    Qs : Volume velocity at the driver.
    Qv : Volume velocity at the port.
    if impedance is True:
        Ze : Total electrical impedance with acoustical impedance

    """
    Za = parallel(Zab, Zap) # total acoustic impedance of the ported enclosure
    Qs = driver.Ps / (driver.Zs + Za)
    Qv = -Qs * Zab / (Zab + Zap)
    v = Qs / driver.Sd
    vp = Qv / Sp

    Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * (Za + driver.Zaf))
    return v, vp, Ze


## New approach
class speakerBox_v2:
    def __init__(self, Vb, frequencyRange=f10(20, 2500, 50),
                 eta=1e-5, c=air.c, rho=air.rho, **kwargs):
        """
        Setup a louspeaker enclosure as a lumped-element object.

        :param Vb: float,
            Back volume (set behind the drive unit)
        :param frequencyRange: array,

        :param eta:
        :param c:
        :param rho:
        """
        # simulation parameters
        self.c = c
        self.rho = rho
        self.eta = eta
        self.frequencyRange = frequencyRange
        self.s = laplace(frequencyRange)


        # mandatory parameter
        # 0. Sealed enclosure
        self.Vb = Vb
        self.config = "sealed"
        if kwargs:
            self.detectConfig(**kwargs)

        ## PROCESS
        self.computeImpedance()

        ## TO STORE RESULTS
        # velocity
        self.v = np.zeros(len(frequencyRange), dtype=complex)
        self.vp = np.zeros(len(frequencyRange), dtype=complex)
        self.vp2 = np.zeros(len(frequencyRange), dtype=complex)
        self.vpr = np.zeros(len(frequencyRange), dtype=complex)
        # volume velocity
        self.Q = np.zeros(len(frequencyRange), dtype=complex)
        self.Qp = np.zeros(len(frequencyRange), dtype=complex)
        self.Qp2 = np.zeros(len(frequencyRange), dtype=complex)
        self.Qpr = np.zeros(len(frequencyRange), dtype=complex)
        # impedance
        self.Ze = np.zeros(len(frequencyRange), dtype=complex)  # total electrical impedance of driver in enclosure

        # reference to driver
        self.whichDriver = None
        self.Nd = 1
        self.wiring = "parallel"

        # acoustic simulation reference
        self.isFEM = False
        self.isBEM = False
        self.ref2bem = None  # is it referenced to bem mesh ?
        self.poly_data = False  # is class from polytech?


    def detectConfig(self, **kwargs):
        """
        Set enclosure config depending on input arguments.
        :param kwargs:
        :return:
        """
        ## 6th order bandpass with ports
        if ("Lp" in kwargs and ("Sp" in kwargs or "rp" in kwargs) and
                "Vf" in kwargs and "Lp2" in kwargs and ("Sp2" in kwargs or "rp2" in kwargs)):
            self.config = "bandpass_2"
            for key, value in kwargs.items():
                setattr(self, key, value)
                if key == "Sp":
                    self.rp = np.sqrt(value / np.pi)
                elif key == "rp":
                    self.Sp = np.pi * value**2
                elif key == "Sp2":
                    self.rp2 = np.sqrt(value / np.pi)
                elif key == "rp2":
                    self.Sp2 = np.pi * value**2

        ## 6th order bandpass with passive radiators
        elif ("Mmd" in kwargs and "Cmd" in kwargs and "Rmd" in kwargs and
              "Vf" in kwargs and "Mmd2" in kwargs and "Cmd2" in kwargs and "Rmd2" in kwargs and
              "Sd" in kwargs and "Sd2" in kwargs):
            self.config = "bandpass_pr_2"
            for key, value in kwargs.items():
                setattr(self, key, value)

        ## 4th order bandpass with port
        elif "Lp" in kwargs and ("Sp" in kwargs or "rp" in kwargs) and "Vf" in kwargs:
            self.config = "bandpass"
            for key, value in kwargs.items():
                setattr(self, key, value)
                if key == "Sp":
                    self.rp = np.sqrt(value / np.pi)
                elif key == "rp":
                    self.Sp = np.pi * value**2

        ## 4th order bandpass with passive radiator
        elif ("Mmd" in kwargs and "Cmd" in kwargs and
              "Rmd" in kwargs and "Vf" in kwargs and "Sd" in kwargs):
            self.config = "bandpass_pr"
            for key, value in kwargs.items():
                setattr(self, key, value)

        ## ported enclosure
        elif "Lp" in kwargs and ("Sp" in kwargs or "rp" in kwargs):
            self.config = "vented"
            for key, value in kwargs.items():
                setattr(self, key, value)
                if key == "Sp":
                    self.rp = np.sqrt(value / np.pi)
                elif key == "rp":
                    self.Sp = np.pi * value**2

        ## passive radiator enclosure
        elif "Mmd" in kwargs and "Cmd" in kwargs and "Rmd" in kwargs and "Sd" in kwargs:
            self.config = "passiveRadiator"
            for key, value in kwargs.items():
                setattr(self, key, value)
        else:
            raise ValueError("Invalid configuration: The provided kwargs do not match any valid configuration.\n"
                             "1. Ported enclosure: Lp and (Sp or rp), \n"
                             "2. Passive radiator: Mmd, Cmd and Rmd, \n"
                             "3. 4th order bandpass (port): Vf, Lp and (Sp or rp), \n"
                             "4. 4th order bandpass (passive radiator): Vf, Mmd, Cmd and Rmd, \n"
                             "5. 6th order bandpass (port): Vf, Lp, (Sp or rp), Lp2 and (Sp2 or rp2), \n"
                             "6. 6th order bandpass (passive radiator): Vf, Mmd, Cmd, Rmd, Mmd2, Cmd2 and Rmd2.")
        return None

    ### =============================
    ## ACOUSTICAL IMPEDANCE FUNCTIONS
    def computeImpedance(self):
        """
        Just a small function to send impedance calc to proper function
        :return: None
        """

        if self.config == "sealed":
            self.Za_in = self.sealed_box()

        elif self.config == "vented":
            self.Za_in, self.Zab, self.Zap = self.vented_box()

        elif self.config == "passiveRadiator":
            self.Za_in, self.Zab, self.Zapr = self.pr_box()

        elif self.config == "bandpass":
            self.Za_in, self.Zab, self.Zaf, self.Zap = self.bp4_box()

        elif self.config == "bandpass_2":
            self.Za_in, self.Zab, self.Zaf, self.Zap, self.Zap2 = self.bp6_box()

        elif self.config == "bandpass_pr":
            self.Za_in, self.Zab, self.Zaf, self.Zapr = self.bp4_pr_box()
        return None

    def sealed_box(self,):
        """
        Create a sealed enclosure.
        :return:
        """
        Cab = self.Vb / self.rho / self.c ** 2          # compliance of the enclosure
        Rab = self.rho * self.c / self.eta / self.Vb    # losses through leakage
        Zab = parallel(1 / self.s / Cab, Rab)
        return Zab

    def vented_box(self,):
        """
        Create ported enclosure.
        :return:
        """
        w = 2*np.pi*self.frequencyRange

        # box impedance
        Cab = self.Vb / self.rho / self.c ** 2  # compliance of the enclosure
        Rab = self.rho * self.c / self.eta / self.Vb
        Zab = parallel(1 / self.s / Cab, Rab)

        # port impedance
        # if flanged is True:
        #     ll = 8 * self.rp / 3 / np.pi
        # else:
        #     ll = 2 * self.rp / np.pi
        ll = 2 * self.rp / np.pi
        Map = (self.Lp + 0.64 * self.rp) * self.rho / self.Sp
        Rap = np.sqrt(2 * w * self.rho * 1.86e-5) / self.Sp * (ll / self.rp + 0.7)

        # port radiation
        RAR2 = 0.159 * w ** 2 * self.rho / self.c
        MAR2 = 0.270 * self.rho / self.rp
        Zap = self.s * MAR2 + RAR2 + self.s * Map + Rap

        # total enclosure impedance
        Za_in = parallel(Zap, Zab)  # total acoustic impedance
        return Za_in, Zab, Zap

    def pr_box(self,):
        """
        Create a passive radiator enclosure.
        :return:
        """
        w = 2 * np.pi * self.frequencyRange

        # box impedance
        Cab = self.Vb / self.rho / self.c ** 2  # compliance of the enclosure
        Rab = self.rho * self.c / self.eta / self.Vb
        Zab = parallel(1 / self.s / Cab, Rab)

        # passive radiator's impedance
        RAR = 0.159 * w ** 2 * self.rho / self.c         # piston's radiation impedance - air resistance
        MAR = 0.270 * self.rho / np.sqrt(self.Sd/np.pi)  # piston's radiation impedance - air mass
        Zapr = ((self.Rmd + self.s * self.Mmd + 1 / (self.s * self.Cmd)) / self.Sd ** 2) + self.s * MAR + RAR

        # total enclosure impedance
        Za_in = parallel(Zapr, Zab)  # total acoustic impedance
        return Za_in, Zab, Zapr

    def bp4_box(self,):
        """
        Create a 4th order bandpass enclosure (ports).
        :return:
        """

        w = 2 * np.pi * self.frequencyRange

        # back enclosure impedance
        Cab = self.Vb / self.rho / self.c ** 2  # compliance of the enclosure
        Rab = self.rho * self.c / self.eta / self.Vb
        Zab = parallel(1 / self.s / Cab, Rab)

        # front enclosure impedance
        Caf = self.Vf / self.rho / self.c**2
        Raf = self.rho * self.c / self.eta / self.Vf
        Zaf = parallel(1/self.s/Caf, Raf)
        ll = 2 * self.rp / np.pi
        Map = (self.Lp + 0.64 * self.rp) * self.rho / self.Sp
        Rap = np.sqrt(2 * w * self.rho * 1.86e-5) / self.Sp * (ll / self.rp + 0.7)

        # port radiation
        RAR = 0.159 * w ** 2 * self.rho / self.c
        MAR = 0.270 * self.rho / self.rp
        Zap = self.s * MAR + RAR + self.s * Map + Rap
        Za_front = parallel(Zap, Zaf)
        Za_in = parallel(Za_front, Zab)
        return Za_in, Zab, Zaf, Zap

    def bp6_box(self,):
        """
        Create a 6th order bandpass enclosure (ports).
        :return:
        """
        w = 2 * np.pi * self.frequencyRange

        # back enclosure impedance
        Cab = self.Vb / self.rho / self.c ** 2  # compliance of the enclosure
        Rab = self.rho * self.c / self.eta / self.Vb
        Zab = parallel(1 / self.s / Cab, Rab)
        ll2 = 2 * self.rp2 / np.pi
        Map2 = (self.Lp2 + 0.64 * self.rp2) * self.rho / self.Sp2
        Rap2 = np.sqrt(2 * w * self.rho * 1.86e-5) / self.Sp2 * (ll2 / self.rp2 + 0.7)
        RAR2 = 0.159 * w ** 2 * self.rho / self.c
        MAR2 = 0.270 * self.rho / self.rp2
        Zap2 = self.s * MAR2 + RAR2 + self.s * Map2 + Rap2

        # front enclosure impedance
        Caf = self.Vf / self.rho / self.c**2
        Raf = self.rho * self.c / self.eta / self.Vf
        Zaf = parallel(1/self.s/Caf, Raf)
        ll = 2 * self.rp / np.pi
        Map = (self.Lp + 0.64 * self.rp) * self.rho / self.Sp
        Rap = np.sqrt(2 * w * self.rho * 1.86e-5) / self.Sp * (ll / self.rp + 0.7)

        # port radiation
        RAR = 0.159 * w ** 2 * self.rho / self.c
        MAR = 0.270 * self.rho / self.rp
        Zap = self.s * MAR + RAR + self.s * Map + Rap
        Za_front = parallel(Zap, Zaf)
        Za_back = parallel(Zap2, Zab)
        Za_in = parallel(Za_front, Zab)
        return Za_in, Zab, Zaf, Zap, Zap2

    def bp4_pr_box(self, ):
        """
        Create a 4th order bandpass enclosure (passive radiator).
        :return:
        """

        w = 2 * np.pi * self.frequencyRange

        # box impedance - front
        Caf = self.Vf / self.rho / self.c ** 2  # compliance of the enclosure
        Raf = self.rho * self.c / self.eta / self.Vf
        Zaf = parallel(1 / self.s / Caf, Raf)

        # box impedance - back
        Cab = self.Vb / self.rho / self.c ** 2  # compliance of the enclosure
        Rab = self.rho * self.c / self.eta / self.Vb
        Zab = parallel(1 / self.s / Cab, Rab)

        # passive radiator's impedance
        RAR = 0.159 * w ** 2 * self.rho / self.c         # piston's radiation impedance - air resistance
        MAR = 0.270 * self.rho / np.sqrt(self.Sd/np.pi)  # piston's radiation impedance - air mass
        Zapr = ((self.Rmd + self.s * self.Mmd + 1 / (self.s * self.Cmd)) / self.Sd ** 2) + self.s * MAR + RAR

        # total enclosure impedance
        Za_front = parallel(Zapr, Zaf)
        Za_in = parallel(Za_front, Zab)
        return Za_in, Zab, Zaf, Zapr

    ### ======================================
    ## IMPEDANCE AND VOLUME VELOCITY FUNCTIONS
    def sealed_Q(self, driver):
        Q = driver.Ps / (driver.Zs + self.Za_in)
        Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * (self.Za_in + driver.Zaf))
        
        # load parameters (complexify a lot for a simple sealed enclosure)
        # U = driver.U
        # Le = driver.Le
        # Re = driver.Re
        # Rms = driver.Rms
        # Cms = driver.Cms
        # Mms = driver.Mms
        # Sd = driver.Sd
        # Bl = driver.Bl
        
        # from electroacPy.circuitSolver.circuit import circuit
        # import electroacPy.circuitSolver.components.acoustic as compa
        
        # drv = compa.EAD(1, 2, 3, U, Le, Re, Cms, Mms, Rms, Bl, Sd)
        # box = compa.cavity(3, 0, self.Vb, c=self.c, rho=self.rho)
        # rad = compa.radiator(2, 0, Sd, c=self.c, rho=self.rho)
        # cir = circuit(self.frequencyRange)
        # cir.addComponent([drv, box, rad])
        # cir.run()
        
        # Q = cir.getPotential(2) * rad.G
        # Ze = 1
        return Q, Ze

    def vented_Q(self, driver):
        Q = driver.Ps / (driver.Zs + self.Za_in)
        Qp = - Q * self.Zab / (self.Zab + self.Zap)
        Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * (self.Za_in + driver.Zaf))
        return Q, Qp, Ze

    def pr_Q(self, driver):
        Q = driver.Ps / (driver.Zs + self.Za_in)
        Qpr = - Q * self.Zab / (self.Zab + self.Zapr)
        Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * (self.Za_in + driver.Zaf))
        return Q, Qpr, Ze

    def bp4_Q(self, driver):
        Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 *
                                           (self.Za_in + driver.Zac + driver.Zas))
        # MNA approach -- we solve AX=Z
        M = 1 # number of independant voltage source (1 driver)
        N = 4 # number of nodes

        # matrix initialization (see QUCS' online documentation)
        A = np.zeros([M+N, M+N], dtype=complex)
        G = np.zeros([N, N], dtype=complex)
        Q = np.zeros(len(self.frequencyRange), dtype=complex)
        Qp = np.zeros(len(self.frequencyRange), dtype=complex)

        # let's solve the system
        for i in range(len(self.frequencyRange)):
            Zad = driver.Zac[i]+driver.Zas[i]
            Zap, Zab, Zaf = self.Zap[i], self.Zab[i], self.Zaf[i]
            Ps = driver.Ps[i]
            G[0, 0] = 1/Zad
            G[1, 1] = 1/Zad + 1/Zaf + 1/Zap
            G[2, 2] = 1/Zad
            G[3, 3] = 1/Zad + 1/Zab
            G[0, 1], G[1, 0] = -1/Zad, -1/Zad
            G[2, 3], G[3, 2] = -1/Zad, -1/Zad
            B = np.array([[1, 0, -1, 0]]).T
            C = B.T
            D = 0
            Z = np.array([[0, 0, 0, 0, Ps]]).T
            A[0:4, 0:4] = G
            A[0:4, -1:] = B
            A[-1:, 0:4] = C
            A[-1, -1] = D
            # print(A)
            X = np.linalg.inv(A)@Z
            Qp[i] = X[1, 0] / Zap
            Q[i] = X[-1, 0]
        return Q, Qp, Ze

    def bp6_Q(self, driver):
        Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 *
                                           (self.Za_in + driver.Zac + driver.Zas))

        # MNA approach -- we solve AX=Z
        M = 1  # number of independant voltage source (1 driver)
        N = 4  # number of nodes

        # matrix initialization (see QUCS' online documentation)
        A = np.zeros([M + N, M + N], dtype=complex)
        G = np.zeros([N, N], dtype=complex)
        Q = np.zeros(len(self.frequencyRange), dtype=complex)
        Qp = np.zeros(len(self.frequencyRange), dtype=complex)
        Qp2 = np.zeros(len(self.frequencyRange), dtype=complex)

        # let's solve the system
        for i in range(len(self.frequencyRange)):
            Zad = driver.Zac[i] + driver.Zas[i]
            Zap, Zap2, Zab, Zaf = self.Zap[i], self.Zap2[i], self.Zab[i], self.Zaf[i]
            Ps = driver.Ps[i]
            G[0, 0] = 1 / Zad
            G[1, 1] = 1 / Zad + 1 / Zaf + 1 / Zap
            G[2, 2] = 1 / Zad
            G[3, 3] = 1 / Zad + 1 / Zab + 1 / Zap2
            G[0, 1], G[1, 0] = -1 / Zad, -1 / Zad
            G[2, 3], G[3, 2] = -1 / Zad, -1 / Zad
            B = np.array([[1, 0, -1, 0]]).T
            C = B.T
            D = 0
            Z = np.array([[0, 0, 0, 0, Ps]]).T
            A[0:4, 0:4] = G
            A[0:4, -1:] = B
            A[-1:, 0:4] = C
            A[-1, -1] = D

            X = np.linalg.inv(A) @ Z
            Qp[i] = X[1, 0] / Zap
            Qp2[i] = X[3, 0] / Zap2
            Q[i] = X[-1, 0]

        return Q, Qp, Qp2, Ze


    def bp4_pr_Q(self, driver):
        # MNA approach -- we solve AX=Z
        M = 1  # number of independant voltage source (1 driver)
        N = 4  # number of nodes

        # matrix initialization (see QUCS' online documentation)
        A = np.zeros([M + N, M + N], dtype=complex)
        G = np.zeros([N, N], dtype=complex)
        Q = np.zeros(len(self.frequencyRange), dtype=complex)
        Qpr = np.zeros(len(self.frequencyRange), dtype=complex)

        # let's solve the system
        for i in range(len(self.frequencyRange)):
            Zad = driver.Zac[i] + driver.Zas[i]
            Zapr, Zab, Zaf = self.Zapr[i], self.Zab[i], self.Zaf[i]
            Ps = driver.Ps[i]
            G[0, 0] = 1 / Zad
            G[1, 1] = 1 / Zad + 1 / Zaf + 1 / Zapr
            G[2, 2] = 1 / Zad
            G[3, 3] = 1 / Zad + 1 / Zab
            G[0, 1], G[1, 0] = -1 / Zad, -1 / Zad
            G[2, 3], G[3, 2] = -1 / Zad, -1 / Zad
            B = np.array([[1, 0, -1, 0]]).T
            C = B.T
            D = 0
            Z = np.array([[0, 0, 0, 0, Ps]]).T
            A[0:4, 0:4] = G
            A[0:4, -1:] = B
            A[-1:, 0:4] = C
            A[-1, -1] = D
            # print(A)
            X = np.linalg.inv(A) @ Z
            Qpr[i] = X[1, 0] / Zapr
            Q[i] = X[-1, 0]

        Za_in = Ps / Q
        Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 *
                                           (Za_in + driver.Zac + driver.Zas))
        return Q, Qpr, Ze

    ## SET DRIVERS IN ENCLOSURE
    def getDriverResponse(self, driver, Nd=1, wiring='parallel'):
        d = driver
        n = Nd
        self.Nd = Nd
        self.wiring = wiring

        if wiring == 'parallel':
            drvTmp = ead(d.U, d.Le / n, d.Re / n, d.Cms / n, d.Mms * n, d.Rms * n, d.Bl, d.Sd * n, d.f_array, d.c,
                         d.rho)
        elif wiring == 'series':
            drvTmp = ead(d.U, d.Le * n, d.Re * n, d.Cms / n, d.Mms * n, d.Rms * n, d.Bl * n, d.Sd * n, d.f_array, d.c,
                         d.rho)
        else:
            raise ValueError("Speaker wiring not understood. Accepted values are: 'parallel', 'series'.")

        ## Compute volume / particule velocity -> add enclosure config here
        if self.config == "sealed":
            self.Q, self.Ze = self.sealed_Q(drvTmp)
            self.v = self.Q / d.Sd
        elif self.config == "vented":
            self.Q, self.Qp, self.Ze = self.vented_Q(drvTmp)
            self.v, self.vp = self.Q / drvTmp.Sd, self.Qp / self.Sp
        elif self.config == "passiveRadiator":
            self.Q, self.Qpr, self.Ze = self.pr_Q(drvTmp)
            self.v, self.vpr = self.Q / drvTmp.Sd, self.Qpr / self.Sd
        elif self.config == "bandpass":
            self.Q, self.Qp, self.Ze = self.bp4_Q(drvTmp)
            self.v, self.vp = self.Q / drvTmp.Sd, self.Qp / self.Sp
        elif self.config == "bandpass_2":
            self.Q, self.Qp, self.Qp2, self.Ze = self.bp6_Q(drvTmp)
            self.v, self.vp, self.vp2  = self.Q / drvTmp.Sd, self.Qp / self.Sp, self.Qp2 / self.Sp2
        elif self.config == "bandpass_pr":
            self.Q, self.Qpr, self.Ze = self.bp4_pr_Q(drvTmp)
            self.v, self.vpr = self.Q / drvTmp.Sd, self.Qpr / self.Sd
        return None


    def plotZe(self):
        """
        Plot the electrical impedance ZeTot in both modulus and phase of driver in enclosure.

        Returns
        -------
        Matplotlib figure

        """
        try:
            fig, ax = plt.subplots()
            ax.semilogx(self.frequencyRange, np.abs(self.Ze),
                        label="modulus")
            ax.legend(loc='upper left')
            ax.grid(linestyle='dotted', which='both')
            ax.set(ylabel="|Z| [Ohm]", xlabel="Frequency [Hz]")
            ax2 = ax.twinx()
            ax2.semilogx(self.frequencyRange, np.angle(self.Ze), '--',
                         label="phase")
            ax2.legend(loc='upper right')
            ax2.set(ylabel="phase [rads]", ylim=[-pi / 2, pi / 2])
            ax2.grid()
            plt.tight_layout()
        except:
            print('No driver set in specified enclosure.')
        return plt.show()

    def plotXVA(self):
        """
        Plot displacement, velocity and acceleration response of driver in enclosure

        Returns
        -------
        Matplotlib figure
        """
        w = self.frequencyRange * 2 * pi
        try:
            fig, ax = plt.subplots(3, 1)
            ax[0].semilogx(self.frequencyRange, np.abs(self.v/1j/w * 1e3), label='Displacement')
            ax[1].semilogx(self.frequencyRange, np.abs(self.v), label='Velocity')
            ax[2].semilogx(self.frequencyRange, np.abs(self.v*1j*w), label='Acceleration')
            ax[2].set(xlabel="Frequency [Hz]")
            ax[0].set(ylabel="mm")
            ax[1].set(ylabel="m/s")
            ax[2].set(ylabel="m/s^2")
            for i in range(3):
                ax[i].grid(which='both')
                ax[i].legend(loc='best')
            plt.tight_layout()
        except:
            print('No driver set in specified enclosure.')
        return None




## Previous functions -- legacy
# def sealedBoxQ(driver, Vb, eta=0, velocity=False, impedance=False, c=air.c, rho=air.rho):
#     """
#     Compute the volume velocity at the driver output (sealed enclosure). Takes Vb as the enclosure volume if it is float
#      or int. Takes Vb as the acoustical impedance if it is passed as a ndarray (e.g. Zab from FEM calculations)
#
#     Parameters
#     ----------
#     driver : electro_acoustic_driver class.
#     Vb : float, int, complex
#         if float or int: Volume of the enclosure. m**3
#         if complex: acoustical impedance seen by the driver
#     eta: float, optional
#         acoustical loss factor in enclosure. Default is 0
#     c : TYPE, optional
#         celerity of sound in air. The default is self.c.
#     rho : TYPE, optional
#         density of air. The default is self.c.
#
#     Returns
#     -------
#     Qs : Volume velocity at the driver.
#     if impedance if True:
#         Ze : Total electrical impedance of the speaker  ()
#     """
#     # data from electro_acoustic_driver class
#
#     if driver.identifier == 'EAC':
#         freq = driver.f_array
#         w = 2 * pi * freq
#         Zac = driver.Zac
#         Zas = driver.Zas
#         Zaf = driver.Zaf
#         Ps = driver.Ps
#         if type(Vb) == float or type(Vb) == int or type(Vb) == np.float64:
#             Cab = Vb/rho/c**2  # compliance of the enclosure
#             if eta == 0:
#                 Zab = 1/(1j*w*Cab)
#             else:
#                 Rab = rho * c / eta / Vb
#                 Zab = 1/((1j*w*Cab) + 1/Rab)
#
#         elif type(Vb) == np.ndarray:
#             Zab = Vb
#
#         # print('Zab = ', Zab)
#         # compute volume velocities
#         Qs = Ps / (Zac+Zas+Zaf+Zab)
#
#         # Total impedance
#         Ze = driver.Ze+driver.Bl**2 /(driver.Zms+driver.Sd**2*(Zab+driver.Zaf))
#
#         # output variables
#         out = []
#         if velocity is False:
#             if impedance is False:
#                 out = Qs
#             else:
#                 out.append(Qs)
#                 out.append(Ze)
#         elif velocity is True:
#             if impedance is False:
#                 out = Qs / driver.Sd
#             else:
#                 out.append(Qs/driver.Sd)
#                 out.append(Ze)
#
#     elif driver.identifier == "PLV":
#         raise TypeError('If enclosure is sealed, no need to get LEM data (in-enclosure acceleration data).')
#     return out
#
#
# def portedBoxQ(driver, Vb, Lp, Sp, velocity=False, impedance=False, c=air.c, rho=air.rho, eta=1e-5):
#     """
#     Compute the volume velocities at the driver and port of a bass-reflex
#     enclosure. Will also output the total impedance expressed in the electrical
#     domain.
#
#     Parameters
#     ----------
#     driver : electro_acoustic_driver class.
#     Vb : Volume of the enclosure. m**3
#     Lp : Length of the port. m
#     Sp : Cross-section area of the port. m**2
#     c : TYPE, optional
#         celerity of sound in air. The default is self.c.
#     rho : TYPE, optional
#         density of air. The default is self.c.
#     eta : FLOAT, optional
#         losse factor inside enclosure
#     Returns
#     -------
#     Qs : Volume velocity at the driver.
#     Qv : Volume velocity at the port.
#     if impedance is True:
#         Ze : Total electrical impedance with acoustical impedance
#
#     """
#
#     if driver.identifier == 'EAC':
#         # data from electro_acoustic_driver class
#         freq = driver.f_array
#         w = 2 * pi * freq
#         k = w / c
#         s = 1j * w
#
#         # box impedance
#         Cab = Vb / rho / c ** 2  # compliance of the enclosure
#         Rab = rho * c / eta / Vb
#         Zbox = parallel(1 / s / Cab, Rab)
#
#         # port impedance
#         rp = np.sqrt(Sp / pi)  # radius of the port
#         Map = rho * Lp / Sp  # acoustical mass of the port
#         Mal = 0.85 * 2 * rp  # length correction
#         # Mal = 0.54*rho/2/rp
#         Map += Mal
#
#         # radiation impedance (port)
#         Pp = 2 * np.pi * rp  # perimeter of the port
#         alpha = np.sqrt(freq) * (0.95e-5 + 2.03e-5) * Pp / 2 / Sp
#         kl = k * (1 + alpha * (1 - 1j))
#         d0 = (0.6133 + 0.85) * rp  # length correction
#         Zrad = rho * c / Sp * (1 / 4 * (kl * rp) ** 2 + 1j * kl * d0)
#         Zp = s * Map + Zrad  # total impedance of the port (tube + radiation)
#
#         # enclosure impedance
#         Zab = parallel(1 / s / Cab, Rab, s * Map + Zrad)
#
#         # total system impedance calculation
#         ZaTot = driver.Zs
#         Ps = driver.Ps
#
#         # compute volume velocities
#         Qs = Ps / (ZaTot + Zab)
#         Qv = -Qs * Zbox / (Zbox + Zp)
#
#         # total electrical impedance
#         Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * (Zab + driver.Zaf))
#
#         # output variables
#         out = []
#         if velocity is False:
#             if impedance is False:
#                 out.append(Qs)
#                 out.append(Qv)
#             else:
#                 out.append(Qs)
#                 out.append(Qv)
#                 out.append(Ze)
#         elif velocity is True:
#             if impedance is False:
#                 out.append(Qs / driver.Sd)
#                 out.append(Qv / Sp)
#             else:
#                 out.append(Qs / driver.Sd)
#                 out.append(Qv / Sp)
#                 out.append(Ze)
#
#     # acceleration data from polytech / Klippel(?) laser
#     elif driver.identifier == 'PLV':
#         freq = driver.freq_dec
#         w = 2 * pi * freq
#         k = w / c
#         s = 1j * w
#
#         # box impedance
#         Cab = Vb / rho / c ** 2  # compliance of the enclosure
#         Rab = rho * c / eta / Vb
#         Zbox = parallel(1 / s / Cab, Rab)
#
#         # port impedance
#         rp = np.sqrt(Sp / pi)  # radius of the port
#         Map = rho * Lp / Sp  # acoustical mass of the port
#         Mal = 0.85 * 2 * rp  # length correction
#         # Mal = 0.54*rho/2/rp
#         Map += Mal
#
#         # radiation impedance (port)
#         Pp = 2 * np.pi * rp  # perimeter of the port
#         alpha = np.sqrt(freq) * (0.95e-5 + 2.03e-5) * Pp / 2 / Sp
#         kl = k * (1 + alpha * (1 - 1j))
#         d0 = (0.6133 + 0.85) * rp  # length correction
#         Zrad = rho * c / Sp * (1 / 4 * (kl * rp) ** 2 + 1j * kl * d0)
#         Zp = s * Map + Zrad  # total impedance of the port (tube + radiation)
#
#         vs = driver.v_mean
#         Sd = ccpcsa(driver.point_cloud)
#         Qs = vs * Sd
#         Qv = -Qs * Zbox / (Zbox + Zp)
#         Ze = np.zeros(len(freq))
#         out = []
#         if velocity is False:
#             if impedance is False:
#                 out.append(Qv)
#             else:
#                 out.append(Qv)
#                 out.append(Ze)
#         elif velocity is True:
#             if impedance is False:
#                 out.append(Qv / Sp)
#             else:
#                 out.append(Qv / Sp)
#                 out.append(Ze)
#     return out
