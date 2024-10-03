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
import matplotlib.pyplot as plt
import numpy as np

# circuit solver
import electroacPy.circuitSolver.components.electric as compe
import electroacPy.circuitSolver.components.acoustic as compa
from electroacPy.circuitSolver.blocks.electrodynamic import EAD
from electroacPy.circuitSolver.solver import circuit

pi = np.pi

## New approach
class speakerBox:
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
        # self.computeImpedance()

        ## TO STORE RESULTS
        # velocity
        self.v    = np.zeros(len(frequencyRange), dtype=complex)
        self.vp   = np.zeros(len(frequencyRange), dtype=complex)
        self.vp2  = np.zeros(len(frequencyRange), dtype=complex)
        self.vpr  = np.zeros(len(frequencyRange), dtype=complex)
        self.vpr2 = np.zeros(len(frequencyRange), dtype=complex)
        # volume velocity
        self.Q    = np.zeros(len(frequencyRange), dtype=complex)
        self.Qp   = np.zeros(len(frequencyRange), dtype=complex)
        self.Qp2  = np.zeros(len(frequencyRange), dtype=complex)
        self.Qpr  = np.zeros(len(frequencyRange), dtype=complex)
        self.Qpr2 = np.zeros(len(frequencyRange), dtype=complex)
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

    # def sealed_box(self,):
    #     """
    #     Create a sealed enclosure.
    #     :return:
    #     """
    #     # Cab = self.Vb / self.rho / self.c ** 2          # compliance of the enclosure
    #     # Rab = self.rho * self.c / self.eta / self.Vb    # losses through leakage
    #     # Zab = parallel(1 / self.s / Cab, Rab)
        
    #     enclosure = circuit(self.frequencyRange)
    #     Ps  = compa.pressureSource(1, 0, 1)
    #     box = compa.cavity(1, 0, self.Vb)
        
    #     enclosure.addComponent(Ps, box)
    #     enclosure.run()
        
        
    #     return Zab

    # def vented_box(self,):
    #     """
    #     Create ported enclosure.
    #     :return:
    #     """
    #     w = 2*np.pi*self.frequencyRange

    #     # box impedance
    #     Cab = self.Vb / self.rho / self.c ** 2  # compliance of the enclosure
    #     Rab = self.rho * self.c / self.eta / self.Vb
    #     Zab = parallel(1 / self.s / Cab, Rab)

    #     # port impedance
    #     # if flanged is True:
    #     #     ll = 8 * self.rp / 3 / np.pi
    #     # else:
    #     #     ll = 2 * self.rp / np.pi
    #     ll = 2 * self.rp / np.pi
    #     Map = (self.Lp + 0.64 * self.rp) * self.rho / self.Sp
    #     Rap = np.sqrt(2 * w * self.rho * 1.86e-5) / self.Sp * (ll / self.rp + 0.7)

    #     # port radiation
    #     RAR2 = 0.159 * w ** 2 * self.rho / self.c
    #     MAR2 = 0.270 * self.rho / self.rp
    #     Zap = self.s * MAR2 + RAR2 + self.s * Map + Rap

    #     # total enclosure impedance
    #     Za_in = parallel(Zap, Zab)  # total acoustic impedance
    #     return Za_in, Zab, Zap

    # def pr_box(self,):
    #     """
    #     Create a passive radiator enclosure.
    #     :return:
    #     """
    #     w = 2 * np.pi * self.frequencyRange

    #     # box impedance
    #     Cab = self.Vb / self.rho / self.c ** 2  # compliance of the enclosure
    #     Rab = self.rho * self.c / self.eta / self.Vb
    #     Zab = parallel(1 / self.s / Cab, Rab)

    #     # passive radiator's impedance
    #     RAR = 0.159 * w ** 2 * self.rho / self.c         # piston's radiation impedance - air resistance
    #     MAR = 0.270 * self.rho / np.sqrt(self.Sd/np.pi)  # piston's radiation impedance - air mass
    #     Zapr = ((self.Rmd + self.s * self.Mmd + 1 / (self.s * self.Cmd)) / self.Sd ** 2) + self.s * MAR + RAR

    #     # total enclosure impedance
    #     Za_in = parallel(Zapr, Zab)  # total acoustic impedance
    #     return Za_in, Zab, Zapr

    # def bp4_box(self,):
    #     """
    #     Create a 4th order bandpass enclosure (ports).
    #     :return:
    #     """

    #     w = 2 * np.pi * self.frequencyRange

    #     # back enclosure impedance
    #     Cab = self.Vb / self.rho / self.c ** 2  # compliance of the enclosure
    #     Rab = self.rho * self.c / self.eta / self.Vb
    #     Zab = parallel(1 / self.s / Cab, Rab)

    #     # front enclosure impedance
    #     Caf = self.Vf / self.rho / self.c**2
    #     Raf = self.rho * self.c / self.eta / self.Vf
    #     Zaf = parallel(1/self.s/Caf, Raf)
    #     ll = 2 * self.rp / np.pi
    #     Map = (self.Lp + 0.64 * self.rp) * self.rho / self.Sp
    #     Rap = np.sqrt(2 * w * self.rho * 1.86e-5) / self.Sp * (ll / self.rp + 0.7)

    #     # port radiation
    #     RAR = 0.159 * w ** 2 * self.rho / self.c
    #     MAR = 0.270 * self.rho / self.rp
    #     Zap = self.s * MAR + RAR + self.s * Map + Rap
    #     Za_front = parallel(Zap, Zaf)
    #     Za_in = parallel(Za_front, Zab)
    #     return Za_in, Zab, Zaf, Zap

    # def bp6_box(self,):
    #     """
    #     Create a 6th order bandpass enclosure (ports).
    #     :return:
    #     """
    #     w = 2 * np.pi * self.frequencyRange

    #     # back enclosure impedance
    #     Cab = self.Vb / self.rho / self.c ** 2  # compliance of the enclosure
    #     Rab = self.rho * self.c / self.eta / self.Vb
    #     Zab = parallel(1 / self.s / Cab, Rab)
    #     ll2 = 2 * self.rp2 / np.pi
    #     Map2 = (self.Lp2 + 0.64 * self.rp2) * self.rho / self.Sp2
    #     Rap2 = np.sqrt(2 * w * self.rho * 1.86e-5) / self.Sp2 * (ll2 / self.rp2 + 0.7)
    #     RAR2 = 0.159 * w ** 2 * self.rho / self.c
    #     MAR2 = 0.270 * self.rho / self.rp2
    #     Zap2 = self.s * MAR2 + RAR2 + self.s * Map2 + Rap2

    #     # front enclosure impedance
    #     Caf = self.Vf / self.rho / self.c**2
    #     Raf = self.rho * self.c / self.eta / self.Vf
    #     Zaf = parallel(1/self.s/Caf, Raf)
    #     ll = 2 * self.rp / np.pi
    #     Map = (self.Lp + 0.64 * self.rp) * self.rho / self.Sp
    #     Rap = np.sqrt(2 * w * self.rho * 1.86e-5) / self.Sp * (ll / self.rp + 0.7)

    #     # port radiation
    #     RAR = 0.159 * w ** 2 * self.rho / self.c
    #     MAR = 0.270 * self.rho / self.rp
    #     Zap = self.s * MAR + RAR + self.s * Map + Rap
    #     Za_front = parallel(Zap, Zaf)
    #     Za_back = parallel(Zap2, Zab)
    #     Za_in = parallel(Za_front, Zab)
    #     return Za_in, Zab, Zaf, Zap, Zap2

    # def bp4_pr_box(self, ):
    #     """
    #     Create a 4th order bandpass enclosure (passive radiator).
    #     :return:
    #     """

    #     w = 2 * np.pi * self.frequencyRange

    #     # box impedance - front
    #     Caf = self.Vf / self.rho / self.c ** 2  # compliance of the enclosure
    #     Raf = self.rho * self.c / self.eta / self.Vf
    #     Zaf = parallel(1 / self.s / Caf, Raf)

    #     # box impedance - back
    #     Cab = self.Vb / self.rho / self.c ** 2  # compliance of the enclosure
    #     Rab = self.rho * self.c / self.eta / self.Vb
    #     Zab = parallel(1 / self.s / Cab, Rab)

    #     # passive radiator's impedance
    #     RAR = 0.159 * w ** 2 * self.rho / self.c         # piston's radiation impedance - air resistance
    #     MAR = 0.270 * self.rho / np.sqrt(self.Sd/np.pi)  # piston's radiation impedance - air mass
    #     Zapr = ((self.Rmd + self.s * self.Mmd + 1 / (self.s * self.Cmd)) / self.Sd ** 2) + self.s * MAR + RAR

    #     # total enclosure impedance
    #     Za_front = parallel(Zapr, Zaf)
    #     Za_in = parallel(Za_front, Zab)
    #     return Za_in, Zab, Zaf, Zapr

    ### ======================================
    ## IMPEDANCE AND VOLUME VELOCITY FUNCTIONS
    def sealed_box(self, driver):
        # Q = driver.Ps / (driver.Zs + self.Za_in)
        # Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * (self.Za_in + driver.Zaf))
        
        # define components
        enclosure = circuit(self.frequencyRange)
        U   = compe.voltageSource(1, 0, driver.U)
        DRV = EAD(1, 0, 2, 3, driver.Le, driver.Re, 
                  driver.Cms, driver.Mms, driver.Rms, 
                  driver.Bl, driver.Sd, v_probe="v") 
        BOX = compa.cavity(3, 0, self.Vb, self.eta, self.rho, self.c)
        RAD = compa.radiator(2, 0, driver.Sd, self.rho, self.c)
        
        # setup and run
        enclosure.addComponent(U, BOX, RAD)
        enclosure.addBlock(DRV)
        enclosure.run()
        
        # extract data
        Q  = enclosure.getPotential(2) * RAD.Gs
        v  = enclosure.getFlow("v")
        Ze = enclosure.getPotential(1) / enclosure.getFlow(1)
        return Q, v, Ze

    def vented_box(self, driver):
        # Q = driver.Ps / (driver.Zs + self.Za_in)
        # Qp = - Q * self.Zab / (self.Zab + self.Zap)
        # Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * (self.Za_in + driver.Zaf))
        
        # define component
        enclosure = circuit(self.frequencyRange)
        U    = compe.voltageSource(1, 0, driver.U)
        DRV  = EAD(1, 0, 2, 3, driver.Le, driver.Re, 
                   driver.Cms, driver.Mms, driver.Rms, 
                   driver.Bl, driver.Sd, v_probe="v") 
        BOX  = compa.cavity(3, 0, self.Vb, self.eta, self.rho, self.c)
        RAD  = compa.radiator(2, 0, driver.Sd, self.rho, self.c)
        PORT = compa.port(3, 4, self.Lp, self.rp, self.rho, self.c)
        RADP = compa.radiator(4, 0, self.Sp, self.rho, self.c)
        
        # setup and run
        enclosure.addComponent(U, BOX, RAD, PORT, RADP)
        enclosure.addBlock(DRV)
        enclosure.run()
        
        # extract data
        Q  = enclosure.getPotential(2) * RAD.Gs
        Qp = enclosure.getPotential(4) * RADP.Gs
        v  = enclosure.getFlow("v")
        vp =  Qp / self.Sp
        Ze = enclosure.getPotential(1) / enclosure.getFlow(1)
        return Q, Qp, v, vp, Ze

    def passive_radiator(self, driver):
        # Q = driver.Ps / (driver.Zs + self.Za_in)
        # Qpr = - Q * self.Zab / (self.Zab + self.Zapr)
        # Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * (self.Za_in + driver.Zaf))
        
        # define components
        enclosure = circuit(self.frequencyRange)
        U     = compe.voltageSource(1, 0, driver.U)
        DRV   = EAD(1, 0, 2, 3, driver.Le, driver.Re, 
                    driver.Cms, driver.Mms, driver.Rms, 
                    driver.Bl, driver.Sd, v_probe="v") 
        BOX   = compa.cavity(3, 0, self.Vb, self.eta, self.rho, self.c)
        RAD   = compa.radiator(2, 0, driver.Sd, self.rho, self.c)
        PR    = compa.membrane(3, 4, self.Cmd, self.Mmd, self.Rmd, self.Sd, self.rho, self.c)
        RADPR = compa.radiator(4, 0, self.Sd, self.rho, self.c)       
        
        # setup and run
        enclosure.addComponent(U, BOX, RAD, PR, RADPR)
        enclosure.addBlock(DRV)
        enclosure.run()
        
        # extract data
        Q   = enclosure.getPotential(2) * RAD.Gs
        Qpr = enclosure.getPotential(4) * RADPR.Gs
        v   = enclosure.getFlow("v")
        vpr =  Qpr / self.Sd
        Ze  = enclosure.getPotential(1) / enclosure.getFlow(1)
        return Q, Qpr, v, vpr, Ze

    def bandpass4_port(self, driver):
    #     Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 *
    #                                        (self.Za_in + driver.Zac + driver.Zas))
    #     # MNA approach -- we solve AX=Z
    #     M = 1 # number of independant voltage source (1 driver)
    #     N = 4 # number of nodes

    #     # matrix initialization (see QUCS' online documentation)
    #     A = np.zeros([M+N, M+N], dtype=complex)
    #     G = np.zeros([N, N], dtype=complex)
    #     Q = np.zeros(len(self.frequencyRange), dtype=complex)
    #     Qp = np.zeros(len(self.frequencyRange), dtype=complex)

    #     # let's solve the system
    #     for i in range(len(self.frequencyRange)):
    #         Zad = driver.Zac[i]+driver.Zas[i]
    #         Zap, Zab, Zaf = self.Zap[i], self.Zab[i], self.Zaf[i]
    #         Ps = driver.Ps[i]
    #         G[0, 0] = 1/Zad
    #         G[1, 1] = 1/Zad + 1/Zaf + 1/Zap
    #         G[2, 2] = 1/Zad
    #         G[3, 3] = 1/Zad + 1/Zab
    #         G[0, 1], G[1, 0] = -1/Zad, -1/Zad
    #         G[2, 3], G[3, 2] = -1/Zad, -1/Zad
    #         B = np.array([[1, 0, -1, 0]]).T
    #         C = B.T
    #         D = 0
    #         Z = np.array([[0, 0, 0, 0, Ps]]).T
    #         A[0:4, 0:4] = G
    #         A[0:4, -1:] = B
    #         A[-1:, 0:4] = C
    #         A[-1, -1] = D
    #         # print(A)
    #         X = np.linalg.inv(A)@Z
    #         Qp[i] = X[1, 0] / Zap
    #         Q[i] = X[-1, 0]
    
        # define component
        enclosure = circuit(self.frequencyRange)
        U     = compe.voltageSource(1, 0, driver.U)
        DRV   = EAD(1, 0, 2, 4, driver.Le, driver.Re, 
                    driver.Cms, driver.Mms, driver.Rms, 
                    driver.Bl, driver.Sd, v_probe="v") 
        BOXF  = compa.cavity(2, 0, self.Vf, self.eta, self.rho, self.c)
        PORTF = compa.port(2, 3, self.Lp, self.rp, self.rho, self.c)
        RADPF = compa.radiator(3, 0, self.Sp, self.rho, self.c)
        BOXB  = compa.cavity(4, 0, self.Vb, self.eta, self.rho, self.c)
    
        # setup and run
        enclosure.addComponent(U, BOXF, PORTF, RADPF, BOXB)
        enclosure.addBlock(DRV)
        enclosure.run()
        
        # extract data
        Qp = enclosure.getPotential(3) * RADPF.Gs
        vp = Qp / self.Sp
        v  = enclosure.getFlow("v")
        Ze = enclosure.getPotential(1) / enclosure.getFlow(1)
        return Qp, vp, v, Ze

    def bandpass6_port(self, driver):
        # Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 *
        #                                    (self.Za_in + driver.Zac + driver.Zas))

        # # MNA approach -- we solve AX=Z
        # M = 1  # number of independant voltage source (1 driver)
        # N = 4  # number of nodes

        # # matrix initialization (see QUCS' online documentation)
        # A = np.zeros([M + N, M + N], dtype=complex)
        # G = np.zeros([N, N], dtype=complex)
        # Q = np.zeros(len(self.frequencyRange), dtype=complex)
        # Qp = np.zeros(len(self.frequencyRange), dtype=complex)
        # Qp2 = np.zeros(len(self.frequencyRange), dtype=complex)

        # # let's solve the system
        # for i in range(len(self.frequencyRange)):
        #     Zad = driver.Zac[i] + driver.Zas[i]
        #     Zap, Zap2, Zab, Zaf = self.Zap[i], self.Zap2[i], self.Zab[i], self.Zaf[i]
        #     Ps = driver.Ps[i]
        #     G[0, 0] = 1 / Zad
        #     G[1, 1] = 1 / Zad + 1 / Zaf + 1 / Zap
        #     G[2, 2] = 1 / Zad
        #     G[3, 3] = 1 / Zad + 1 / Zab + 1 / Zap2
        #     G[0, 1], G[1, 0] = -1 / Zad, -1 / Zad
        #     G[2, 3], G[3, 2] = -1 / Zad, -1 / Zad
        #     B = np.array([[1, 0, -1, 0]]).T
        #     C = B.T
        #     D = 0
        #     Z = np.array([[0, 0, 0, 0, Ps]]).T
        #     A[0:4, 0:4] = G
        #     A[0:4, -1:] = B
        #     A[-1:, 0:4] = C
        #     A[-1, -1] = D

        #     X = np.linalg.inv(A) @ Z
        #     Qp[i] = X[1, 0] / Zap
        #     Qp2[i] = X[3, 0] / Zap2
        #     Q[i] = X[-1, 0]

        # define component
        enclosure = circuit(self.frequencyRange)
        U     = compe.voltageSource(1, 0, driver.U)
        DRV   = EAD(1, 0, 2, 4, driver.Le, driver.Re, 
                    driver.Cms, driver.Mms, driver.Rms, 
                    driver.Bl, driver.Sd, v_probe="v") 
        BOXF  = compa.cavity(2, 0, self.Vf, self.eta, self.rho, self.c)
        PORTF = compa.port(2, 3, self.Lp, self.rp, self.rho, self.c)
        RADPF = compa.radiator(3, 0, self.Sp, self.rho, self.c)
        BOXB  = compa.cavity(4, 0, self.Vb, self.eta, self.rho, self.c)
        PORTB = compa.port(4, 5, self.Lp2, self.rp2, self.rho, self.c)
        RADPB = compa.radiator(5, 0, self.Sp2, self.rho, self.c)
    
        # setup and run
        enclosure.addComponent(U, BOXF, PORTF, RADPF, BOXB, PORTB, RADPB)
        enclosure.addBlock(DRV)
        enclosure.run()
        
        
        # extract data
        Qp = enclosure.getPotential(3) * RADPF.Gs
        vp = Qp / self.Sp
        Qp2 = enclosure.getPotential(5) * RADPB.Gs
        vp2 = Qp2 / self.Sp2
        v  = enclosure.getFlow("v")
        Ze = enclosure.getPotential(1) / enclosure.getFlow(1)
        return Qp, vp, Qp2, vp2, v, Ze


    def bandpass4_passive_radiator(self, driver):
        # # MNA approach -- we solve AX=Z
        # M = 1  # number of independant voltage source (1 driver)
        # N = 4  # number of nodes

        # # matrix initialization (see QUCS' online documentation)
        # A = np.zeros([M + N, M + N], dtype=complex)
        # G = np.zeros([N, N], dtype=complex)
        # Q = np.zeros(len(self.frequencyRange), dtype=complex)
        # Qpr = np.zeros(len(self.frequencyRange), dtype=complex)

        # # let's solve the system
        # for i in range(len(self.frequencyRange)):
        #     Zad = driver.Zac[i] + driver.Zas[i]
        #     Zapr, Zab, Zaf = self.Zapr[i], self.Zab[i], self.Zaf[i]
        #     Ps = driver.Ps[i]
        #     G[0, 0] = 1 / Zad
        #     G[1, 1] = 1 / Zad + 1 / Zaf + 1 / Zapr
        #     G[2, 2] = 1 / Zad
        #     G[3, 3] = 1 / Zad + 1 / Zab
        #     G[0, 1], G[1, 0] = -1 / Zad, -1 / Zad
        #     G[2, 3], G[3, 2] = -1 / Zad, -1 / Zad
        #     B = np.array([[1, 0, -1, 0]]).T
        #     C = B.T
        #     D = 0
        #     Z = np.array([[0, 0, 0, 0, Ps]]).T
        #     A[0:4, 0:4] = G
        #     A[0:4, -1:] = B
        #     A[-1:, 0:4] = C
        #     A[-1, -1] = D
        #     # print(A)
        #     X = np.linalg.inv(A) @ Z
        #     Qpr[i] = X[1, 0] / Zapr
        #     Q[i] = X[-1, 0]

        # Za_in = Ps / Q
        # Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 *
        #                                    (Za_in + driver.Zac + driver.Zas))
        
        # define component
        enclosure = circuit(self.frequencyRange)
        U     = compe.voltageSource(1, 0, driver.U)
        DRV   = EAD(1, 0, 2, 4, driver.Le, driver.Re, 
                    driver.Cms, driver.Mms, driver.Rms, 
                    driver.Bl, driver.Sd, v_probe="v") 
        BOXF  = compa.cavity(2, 0, self.Vf, self.eta, self.rho, self.c)
        PRF   = compa.membrane(2, 3, self.Cmd, self.Mmd, self.Rmd, self.Sd, self.rho, self.c)
        RADPF = compa.radiator(3, 0, self.Sd, self.rho, self.c)
        BOXB  = compa.cavity(4, 0, self.Vb, self.eta, self.rho, self.c)
    
        # setup and run
        enclosure.addComponent(U, BOXF, PRF, RADPF, BOXB)
        enclosure.addBlock(DRV)
        enclosure.run()
        
        # extract data
        Qpr = enclosure.getPotential(3) * RADPF.Gs
        vpr = Qpr / self.Sd
        v   = enclosure.getFlow("v")
        Ze  = enclosure.getPotential(1) / enclosure.getFlow(1)
        return Qpr, vpr, v, Ze
    
    
    def bandpass6_passive_radiator(self, driver):     
        # define component
        enclosure = circuit(self.frequencyRange)
        U     = compe.voltageSource(1, 0, driver.U)
        DRV   = EAD(1, 0, 2, 4, driver.Le, driver.Re, 
                    driver.Cms, driver.Mms, driver.Rms, 
                    driver.Bl, driver.Sd, v_probe="v") 
        BOXF  = compa.cavity(2, 0, self.Vf, self.eta, self.rho, self.c)
        PRF   = compa.membrane(2, 3, self.Cmd, self.Mmd, self.Rmd, self.Sd, self.rho, self.c)
        RADPF = compa.radiator(3, 0, self.Sd, self.rho, self.c)
        BOXB  = compa.cavity(4, 0, self.Vb, self.eta, self.rho, self.c)
        PRB   = compa.membrane(4, 5, self.Cmd2, self.Mmd2, self.Rmd2, self.Sd2, self.rho, self.c) 
        RADPB = compa.radiator(5, 0, self.Sd2, self.rho, self.c)
        
        # setup and run
        enclosure.addComponent(U, BOXF, PRF, RADPF, BOXB, PRB, RADPB)
        enclosure.addBlock(DRV)
        enclosure.run()
        
        # extract data
        Qpr  = enclosure.getPotential(3) * RADPF.Gs
        vpr  = Qpr / self.Sd
        Qpr2 = enclosure.getPotential(5) * RADPB.Gs
        vpr2 = Qpr2 / self.Sd2
        v    = enclosure.getFlow("v")
        Ze   = enclosure.getPotential(1) / enclosure.getFlow(1)
        return Qpr, vpr, Qpr2, vpr2, v, Ze


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
            self.Q, self.v, self.Ze = self.sealed_box(drvTmp)
        elif self.config == "vented":
            self.Q, self.Qp, self.v, self.vp, self.Ze = self.vented_box(drvTmp)
        elif self.config == "passiveRadiator":
            self.Q, self.Qpr, self.v, self.vpr, self.Ze = self.passive_radiator(drvTmp)
        elif self.config == "bandpass":
            self.Qp, self.vp, self.v, self.Ze = self.bandpass4_port(drvTmp)
        elif self.config == "bandpass_2":
            self.Qp, self.vp, self.Qp2, self.vp2, self.v, self.Ze = self.bandpass6_port(drvTmp)
        elif self.config == "bandpass_pr":
            self.Qpr, self.vpr, self.v, self.Ze = self.bandpass4_passive_radiator(drvTmp)
        elif self.config == "bandpass_pr_2":
            self.Qpr, self.vpr, self.Qpr2, self.vpr2, self.v, self.Ze = self.bandpass6_passive_radiator(drvTmp)
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
            ax2.set(ylabel="phase [rads]", ylim=[-pi, pi])
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

