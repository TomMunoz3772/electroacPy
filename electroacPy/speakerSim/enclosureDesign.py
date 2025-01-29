#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:56:03 2023

@author: tom.munoz
"""

from electroacPy.global_ import air
from electroacPy.speakerSim.electroAcousticDriver import electroAcousticDriver as ead
from electroacPy.general.freqop import freq_log10, laplace
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
    def __init__(self, Vb, frequencyRange=freq_log10(20, 2500, 50),
                 Qab=120, Qal=30, c=air.c, rho=air.rho, **kwargs):
        """
        Setup a louspeaker enclosure as a lumped-element object.

        :param Vb: float,
            Back volume (set behind the drive unit)
        :param frequencyRange: array,

        :param Qab: float,
            Quality factor linked to losses in enclosure (damping). Default is 120
        :param Qal: float,
            Quality factor linked to losses out of enclosure (leakage). Default is 30
        :param c: float,
            Speed of sound in air. Default is 343 m/s
        :param rho: float,
            Air density. Default is 1.22 kg/m^3
        """
        # simulation parameters
        self.c = c
        self.rho = rho
        self.Qab = Qab
        self.Qal = Qal        
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
            self.flange = "single"
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
                elif key == "flange":
                    self.flange = value

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
            self.flange = "single"
            for key, value in kwargs.items():
                setattr(self, key, value)
                if key == "Sp":
                    self.rp = np.sqrt(value / np.pi)
                elif key == "rp":
                    self.Sp = np.pi * value**2
                elif key == "flange":
                    self.flange = value

        ## 4th order bandpass with passive radiator
        elif ("Mmd" in kwargs and "Cmd" in kwargs and
              "Rmd" in kwargs and "Vf" in kwargs and "Sd" in kwargs):
            self.config = "bandpass_pr"
            for key, value in kwargs.items():
                setattr(self, key, value)

        ## ported enclosure
        elif "Lp" in kwargs and ("Sp" in kwargs or "rp" in kwargs):
            self.config = "vented"
            self.flange = "single"
            for key, value in kwargs.items():
                setattr(self, key, value)
                if key == "Sp":
                    self.rp = np.sqrt(value / np.pi)
                elif key == "rp":
                    self.Sp = np.pi * value**2
                elif key == "flange":
                    self.flange = value

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

    ### ======================================
    ## IMPEDANCE AND VOLUME VELOCITY FUNCTIONS
    def sealed_box(self, driver):   
        # prepare some values for components
        Cab = self.Vb / self.rho / self.c**2    # box volume
        wc  = 2 * np.pi * driver.Fs * np.sqrt(1 + driver.Vas / self.Vb)
        Rab = 1 / wc / Cab / self.Qab    # internal losses  
        Ral = self.Qal / wc / Cab        # leakage 
        
        # define components
        enclosure = circuit(self.frequencyRange)
        U   = compe.voltageSource(1, 0, driver.U)
        DRV = EAD(1, 0, 2, 3, driver.Le, driver.Re, 
                  driver.Cms, driver.Mms, driver.Rms, 
                  driver.Bl, driver.Sd, v_probe="v") 
        RAL = compe.resistance(3, 0, Ral)
        RAB = compe.resistance(3, 4, Rab)
        CAB = compe.capacitance(4, 0, Cab)
        RAD = compa.radiator(2, 0, driver.Sd, self.rho, self.c)
        
        # setup and run
        enclosure.addComponent(U, RAL, RAB, CAB, RAD)
        enclosure.addBlock(DRV)
        enclosure.run()
        
        # extract data
        Q  = enclosure.getPotential(2) * RAD.Gs
        v  = enclosure.getFlow("v")
        Ze = -enclosure.getPotential(1) / enclosure.getFlow(1)
        return Q, v, Ze

    def vented_box(self, driver):  
        # prepare some values for components
        Cab = self.Vb / self.rho / self.c**2    # box volume
        wc  = 2 * np.pi * driver.Fs * np.sqrt(1 + driver.Vas / self.Vb)
        Rab = 1 / wc / Cab / self.Qab    # internal losses  
        Ral = self.Qal / wc / Cab        # leakage
        
        # define component
        enclosure = circuit(self.frequencyRange)
        U    = compe.voltageSource(1, 0, driver.U)
        DRV  = EAD(1, 0, 2, 3, driver.Le, driver.Re, 
                   driver.Cms, driver.Mms, driver.Rms, 
                   driver.Bl, driver.Sd, v_probe="v") 
        
        RAL = compe.resistance(3, 0, Ral)
        RAB = compe.resistance(3, 4, Rab)
        CAB = compe.capacitance(4, 0, Cab)

        RAD  = compa.radiator(2, 0, driver.Sd, self.rho, self.c)
        PORT = compa.port(3, 5, self.Lp, self.rp, 
                          self.flange, rho=self.rho, c=self.c)
        RADP = compa.radiator(5, 0, self.Sp, self.rho, self.c)
        
        # setup and run
        enclosure.addComponent(U, RAL, RAB, CAB, RAD, PORT, RADP)
        enclosure.addBlock(DRV)
        enclosure.run()
        
        # extract data
        Q  = enclosure.getPotential(2) * RAD.Gs
        Qp = enclosure.getPotential(5) * RADP.Gs
        v  = enclosure.getFlow("v")
        vp =  Qp / self.Sp
        Ze = -enclosure.getPotential(1) / enclosure.getFlow(1)
        return Q, Qp, v, vp, Ze

    def passive_radiator(self, driver):
        # prepare some values for components
        Cab = self.Vb / self.rho / self.c**2    # box volume
        wc  = 2 * np.pi * driver.Fs * np.sqrt(1 + driver.Vas / self.Vb)
        Rab = 1 / wc / Cab / self.Qab    # internal losses  
        Ral = self.Qal / wc / Cab 
        
        # define components
        enclosure = circuit(self.frequencyRange)
        U     = compe.voltageSource(1, 0, driver.U)
        DRV   = EAD(1, 0, 2, 3, driver.Le, driver.Re, 
                    driver.Cms, driver.Mms, driver.Rms, 
                    driver.Bl, driver.Sd, v_probe="v") 

        # box
        RAL = compe.resistance(3, 0, Ral)        
        RAB = compe.resistance(3, 4, Rab)
        CAB = compe.capacitance(4, 0, Cab)
        
        RAD   = compa.radiator(2, 0, driver.Sd, self.rho, self.c)
        PR    = compa.membrane(3, 5, self.Cmd, self.Mmd, self.Rmd, self.Sd, self.rho, self.c)
        RADPR = compa.radiator(5, 0, self.Sd, self.rho, self.c)       
        
        # setup and run
        enclosure.addComponent(U, RAL, RAB, CAB, RAD, PR, RADPR)
        enclosure.addBlock(DRV)
        enclosure.run()
        
        # extract data
        Q   = enclosure.getPotential(2) * RAD.Gs
        Qpr = enclosure.getPotential(5) * RADPR.Gs
        v   = enclosure.getFlow("v")
        vpr =  Qpr / self.Sd
        Ze  = -enclosure.getPotential(1) / enclosure.getFlow(1)
        return Q, Qpr, v, vpr, Ze

    def bandpass4_port(self, driver):
        # for simplicity, we consider same losses for front and back enclosure
        # prepare some values for components
        Cab = self.Vb / self.rho / self.c**2    # box volume
        wc  = 2 * np.pi * driver.Fs * np.sqrt(1 + driver.Vas / self.Vb)
        Rab = 1 / wc / Cab / self.Qab    # internal losses  
        Ral = self.Qal / wc / Cab        # leakage 
        
        Cabf = self.Vf / self.rho / self.c**2
        wcf  = 2 * np.pi * driver.Fs * np.sqrt(1 + driver.Vas / self.Vf)
        Rabf = 1 / wcf / Cabf / self.Qab    # internal losses (front)
        Ralf = self.Qal / wcf / Cabf        # leakage (front)

        # define component
        enclosure = circuit(self.frequencyRange)
        U     = compe.voltageSource(1, 0, driver.U)
        DRV   = EAD(1, 0, 2, 5, driver.Le, driver.Re, 
                    driver.Cms, driver.Mms, driver.Rms, 
                    driver.Bl, driver.Sd, v_probe="v") 
    
        RALF = compe.resistance(2, 0, Ralf)
        RABF = compe.resistance(2, 3, Rabf)
        CABF = compe.capacitance(3, 0, Cabf)
        PORTF = compa.port(2, 4, self.Lp, self.rp, 
                           self.flange, rho=self.rho, c=self.c)
        RADPF = compa.radiator(4, 0, self.Sp, self.rho, self.c)
        
        RAL = compe.resistance(5, 0, Ral)        
        RAB = compe.resistance(5, 6, Rab)
        CAB = compe.capacitance(6, 0, Cab)    
    
        # setup and run
        enclosure.addComponent(U, RALF, RABF, CABF, 
                       RAL, RAB, CAB, PORTF, RADPF)
        enclosure.addBlock(DRV)
        enclosure.run()
        
        # extract data
        Qp = enclosure.getPotential(4) * RADPF.Gs
        vp = Qp / self.Sp
        v  = enclosure.getFlow("v")
        Ze = -enclosure.getPotential(1) / enclosure.getFlow(1)
        return Qp, vp, v, Ze

    def bandpass6_port(self, driver):
        # for simplicity, we consider same losses for front and back enclosure
        # prepare some values for components
        Cab = self.Vb / self.rho / self.c**2    # box volume
        wc  = 2 * np.pi * driver.Fs * np.sqrt(1 + driver.Vas / self.Vb)
        Rab = 1 / wc / Cab / self.Qab    # internal losses  
        Ral = self.Qal / wc / Cab        # leakage 
        
        Cabf = self.Vf / self.rho / self.c**2
        wcf  = 2 * np.pi * driver.Fs * np.sqrt(1 + driver.Vas / self.Vf)
        Rabf = 1 / wcf / Cabf / self.Qab    # internal losses (front)
        Ralf = self.Qal / wcf / Cabf        # leakage (front)
        
        # define component
        enclosure = circuit(self.frequencyRange)
        U     = compe.voltageSource(1, 0, driver.U)
        DRV   = EAD(1, 0, 2, 5, driver.Le, driver.Re, 
                    driver.Cms, driver.Mms, driver.Rms, 
                    driver.Bl, driver.Sd, v_probe="v") 

        RALF  = compe.resistance(2, 0, Ralf)
        RABF  = compe.resistance(2, 3, Rabf)
        CABF  = compe.capacitance(3, 0, Cabf)
        PORTF = compa.port(2, 4, self.Lp, self.rp, 
                           self.flange, rho=self.rho, c=self.c)
        RADPF = compa.radiator(4, 0, self.Sp, self.rho, self.c)
        
        RAL   = compe.resistance(5, 0, Ral)
        RAB   = compe.resistance(5, 6, Rab)
        CAB   = compe.capacitance(6, 0, Cab)
        PORTB = compa.port(5, 7, self.Lp2, self.rp2, 
                           self.flange, rho=self.rho, c=self.c)
        RADPB = compa.radiator(7, 0, self.Sp2, self.rho, self.c)
        
        # setup and run
        enclosure.addComponent(U, RALF, RABF, CABF,
                       RAL, RAB, CAB, PORTF, RADPF, 
                       PORTB, RADPB)
        enclosure.addBlock(DRV)
        enclosure.run()
        

        # extract data
        Qp = enclosure.getPotential(4) * RADPF.Gs  # front 
        vp = Qp / self.Sp
        Qp2 = enclosure.getPotential(7) * RADPB.Gs # back
        vp2 = Qp2 / self.Sp2
        v  = enclosure.getFlow("v")
        Ze = -enclosure.getPotential(1) / enclosure.getFlow(1)
        return Qp, vp, Qp2, vp2, v, Ze


    def bandpass4_passive_radiator(self, driver):
        # for simplicity, we consider same losses for front and back enclosure
        # prepare some values for components
        Cab = self.Vb / self.rho / self.c**2    # box volume
        wc  = 2 * np.pi * driver.Fs * np.sqrt(1 + driver.Vas / self.Vb)
        Rab = 1 / wc / Cab / self.Qab    # internal losses  
        Ral = self.Qal / wc / Cab        # leakage 
        
        Cabf = self.Vf / self.rho / self.c**2
        wcf  = 2 * np.pi * driver.Fs * np.sqrt(1 + driver.Vas / self.Vf)
        Rabf = 1 / wcf / Cabf / self.Qab    # internal losses (front)
        Ralf = self.Qal / wcf / Cabf        # leakage (front)
        
        # define component
        enclosure = circuit(self.frequencyRange)
        U     = compe.voltageSource(1, 0, driver.U)
        DRV   = EAD(1, 0, 2, 5, driver.Le, driver.Re, 
                    driver.Cms, driver.Mms, driver.Rms, 
                    driver.Bl, driver.Sd, v_probe="v") 

        RALF  = compe.resistance(2, 0, Ralf)
        RABF  = compe.resistance(2, 3, Rabf)
        CABF  = compe.capacitance(3, 0, Cabf)
        PRF   = compa.membrane(2, 4, self.Cmd, self.Mmd, self.Rmd, self.Sd, self.rho, self.c)
        RADPF = compa.radiator(4, 0, self.Sd, self.rho, self.c)
        
        RAL = compe.resistance(5, 0, Ral)        
        RAB = compe.resistance(5, 6, Rab)
        CAB = compe.capacitance(6, 0, Cab)
    
        # setup and run
        enclosure.addComponent(U, RALF, RABF, CABF, PRF, RADPF,
                               RAL, RAB, CAB)
        enclosure.addBlock(DRV)
        enclosure.run()
        
        # extract data
        Qpr = enclosure.getPotential(4) * RADPF.Gs
        vpr = Qpr / self.Sd
        v   = enclosure.getFlow("v")
        Ze  = -enclosure.getPotential(1) / enclosure.getFlow(1)
        return Qpr, vpr, v, Ze
    
    
    def bandpass6_passive_radiator(self, driver):     
        # for simplicity, we consider same losses for front and back enclosure
        # prepare some values for components
        Cab = self.Vb / self.rho / self.c**2    # box volume
        wc  = 2 * np.pi * driver.Fs * np.sqrt(1 + driver.Vas / self.Vb)
        Rab = 1 / wc / Cab / self.Qab    # internal losses  
        Ral = self.Qal / wc / Cab        # leakage 
        
        Cabf = self.Vf / self.rho / self.c**2
        wcf  = 2 * np.pi * driver.Fs * np.sqrt(1 + driver.Vas / self.Vf)
        Rabf = 1 / wcf / Cabf / self.Qab    # internal losses (front)
        Ralf = self.Qal / wcf / Cabf        # leakage (front)
        
        # define component
        enclosure = circuit(self.frequencyRange)
        U     = compe.voltageSource(1, 0, driver.U)
        DRV   = EAD(1, 0, 2, 4, driver.Le, driver.Re, 
                    driver.Cms, driver.Mms, driver.Rms, 
                    driver.Bl, driver.Sd, v_probe="v") 
        
        RALF  = compe.resistance(2, 0, Ralf)
        RABF  = compe.resistance(2, 3, Rabf)
        CABF  = compe.capacitance(3, 0, Cabf)
        PRF   = compa.membrane(2, 4, self.Cmd, self.Mmd, self.Rmd, self.Sd, self.rho, self.c)
        RADPF = compa.radiator(4, 0, self.Sd, self.rho, self.c)
        
        RAL   = compe.resistance(5, 0, Ral)
        RAB   = compe.resistance(5, 6, Rab)
        CAB   = compe.capacitance(6, 0, Cab)
        PRB   = compa.membrane(5, 7, self.Cmd2, self.Mmd2, self.Rmd2, self.Sd2, self.rho, self.c) 
        RADPB = compa.radiator(7, 0, self.Sd2, self.rho, self.c)
        
        
        # setup and run
        # enclosure.addComponent(U, BOXF, PRF, RADPF, BOXB, PRB, RADPB)
        enclosure.addComponent(U, RALF, RABF, CABF, PRF, RADPF,
                       RAL, RAB, CAB, PRB, RADPB)
        enclosure.addBlock(DRV)
        enclosure.run()
        
        # extract data
        Qpr  = enclosure.getPotential(4) * RADPF.Gs
        vpr  = Qpr / self.Sd
        Qpr2 = enclosure.getPotential(7) * RADPB.Gs
        vpr2 = Qpr2 / self.Sd2
        v    = enclosure.getFlow("v")
        Ze   = -enclosure.getPotential(1) / enclosure.getFlow(1)
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


    def plotZe(self, **kwargs):
        """
        Plot the electrical impedance ZeTot in both modulus and phase.
    
        Returns
        -------
        None
    
        """
        
        if "figsize" in kwargs:
            size=kwargs["figsize"]
        else:
            size=None
        
        fig, ax = plt.subplots(2, 1, figsize=size)
        ax[0].semilogx(self.frequencyRange, np.abs(self.Ze))
        ax[0].set(ylabel="Magnitude [Ohm]")
        
        ax[1].semilogx(self.frequencyRange, np.angle(self.Ze))
        ax[1].set(xlabel="Frequency [Hz]", ylabel="Phase [rad]")
        for i in range(2):
            ax[i].grid(which="both", linestyle="dotted")
        plt.tight_layout()
        
        if "savefig" in kwargs:
            path = kwargs["savefig"]
            plt.savefig(path)
      
        return plt.show()
    
    def plotXVA(self, **kwargs):
        """
        Plot the displacement, velocity, and acceleration frequency responses.
    
        Returns
        -------
        None
    
        """
        
        if "figsize" in kwargs:
            size=kwargs["figsize"]
        else:
            size=None
        
        
        # convert velocity into displacement and acceleration
        x = self.v / laplace(self.frequencyRange)
        a = self.v * laplace(self.frequencyRange)
        
        xp   = self.vp / laplace(self.frequencyRange)
        ap   = self.vp * laplace(self.frequencyRange)
        xp2  = self.vp2 / laplace(self.frequencyRange)
        ap2  = self.vp2 * laplace(self.frequencyRange)
        xpr  = self.vpr / laplace(self.frequencyRange)
        apr  = self.vpr * laplace(self.frequencyRange)
        xpr2 = self.vpr2 / laplace(self.frequencyRange)
        apr2 = self.vpr2 * laplace(self.frequencyRange)
        
        fig, ax = plt.subplots(3, 1, figsize=size)
        ax[0].semilogx(self.frequencyRange, np.abs(x), label='Displacement')
        ax[1].semilogx(self.frequencyRange, np.abs(self.v), label='Velocity')
        ax[2].semilogx(self.frequencyRange, np.abs(a), label='Acceleration')
        ax[2].set(xlabel="Frequency [Hz]")
        ax[0].set(ylabel="mm", title="Driver")
        ax[1].set(ylabel="m/s")
        ax[2].set(ylabel="m/s^2", xlabel="Frequency [Hz]")
        for i in range(3):
            ax[i].grid(which='both', linestyle="dotted")
            ax[i].legend(loc='best')
        plt.tight_layout()
        
        if "savefig" in kwargs:
            path = kwargs["savefig"]
            plt.savefig(path)
      
        
        if np.all(self.vp != 0):
            fig, ax = plt.subplots(3, 1, figsize=size)
            ax[0].semilogx(self.frequencyRange, np.abs(xp), label='Displacement')
            ax[1].semilogx(self.frequencyRange, np.abs(self.vp), label='Velocity')
            ax[2].semilogx(self.frequencyRange, np.abs(ap), label='Acceleration')
            ax[0].set(ylabel="mm", title="Port")
            ax[1].set(ylabel="m/s")
            ax[2].set(ylabel="m/s^2", xlabel="Frequency [Hz]")
            for i in range(3):
                ax[i].grid(which='both', linestyle="dotted")
                ax[i].legend(loc='best')
            plt.tight_layout()
            if "savefig" in kwargs:
                path = kwargs["savefig"]
                extension = path[-4:]
                plt.savefig(path[:-4] + "_vp" + extension)       
            
        if np.all(self.vp2 != 0):
            fig, ax = plt.subplots(3, 1, figsize=size)
            ax[0].semilogx(self.frequencyRange, np.abs(xp2), label='Displacement')
            ax[1].semilogx(self.frequencyRange, np.abs(self.vp2), label='Velocity')
            ax[2].semilogx(self.frequencyRange, np.abs(ap2), label='Acceleration')
            ax[0].set(ylabel="mm", title="Port - 2")
            ax[1].set(ylabel="m/s")
            ax[2].set(ylabel="m/s^2", xlabel="Frequency [Hz]")
            for i in range(3):
                ax[i].grid(which='both', linestyle="dotted")
                ax[i].legend(loc='best')
            plt.tight_layout()
            if "savefig" in kwargs:
                path = kwargs["savefig"]
                extension = path[-4:]
                plt.savefig(path[:-4] + "_vp2" + extension)
            
        if np.all(self.vpr != 0):
            fig, ax = plt.subplots(3, 1, figsize=size)
            ax[0].semilogx(self.frequencyRange, np.abs(xpr), label='Displacement')
            ax[1].semilogx(self.frequencyRange, np.abs(self.vpr), label='Velocity')
            ax[2].semilogx(self.frequencyRange, np.abs(apr), label='Acceleration')
            ax[0].set(ylabel="mm", title="Passive-radiator")
            ax[1].set(ylabel="m/s")
            ax[2].set(ylabel="m/s^2", xlabel="Frequency [Hz]")
            for i in range(3):
                ax[i].grid(which='both', linestyle="dotted")
                ax[i].legend(loc='best')  
            plt.tight_layout()
            if "savefig" in kwargs:
                path = kwargs["savefig"]
                extension = path[-4:]
                plt.savefig(path[:-4] + "_vpr" + extension)
            
        if np.all(self.vpr2 != 0):
            fig, ax = plt.subplots(3, 1, figsize=size)
            ax[0].semilogx(self.frequencyRange, np.abs(xpr2), label='Displacement')
            ax[1].semilogx(self.frequencyRange, np.abs(self.vpr2), label='Velocity')
            ax[2].semilogx(self.frequencyRange, np.abs(apr2), label='Acceleration')
            ax[0].set(ylabel="mm", title="Passive-radiator - 2")
            ax[1].set(ylabel="m/s")
            ax[2].set(ylabel="m/s^2", xlabel="Frequency [Hz]")
            for i in range(3):
                ax[i].grid(which='both', linestyle="dotted")
                ax[i].legend(loc='best')   
            plt.tight_layout()
            if "savefig" in kwargs:
                path = kwargs["savefig"]
                extension = path[-4:]
                plt.savefig(path[:-4] + "_vpr2" + extension)
        
        return plt.show()
    
    def exportZe(self, folder_name, file_name):
        import os
        
        # Create the folder if it doesn't exist
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        
        module = np.abs(self.Ze)
        phase = np.angle(self.Ze, deg=True)
        path = os.path.join(folder_name, file_name)
        np.savetxt(path, np.array([self.frequencyRange, module, phase]).T,
                   fmt="%.3f",
                   header="Freq[Hz]  Imp[Ohm]  Phase[Deg]",
                   delimiter=',',
                   comments='')
        
        
        
        
        
        
