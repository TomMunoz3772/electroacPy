#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 15:50:43 2023

@author: tom.munoz
"""
import numpy as np
import electroacPy.general as gtb
from electroacPy.global_ import air
from electroacPy.general import lp_loaders as lpl
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
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


        # Ref Signals
        # Velocity
        self.Hv = Bl/self.Ze  / (self.Zms + Bl**2 / self.Ze) * self.U

        # Displacement
        self.Hx = self.Hv / s

        # Acceleration
        self.Ha = self.Hv * s

        # Acoustic simulation reference
        self.ref2bem = None

        # in box velocity and impedance
        self.inBox    = False
        self.isPorted = False  # easier to manage if radiator in study_ is not a speakerBox object
        self.v        = self.Hv 
        self.Q        = self.v * self.Sd

        # references
        self.ref2bem   = False
        self.poly_data = False  # is class from polytech?

        
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
        ax[0].semilogx(self.f_array, np.abs(self.ZeTot))
        ax[0].set(ylabel="Magnitude [Ohm]")
        
        ax[1].semilogx(self.f_array, np.angle(self.ZeTot))
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
        
        fig, ax = plt.subplots(3, 1, figsize=size)
        ax[0].semilogx(self.f_array, np.abs(self.Hx*1e3), label='Displacement')
        ax[1].semilogx(self.f_array, np.abs(self.Hv), label='Velocity')
        ax[2].semilogx(self.f_array, np.abs(self.Ha), label='Acceleration')
        ax[2].set(xlabel="Frequency [Hz]")
        ax[0].set(ylabel="mm")
        ax[1].set(ylabel="m/s")
        ax[2].set(ylabel="m/s^2")
        for i in range(3):
            ax[i].grid(which='both', linestyle="dotted")
            ax[i].legend(loc='best')
        plt.tight_layout()
       
        if "savefig" in kwargs:
            path = kwargs["savefig"]
            plt.savefig(path)
      
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
        Compute Volume from Qtc value using Tkinter instead of Matplotlib widgets
    
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
        omega = 2 * np.pi * f_axis
        k = omega / c
    
        Zac = driver.Zac
        Zas = driver.Zas
        Zab = 1 / 1j / omega / Cab
        Ps = driver.Ps
        Qs = Ps / (Zac + Zas + Zab)  # removed Zaf
    
        p = 1j * k * rho * c * Qs * np.exp(-1j * k * 1) / (2 * np.pi * 1)
        Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * (Zab))
    
        # Setup Tkinter window
        root = tk.Tk()
        root.title("Sealed Alignment")
    
        # Create figure (using matplotlib's Figure class, not plt.subplots)
        fig = Figure(figsize=(6, 4))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        
        # Initial plot (SPL and Magnitude)
        ax1.semilogx(f_axis, gtb.gain.SPL(p), label='SPL')
        ax2.semilogx(f_axis, np.abs(Ze), label='Magnitude')
    
        ax2.set(xlabel='Frequency [Hz]', ylabel='Impedance')
        ax1.set(ylabel='SPL [dB] at 1 meter')
        ax1.legend(loc='best')
        ax2.legend(loc='best')
    
        for ax in [ax1, ax2]:
            ax.grid(which='both')
    
        # Embed the plot into Tkinter
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
        # Create Labels and Entry fields
        defaultQtc = round(Qtc, 4)
    
        label_qtc = ttk.Label(root, text="Qtc:")
        label_qtc.pack(side=tk.LEFT, padx=5)
        
        entry_qtc = ttk.Entry(root, width=10)
        entry_qtc.pack(side=tk.LEFT, padx=5)
        entry_qtc.insert(0, str(defaultQtc))
    
        label_vb = ttk.Label(root, text="Vb (L):")
        label_vb.pack(side=tk.LEFT, padx=5)
        
        entry_vb = ttk.Entry(root, width=10)
        entry_vb.pack(side=tk.LEFT, padx=5)
        entry_vb.insert(0, str(round(Vb * 1e3, 4)))
    
        # Function to update the plot based on Qtc entry
        def update_plot():
            try:
                Qtc_value = float(entry_qtc.get())
                fc = Qtc_value / driver.Qts * driver.Fs
                Vb_new = driver.Vas / ((fc / driver.Fs)**2 - 1)
                Cab_new = Vb_new / rho / c**2
                Zab_new = 1 / 1j / omega / Cab_new
    
                Qs_new = Ps / (Zac + Zas + Zab_new)
                p_new = 1j * k * rho * c * Qs_new * np.exp(-1j * k * 1) / (2 * np.pi * 1)
                Ze_new = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd ** 2 * Zab_new)
    
                # Clear and update the plots
                ax1.clear()
                ax2.clear()
    
                ax1.semilogx(f_axis, gtb.gain.SPL(p_new), label='SPL')
                ax2.semilogx(f_axis, np.abs(Ze_new), label='Magnitude')
    
                ax2.set(xlabel='Frequency [Hz]', ylabel='Impedance [Ohm]')
                ax1.set(ylabel='SPL [dB] at 1 meter')
                ax1.legend(loc='best')
                ax2.legend(loc='upper left')
    
                for ax in [ax1, ax2]:
                    ax.grid(which='both')
    
                canvas.draw()
    
                # Update volume value in the entry box
                entry_vb.delete(0, tk.END)
                entry_vb.insert(0, str(round(Vb_new * 1e3, 4)))
    
            except ValueError:
                pass  # Prevent the function from crashing if non-numeric input is entered
    
        # Function to update Qtc based on Volume entry
        def update_qtc():
            try:
                Vb_value = float(entry_vb.get()) * 1e-3  # Convert from liters to cubic meters
                fc_new = driver.Fs * np.sqrt(driver.Vas / Vb_value + 1)
                Qtc_new = fc_new / driver.Fs * driver.Qts
    
                # Update Qtc value in the entry box
                entry_qtc.delete(0, tk.END)
                entry_qtc.insert(0, str(round(Qtc_new, 4)))
    
                update_plot()  # Automatically update plot with new values
    
            except ValueError:
                pass  # Handle non-numeric input
    
        # Bind events to update the plot automatically when the user presses Enter or leaves the entry field
        entry_qtc.bind("<Return>", lambda event: update_plot())
        entry_qtc.bind("<FocusOut>", lambda event: update_plot())
        
        entry_vb.bind("<Return>", lambda event: update_qtc())
        entry_vb.bind("<FocusOut>", lambda event: update_qtc())
    
        root.mainloop()
     
    
    def portedAlignment(self):
        driver = copy(self)
        c = self.c
        rho = self.rho
        f_axis = driver.f_array
        omega = 2 * np.pi * f_axis
        k = omega / c
        s = 1j * omega
        eta = 1e-5
    
        # Default parameters
        self.Vb = copy(driver.Vas)
        self.Lp = 343 / driver.Fs / 100  # Length in meters
        self.rp = self.Lp / 2  # Radius in meters
        self.Sp = np.pi * self.rp ** 2  # Port cross-sectional area
    
        # Create input widgets
        default_volume = str(round(self.Vb * 1e3, 2))
        default_length = str(round(self.Lp * 1e2, 2))
        default_radius = str(round(self.rp * 1e2, 2))
        default_section = str(round(self.Sp * 1e4, 2))
    
        # GUI creation
        root = tk.Tk()
        root.title("Ported Alignment")
    
        # Create a matplotlib figure
        fig = Figure(figsize=(6, 4))
        ax_spl = fig.add_subplot(211)
        ax_imp = fig.add_subplot(212)
        
        
        # Create canvas for the plot and add to the tkinter window
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
        # Create textboxes in a horizontal row at the bottom
        volume_label = ttk.Label(root, text="Vb (L):")
        volume_label.pack(side=tk.LEFT, padx=5)
        volume_entry = ttk.Entry(root, width=10)
        volume_entry.pack(side=tk.LEFT, padx=5)
        volume_entry.insert(0, default_volume)
        
        length_label = ttk.Label(root, text="Lp (cm):")
        length_label.pack(side=tk.LEFT, padx=5)
        length_entry = ttk.Entry(root, width=10)
        length_entry.pack(side=tk.LEFT, padx=5)
        length_entry.insert(0, default_length)
        
        radius_label = ttk.Label(root, text="rp (cm):")
        radius_label.pack(side=tk.LEFT, padx=5)
        radius_entry = ttk.Entry(root, width=10)
        radius_entry.pack(side=tk.LEFT, padx=5)
        radius_entry.insert(0, default_radius)
       
        section_label = ttk.Label(root, text=r"Sp (cm²):")
        section_label.pack(side=tk.LEFT, padx=5)
        section_entry = ttk.Entry(root, width=10)
        section_entry.pack(side=tk.LEFT, padx=5)
        section_entry.insert(0, default_section)
        
        def update_plot():
            try:
                # box impedance
                self.Cab = self.Vb / rho / c ** 2
                self.Rab = rho * c / eta / self.Vb
                self.Zbox = gtb.parallel(1 / s / self.Cab, self.Rab)
        
                # port impedance
                self.Map = rho * self.Lp / self.Sp
                self.Mal = 0.85 * 2 * self.rp
                self.Map += self.Mal
                Pp = 2 * np.pi * self.rp
                alpha = np.sqrt(f_axis) * (0.95e-5 + 2.03e-5) * Pp / 2 / self.Sp
                kl = k * (1 + alpha * (1 - 1j))
                d0 = (0.6133 + 0.85) * self.rp
                self.Zrad = rho * c / self.Sp * (1 / 4 * (kl * self.rp) ** 2 + 1j * kl * d0)
                self.Zp = s * self.Map + self.Zrad
                self.Zab = gtb.parallel(1 / s / self.Cab, self.Rab, s * self.Map + self.Zrad)
    
                # Calculate total system
                ZaTot = driver.Zs
                Ps = driver.Ps
                Qs = Ps / (ZaTot + self.Zab)
                Qp = Qs * self.Zbox / (self.Zbox + self.Zp)
    
                p_s = 1j * k * rho * c * Qs * np.exp(-1j * k * 1) / (2 * np.pi * 1)
                p_p = -1j * k * rho * c * Qp * np.exp(-1j * k * 1) / (2 * np.pi * 1)
                Ze = driver.Ze + driver.Bl ** 2 / (driver.Zms + driver.Sd**2 * (self.Zab + driver.Zaf))
        
                # Clear the axes and plot new data
                ax_spl.clear()
                ax_imp.clear()
                ax_spl.semilogx(f_axis, gtb.gain.SPL(p_s), label='driver')
                ax_spl.semilogx(f_axis, gtb.gain.SPL(p_p), label='port')
                ax_spl.semilogx(f_axis, gtb.gain.SPL(p_s+p_p), label='total')
                ax_imp.semilogx(f_axis, np.abs(Ze), label='Impedance')
                ax_spl.set_ylabel('SPL [dB]')
                ax_spl.set_ylim(np.min(gtb.gain.SPL(p_s+p_p)), 
                                np.max(gtb.gain.SPL(p_s+p_p))+6)
                ax_imp.set_ylabel('Impedance [Ohm]')
                ax_imp.set_xlabel('Frequency [Hz]')
                ax_spl.grid(which='both')
                ax_imp.grid(which='both')
                ax_spl.legend(loc='best')
                ax_imp.legend(loc='best')
                canvas.draw()
            except ValueError:
                pass
                
                
        def update_volume():
            try:
                self.Vb = float(volume_entry.get()) * 1e-3  # Convert L to m^3
                update_plot()
            except ValueError:
                print("Invalid input for volume. Please enter a numeric value.")
        
        def update_length():
            try:
                self.Lp = float(length_entry.get()) * 1e-2  # Convert cm to m
                update_plot()
            except ValueError:
                print("Invalid input for length. Please enter a numeric value.")
        
        def update_radius():
            try:
                self.rp = float(radius_entry.get()) * 1e-2  # Convert cm to m
                self.Sp = np.pi * self.rp ** 2
                section_entry.delete(0, "end")
                section_entry.insert(0, str(round(self.Sp * 1e4, 2)))  # Update section box in cm^2
                update_plot()
            except ValueError:
                print("Invalid input for radius. Please enter a numeric value.")
    
        def update_section():
            try:
                self.Sp = float(section_entry.get()) * 1e-4  # Convert cm² to m²
                self.rp = np.sqrt(self.Sp / np.pi)
                radius_entry.delete(0, "end")
                radius_entry.insert(0, str(round(self.rp * 1e2, 2)))  # Update radius box in cm
                update_plot()
            except ValueError:
                print("Invalid input for section. Please enter a numeric value.")
        
    
        # bind Return and FocusOut
        volume_entry.bind("<Return>", lambda event: update_volume())  # Bind Enter key
        volume_entry.bind("<FocusOut>", lambda event: update_volume())
        length_entry.bind("<Return>", lambda event: update_length())  # Bind Enter key
        length_entry.bind("<FocusOut>", lambda event: update_length())
        radius_entry.bind("<Return>", lambda event: update_radius())
        radius_entry.bind("<FocusOut>", lambda event: update_radius())
        section_entry.bind("<Return>", lambda event: update_section())
        section_entry.bind("<FocusOut>", lambda event: update_section())
      
        # Initial plot
        update_plot()
    
        # Run the tkinter loop
        root.mainloop()
    
    def exportZe(self, folder_name, file_name):
        import os
        
        # Create the folder if it doesn't exist
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        
        module = np.abs(self.ZeTot)
        phase = np.angle(self.ZeTot, deg=True)
        path = os.path.join(folder_name, file_name)
        np.savetxt(path, np.array([self.f_array, module, phase]).T,
                   fmt="%.3f", 
                   header="Freq[Hz]  Imp[Ohm]  Phase[Deg]",
                   delimiter=',',
                   comments='')

    

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
        Le = 1e-12     # otherwise it doesn't work with circuitSolver()
    
    
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
    
    
    
