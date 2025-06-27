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
        from ..speakerSim.enclosureDesign import speakerBox
        
        driver = copy(self)
        driver.Le = 1e-12
        driver_b = copy(self)
        c = self.c
        rho = self.rho
    
        ## box parameters
        self.Vb = driver.Vas
        self.fc = driver.Fs * np.sqrt(driver.Vas / self.Vb + 1)
        self.Qtc = self.fc / driver.Fs * driver.Qts
        self.Qab = 120
        self.Qal = 30
        
        ## radiated pressure at 1 m
        f_axis = driver.f_array
        omega = 2 * np.pi * f_axis
        k = omega / c
        s = 1j*omega
        
        # Setup Tkinter window
        root = tk.Tk()
        root.title("Sealed Alignment")
    
        # Create figure (using matplotlib's Figure class, not plt.subplots)
        fig = Figure()
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)
            
        # Embed the plot into Tkinter
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
        # Create Labels and Entry fields
        defaultQtc = round(self.Qtc, 4)
        defaultVb  = round(self.Vb * 1e3, 4)
        defaultQab = 120
        defaultQal = 30
        
        label_qtc = ttk.Label(root, text="Qtc:")
        label_qtc.pack(side=tk.LEFT, padx=5)
        entry_qtc = ttk.Entry(root, width=10)
        entry_qtc.pack(side=tk.LEFT, padx=5)
        entry_qtc.insert(0, str(defaultQtc))
    
        label_vb = ttk.Label(root, text="Vb (L):")
        label_vb.pack(side=tk.LEFT, padx=5)
        entry_vb = ttk.Entry(root, width=10)
        entry_vb.pack(side=tk.LEFT, padx=5)
        entry_vb.insert(0, str(defaultVb))
        
        label_qab = ttk.Label(root, text="Qab:")
        label_qab.pack(side=tk.LEFT, padx=5)
        entry_qab = ttk.Entry(root, width=10)
        entry_qab.pack(side=tk.LEFT, padx=5)
        entry_qab.insert(0, str(defaultQab))
        
        label_qal = ttk.Label(root, text="Qal:")
        label_qal.pack(side=tk.LEFT, padx=5)
        entry_qal = ttk.Entry(root, width=10)
        entry_qal.pack(side=tk.LEFT, padx=5)
        entry_qal.insert(0, str(defaultQal))
        
        # Function to update the plot based on Qtc entry
        def update_plot():
            try:
                box = speakerBox(self.Vb, frequencyRange=driver.f_array,
                                 Qab=self.Qab, Qal=self.Qal)
                box.getDriverResponse(driver)
                
                p_s = 1j * k * rho * c * box.Q * np.exp(-1j * k * 1) / (2 * np.pi * 1)
                
                box.getDriverResponse(driver_b)
                Ze = box.Ze
                v = box.v
                P = np.abs(box.network.getPotential(1) * box.network.getFlow(1))
                
                # Clear and update the plots
                ax1.clear()
                ax2.clear()
                ax3.clear()
                ax4.clear()
    
                ax1.semilogx(f_axis, gtb.gain.SPL(p_s), "b")
                ax2.semilogx(f_axis, np.abs(Ze), "b")
                ax3.semilogx(f_axis, np.abs(v/s) * 1e3, "b")
                ax4.semilogx(f_axis, P, "b")
                
                ax1.set(ylabel="SPL [dB]")
                ax2.set(ylabel="Impedance [Ohm]")
                ax3.set(xlabel="Frequency [Hz]", ylabel="Excursion [mm]")
                ax4.set(xlabel="Frequency [Hz]", ylabel="Power [W]")
                
                for ax in [ax1, ax2, ax3, ax4]:
                    ax.grid(which='both', linestyle="dotted")
                    
                canvas.draw()
    
            except ValueError:
                pass  # Prevent the function from crashing if non-numeric input is entered
    
        # Function to update Qtc based on Volume entry
        def update_vb():
            try:
                self.Vb = float(entry_vb.get()) * 1e-3
                self.fc = driver.Fs * np.sqrt(driver.Vas / self.Vb + 1)
                self.Qtc = self.fc / driver.Fs * driver.Qts
                
                # update qtc value
                entry_qtc.delete(0, tk.END)
                entry_qtc.insert(0, str(round(self.Qtc, 4)))
                
                update_plot()
                
            except ValueError:
                pass    
            
        def update_qtc():
            try:
                self.Qtc = float(entry_qtc.get())
                self.fc = self.Qtc / driver.Qts * driver.Fs
                self.Vb = driver.Vas / ((self.fc / driver.Fs)**2 - 1)
            
                # Update vb value in the entry box
                entry_vb.delete(0, tk.END)
                entry_vb.insert(0, str(round(self.Vb*1e3, 4)))
    
                update_plot()  # Automatically update plot with new values
    
            except ValueError:
                pass  # Handle non-numeric input
                
        def update_QAB():
            try:
                self.Qab = float(entry_qab.get())
                update_plot()
            except ValueError:
                pass
        
        def update_QAL():
            try:
                self.Qal = float(entry_qal.get())
                update_plot()
            except ValueError:
                pass
    
        # Bind events to update the plot automatically when the user presses Enter or leaves the entry field
        entry_qtc.bind("<Return>", lambda event: update_qtc())
        entry_qtc.bind("<FocusOut>", lambda event: update_qtc())
        
        entry_vb.bind("<Return>", lambda event: update_vb())
        entry_vb.bind("<FocusOut>", lambda event: update_vb())
    
        entry_qab.bind("<Return>", lambda event: update_QAB())
        entry_qab.bind("<FocusOut>", lambda event: update_QAB())
    
        entry_qal.bind("<Return>", lambda event: update_QAL())
        entry_qal.bind("<FocusOut>", lambda event: update_QAL())

        
        # initial plot
        update_plot()
        fig.tight_layout()


        root.mainloop()
     
    
    def portedAlignment(self):
        from ..speakerSim.enclosureDesign import speakerBox
        
        driver = copy(self)
        driver_b = copy(self)
        driver.Le = 1e-12
        c = self.c
        rho = self.rho
        f_axis = driver.f_array
        omega = 2 * np.pi * f_axis
        k = omega / c
        s = 1j * omega
        eta = 1e-5
    
        # Default parameters
        self.driver = driver
        self.Vb = copy(driver.Vas)
        self.Lp = 343 / driver.Fs / 100  # Length in meters
        self.rp = self.Lp / 2  # Radius in meters
        self.Sp = np.pi * self.rp ** 2  # Port cross-sectional area
        self.Qab = 120
        self.Qal = 30
        
        # Create input widgets
        default_volume  = str(round(self.Vb * 1e3, 2))
        default_length  = str(round(self.Lp * 1e2, 2))
        default_radius  = str(round(self.rp * 1e2, 2))
        default_section = str(round(self.Sp * 1e4, 2))
        default_Qab     = str(self.Qab)
        default_Qal     = str(self.Qal)
        
        # GUI creation
        root = tk.Tk()
        root.title("Ported Alignment")
    
        # Create a matplotlib figure
        fig = Figure()
        ax_spl = fig.add_subplot(221)
        ax_imp = fig.add_subplot(222)
        ax_vx  = fig.add_subplot(223)
        ax_px  = fig.add_subplot(224)
        ax_xx  = ax_vx.twinx()

        # fig.tight_layout()
        
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
        
        qab_label = ttk.Label(root, text=r"Qab:")
        qab_label.pack(side=tk.LEFT, padx=5)
        qab_entry = ttk.Entry(root, width=10)
        qab_entry.pack(side=tk.LEFT, padx=5)
        qab_entry.insert(0, default_Qab)
        
        qal_label = ttk.Label(root, text=r"Qal:")
        qal_label.pack(side=tk.LEFT, padx=5)
        qal_entry = ttk.Entry(root, width=10)
        qal_entry.pack(side=tk.LEFT, padx=5)
        qal_entry.insert(0, default_Qal)
        
        def update_plot():
            try:
                # box impedance
                box = speakerBox(self.Vb, self.f_array, 
                                 Qab=self.Qab, Qal=self.Qal, 
                                 Lp=self.Lp, Sp=self.Sp)
                box.getDriverResponse(driver)
                p_s = 1j * k * rho * c * box.Q * np.exp(-1j * k * 1) / (2 * np.pi * 1)
                p_p = 1j * k * rho * c * box.Qp * np.exp(-1j * k * 1) / (2 * np.pi * 1)
                box.getDriverResponse(driver_b)
                Ze = box.Ze
                v  = box.v
                vp = box.vp
                P  = box.network.getPotential(1) * box.network.getFlow(1)
                
                # Clear the axes and plot new data
                ax_spl.clear()
                ax_imp.clear()
                ax_vx.clear()
                ax_xx.clear()
                ax_px.clear()
                ax_spl.semilogx(f_axis, gtb.gain.SPL(p_s), "b", label='driver')
                ax_spl.semilogx(f_axis, gtb.gain.SPL(p_p), "r", label='port')
                ax_spl.semilogx(f_axis, gtb.gain.SPL(p_s+p_p),"k", label='total')
                ax_spl.set_ylim(np.max(gtb.gain.SPL(p_s+p_p))-30, 
                                np.max(gtb.gain.SPL(p_s+p_p))+6)

                ax_imp.semilogx(f_axis, np.abs(Ze), "b", label='Impedance')
                
                # ax_vx.semilogx(f_axis, np.abs(v), label="driver")
                ax_vx.semilogx(f_axis, np.abs(vp), color="red", label="port")
                ax_xx.semilogx(f_axis, np.abs(v/s) * 1e3, color="blue", label="driver")
                
                ax_px.semilogx(f_axis, np.abs(P), "b")
                
                # labels
                ax_spl.set_ylabel('SPL [dB]')
                ax_imp.set_ylabel('Impedance [Ohm]')
                ax_vx.set(xlabel="Frequency [Hz]")
                ax_vx.set_ylabel("Port velocity [m/s]", color="red")
                ax_xx.set_ylabel("Driver excursion [mm]", color="blue")
                ax_px.set(xlabel="Frequency [Hz]", ylabel="Power [W]")
                ax_xx.yaxis.set_label_position("right")
                
                ax_vx.tick_params(axis='y', labelcolor="red")
                ax_xx.tick_params(axis='y', labelcolor="blue")
                
                # ax_xx.set_yticks(np.linspace(ax_xx.get_yticks()[0],
                #                              ax_xx.get_yticks()[-1], 
                #                              len(ax_vx.get_yticks())))
                
                # grids
                ax_spl.grid(which='both', linestyle="dotted")
                ax_imp.grid(which='both', linestyle="dotted")
                # ax_vx.grid(which='both', linestyle="dotted")
                # ax_xx.grid(which='both', linestyle="dotted")
                ax_px.grid(which='both', linestyle="dotted")
                
                # legend location
                ax_spl.legend(loc='best')
                # ax_imp.legend(loc='best')
                # ax_vx.legend(loc='upper right')
                # ax_xx.legend(loc='center right')
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
        
        def update_QAB():
            try:
                self.Qab = float(qab_entry.get())
                update_plot()
            except ValueError:
                print("Invalid input. Please enter a numeric value.")
        
        def update_QAL():
            try:
                self.Qal = float(qal_entry.get())
                update_plot()
            except ValueError:
                print("Invalid input. Please enter a numeric value.")
        
        
        # bind Return and FocusOut
        volume_entry.bind("<Return>", lambda event: update_volume())  # Bind Enter key
        volume_entry.bind("<FocusOut>", lambda event: update_volume())
        length_entry.bind("<Return>", lambda event: update_length())  # Bind Enter key
        length_entry.bind("<FocusOut>", lambda event: update_length())
        radius_entry.bind("<Return>", lambda event: update_radius())
        radius_entry.bind("<FocusOut>", lambda event: update_radius())
        section_entry.bind("<Return>", lambda event: update_section())
        section_entry.bind("<FocusOut>", lambda event: update_section())
        qab_entry.bind("<Return>", lambda event: update_QAB())
        qab_entry.bind("<FocusOut>", lambda event: update_QAB())
        qal_entry.bind("<Return>", lambda event: update_QAL())
        qal_entry.bind("<FocusOut>", lambda event: update_QAL())
        
        # Initial plot
        update_plot()
        fig.tight_layout()
    
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
    
    
    
