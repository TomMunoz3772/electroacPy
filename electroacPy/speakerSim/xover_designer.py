#!/usr/bin/python3
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog
from electroacPy.speakerSim.filterDesign import xovers
from electroacPy.loudspeakerSystem import updateResults, apply_Velocity_From_EAD, apply_Velocity_From_PLV, apply_Velocity_From_SPKBOX
from electroacPy.io import load, save
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
import generalToolbox as gtb


class XoverAppApp:
    def __init__(self, master=None):
        # General properties
        self.data_path = ""
        self.pMic = {}
        self.frequency = 0
        self.observation = {}
        self.angle = {}
        self.xovers = {}
        self.canvas_widget = {}
        self.fig = {}
        self.canvas = {}
        self.fig_sum = {}
        self.canvas_sum = {}
        self.dat = None
        self.study = None

        # build ui
        self.toplevel1 = tk.Tk() if master is None else tk.Toplevel(master)
        self.toplevel1.configure(height=200, width=200)
        self.toplevel1.title("ElectroacPy - crossover designer")
        self.frame1 = ttk.Frame(self.toplevel1, name="frame1")
        self.frame1.configure(height=720, width=1250)
        self.panedWindow = ttk.Panedwindow(self.frame1, orient='horizontal')
        self.panedWindow.configure(height=720, width=1250)
        self.notebook = ttk.Notebook(self.panedWindow, name="notebook")
        self.notebook.configure(height=720, width=250)

        # =======================
        # DATASET TAB
        self.frame_import = ttk.Frame(self.notebook, name="frame_import")
        self.frame_import.configure(height=720, width=250)
        self.button_load_data = ttk.Button(
            self.frame_import, name="button_load_data")
        self.button_load_data.configure(text='Import Data')
        self.button_load_data.place(
            anchor="nw", relx=0.05, rely=0.05, x=0, y=0)
        self.button_load_data.configure(command=self.load_data)

        self.button_export = ttk.Button(
            self.frame_import, name="button_export")
        self.button_export.configure(text='Export Filters')
        self.button_export.place(anchor="nw", relx=0.6, rely=0.05, x=0, y=0)
        self.button_export.configure(command=self.export_data)


        # max / min SPL control box
        self.separator1 = ttk.Separator(self.frame_import)
        self.separator1.configure(orient="horizontal")
        self.separator1.place(anchor="nw", rely=0.14, width=250, x=0, y=0)
        self.label_minSPL = ttk.Label(self.frame_import, name="label_minspl")
        self.label_minSPL.configure(text='min SPL')
        self.label_minSPL.place(anchor="nw", relx=0.05, rely=0.16, x=0, y=0)
        self.label_maxSPL = ttk.Label(self.frame_import, name="label_maxspl")
        self.label_maxSPL.configure(text='max SPL')
        self.label_maxSPL.place(anchor="nw", relx=0.05, rely=0.2, x=0, y=0)
        self.entry_minSPL = ttk.Entry(self.frame_import, name="entry_minspl")
        self.entry_minSPL.place(anchor="nw", relx=0.25, rely=0.16, width=50, x=0, y=0)
        self.entry_maxSPL = ttk.Entry(self.frame_import, name="entry_maxspl")
        self.entry_maxSPL.place(anchor="nw", relx=0.25, rely=0.20, width=50, x=0, y=0)

        # min / max frequency axis
        self.label_minFreq = ttk.Label(self.frame_import, name='label_minfreq')
        self.label_minFreq.configure(text="min freq")
        self.label_minFreq.place(anchor="nw", relx=0.5, rely=0.16, x=0, y=0)
        self.label_maxFreq = ttk.Label(self.frame_import, name='label_maxfreq')
        self.label_maxFreq.configure(text="max freq")
        self.label_maxFreq.place(anchor="nw", relx=0.5, rely=0.2, x=0, y=0)
        self.entry_minFreq = ttk.Entry(self.frame_import, name="entry_minfreq")
        self.entry_minFreq.place(anchor="nw", relx=0.71, rely=0.16, width=50, x=0, y=0)
        self.entry_maxFreq = ttk.Entry(self.frame_import, name="entry_maxfreq")
        self.entry_maxFreq.place(anchor="nw", relx=0.71, rely=0.20, width=50, x=0, y=0)


        # pack dataset frame
        self.frame_import.pack(side="top")
        self.notebook.add(self.frame_import, text='Dataset')

        # ======================================================================
        # DRIVER TAB DICT - number of tab should be set as the number of drivers
        self.filt_index = ['A', 'B', 'C', 'D', 'E', 'F', 'delay']
        self.rely = [0.05, 0.09, 0.13]
        self.rely_jump = 0.14
        self.DRIVER_FRAME = {}
        self.cb_entry = {}
        self.cb_label = {}
        self.freq_entry = {}
        self.freq_label = {}
        self.q_entry = {}
        self.q_label = {}
        self.gain_entry = {}
        self.gain_label = {}
        self.delay_entry = {}
        self.delay_label = {}
        self.delay_label_ms = {}
        self.ggain_entry = {} # general gain
        self.ggain_label = {} # general gain
        self.ggain_label_db = {} # general gain
        self.flipPhase_box = {}
        self.flipPhase_var = {}
        # self.flipPhase_label = {}
        self.update_button = {}
        self.filter_name = {}
        self.filter_label = {}

        # help with the phase flipping
        self.wasFlipped = {}

        # =================================
        # PLOTTING FRAME
        self.notebook.grid(column=0, row=0)
        self.frame_plot = ttk.Frame(self.panedWindow, name="frame_plot")
        self.frame_plot.configure(height=720, width=1000)
        self.frame_plot.grid(column=1, row=0)

        self.notebook_plot = ttk.Notebook(self.frame_plot,
                                          name="notebook_plot")  # each tab will be a notebook_obs entry
        self.notebook_plot.configure(height=720, width=1000)
        self.notebook_obs = {}      # each tab will be a specific driver
        self.PLOT_FRAME = {}        # to store plot frame in sub-notebooks (for each observation)

        # panedwindow setup
        self.panedWindow.add(self.notebook, weight=1)
        self.panedWindow.add(self.frame_plot, weight=1)
        self.panedWindow.grid(column=0, row=1)

        # plot data (for each observation)
        self.theta = {}
        self.pMic_on_axis = {}
        self.pMic_30_deg = {}

        # ==================
        # finalize app frame
        # pack main frame
        self.frame1.pack()

        # Main widget
        self.mainwindow = self.toplevel1

    def run(self):
        self.mainwindow.mainloop()

    def add_driver_tab(self, index):
        str_ind = str(index)

        # ========================================================================================
        #
        #   ADD DATASET TAB
        #
        # ========================================================================================
        self.DRIVER_FRAME[str_ind] = ttk.Frame(self.notebook)
        self.DRIVER_FRAME[str_ind] = ttk.Frame(self.notebook, name="frame_DRV_{}".format(str_ind))
        self.DRIVER_FRAME[str_ind].configure(height=200, width=200)
        self.DRIVER_FRAME[str_ind].pack(side="top")
        self.notebook.add(self.DRIVER_FRAME[str_ind], text="Drv " + str_ind)

        # =======================
        # CREATE FILTER ENTRIES
        y_offset = 0
        def_freq_values = ["100", "200", "400", "800", "1600", "3200"]
        def_q_value = "0.707"
        def_g_value = "0"
        def_delay_val = "0"
        def_filter_name = "DRV_{}".format(str_ind)
        for i in range(6):
            abc = self.filt_index[i]
            filt = str_ind + abc
            # print('addtab ', filt)
            # COMBOBOX
            self.cb_entry[filt] = ttk.Combobox(self.DRIVER_FRAME[str_ind], name="cb_entry_{}".format(filt))
            self.cb_entry[filt].place(anchor="nw", relx=0.25, rely=self.rely[0] + y_offset, x=0, y=0)
            self.cb_entry[filt]['values'] = ("Lowpass", "Highpass", "Bandpass", "Lowshelf", "Highshelf", "Peaking")
            self.cb_entry[filt].set("Peaking")
            self.cb_label[filt] = ttk.Label(self.DRIVER_FRAME[str_ind], name="cb_label_{}".format(filt))
            self.cb_label[filt].configure(text='Filter {}'.format(abc))
            self.cb_label[filt].place(anchor="nw", relx=0.05, rely=self.rely[0] + y_offset, x=0, y=0)

            # FREQ ENTRY
            self.freq_entry[filt] = ttk.Entry(self.DRIVER_FRAME[str_ind], name="freq_entry_{}".format(filt))
            self.freq_entry[filt].place(anchor="nw", relx=0.25, rely=self.rely[1] + y_offset, width=50, x=0, y=0)
            self.freq_entry[filt].insert(0, def_freq_values[i])
            self.freq_label[filt] = ttk.Label(self.DRIVER_FRAME[str_ind], name="freq_label_{}".format(filt))
            self.freq_label[filt].configure(text='Freq')
            self.freq_label[filt].place(anchor="nw", relx=0.05, rely=self.rely[1] + y_offset, x=0, y=0)

            # Q ENTRY
            self.q_entry[filt] = ttk.Entry(self.DRIVER_FRAME[str_ind], name="q_entry_{}".format(filt))
            self.q_entry[filt].place(anchor="nw", relx=0.25, rely=self.rely[2] + y_offset, width=50, x=0, y=0)
            self.q_entry[filt].insert(0, def_q_value)
            self.q_label[filt] = ttk.Label(self.DRIVER_FRAME[str_ind], name="q_label_{}".format(filt))
            self.q_label[filt].configure(text='Q')
            self.q_label[filt].place(anchor="nw", relx=0.05, rely=self.rely[2] + y_offset, x=0, y=0)

            # GAIN ENTRY
            self.gain_entry[filt] = ttk.Entry(self.DRIVER_FRAME[str_ind], name="g_entry_{}".format(filt))
            self.gain_entry[filt].place(anchor="nw", relx=0.70, rely=self.rely[2] + y_offset, width=50, x=0, y=0)
            self.gain_entry[filt].insert(0, def_g_value)
            self.gain_label[filt] = ttk.Label(self.DRIVER_FRAME[str_ind], name='g_label_{}'.format(filt))
            self.gain_label[filt].configure(text='dB')
            self.gain_label[filt].place(anchor="nw", relx=0.6, rely=self.rely[2] + y_offset, x=0, y=0)

            y_offset += self.rely_jump
            # exit loop

        # delay
        self.delay_label[str_ind] = ttk.Label(self.DRIVER_FRAME[str_ind], name="delay_label_{}".format(str_ind))
        self.delay_label[str_ind].configure(text='Delay')
        self.delay_label[str_ind].place(anchor="nw", relx=0.05, rely=0.88, width=50, x=0, y=0)
        self.delay_entry[str_ind] = ttk.Entry(self.DRIVER_FRAME[str_ind], name="delay_entry_{}".format(str_ind))
        self.delay_entry[str_ind].place(anchor='nw', relx=0.25, rely=0.88, width=50, x=0, y=0)
        self.delay_entry[str_ind].insert(0, def_delay_val)
        self.delay_label_ms[str_ind] = ttk.Label(self.DRIVER_FRAME[str_ind], name="delay_ms_{}".format(str_ind))
        self.delay_label_ms[str_ind].configure(text='(ms)')
        self.delay_label_ms[str_ind].place(anchor="nw", relx=0.46, rely=0.88, x=0, y=0)

        # gain
        self.ggain_label[str_ind] = ttk.Label(self.DRIVER_FRAME[str_ind], name="ggain_label_{}".format(str_ind))
        self.ggain_label[str_ind].configure(text='Gain')
        self.ggain_label[str_ind].place(anchor="nw", relx=0.05, rely=0.92, width=50, x=0, y=0)
        self.ggain_entry[str_ind] = ttk.Entry(self.DRIVER_FRAME[str_ind], name="ggain_entry_{}".format(str_ind))
        self.ggain_entry[str_ind].place(anchor='nw', relx=0.25, rely=0.92, width=50, x=0, y=0)
        self.ggain_entry[str_ind].insert(0, def_delay_val)
        self.ggain_label_db[str_ind] = ttk.Label(self.DRIVER_FRAME[str_ind], name="ggain_dB_{}".format(str_ind))
        self.ggain_label_db[str_ind].configure(text='(dB)')
        self.ggain_label_db[str_ind].place(anchor="nw", relx=0.46, rely=0.92, x=0, y=0)

        # flip phase
        # self.flipPhase_var[str_ind] = tk.BooleanVar()
        self.flipPhase_box[str_ind] = ttk.Checkbutton(self.DRIVER_FRAME[str_ind],
                                                      name="flip_checkbox_{}".format(str_ind))
        self.flipPhase_box[str_ind].configure(text='Flip phase')
        self.flipPhase_box[str_ind].place(anchor="nw", relx=0.05, rely=0.96, x=0, y=0)
        self.flipPhase_box[str_ind].state(['!alternate'])

        # update button
        self.update_button[str_ind] = ttk.Button(self.DRIVER_FRAME[str_ind], name="update_filt_{}".format(str_ind),
                                                 command=self.update_xovers)
        self.update_button[str_ind].configure(text='Update Filter')
        self.update_button[str_ind].place(anchor="nw", relx=0.65, rely=0.95, x=0, y=0)

        # filter / driver name
        self.filter_name[str_ind] = ttk.Entry(self.DRIVER_FRAME[str_ind], name='filter_name_{}'.format(str_ind))
        self.filter_name[str_ind].place(anchor="nw", relx=0.25, rely=0.01, x=0, y=0)
        self.filter_name[str_ind].insert(0, def_filter_name)
        self.filter_label[str_ind] = ttk.Label(self.DRIVER_FRAME[str_ind], name='filter_label_{}'.format(str_ind))
        self.filter_label[str_ind].configure(text='Name')
        self.filter_label[str_ind].place(anchor="nw", relx=0.05, rely=0.01, x=0, y=0)
        return None

    def add_observation_tab(self, obsName, indexDriver):
        str_ind = obsName + "_" + str(indexDriver)
        text_str = "DRV_" + str(indexDriver)
        # ========================================================================================
        #
        #   ADD PLOTTING TAB
        #
        # ========================================================================================
        self.PLOT_FRAME[str_ind] = ttk.Frame(self.notebook_obs[obsName], name="frame_plot_{}".format(str_ind))
        self.PLOT_FRAME[str_ind].configure(height=200, width=200)
        self.PLOT_FRAME[str_ind].pack(side="top")
        self.notebook_obs[obsName].add(self.PLOT_FRAME[str_ind], text=text_str)  # self.filter_name[str_ind].get()

        # plot data
        # dataset = self.pMic[obsName]
        # print("add_pressure ", str_ind)
        self.plot_data(self.PLOT_FRAME[str_ind], obsName, indexDriver - 1, 0)
        return None

    def add_summed_pressure_tab(self, obsName):
        str_ind = obsName + "_all"
        text_str = "DRV_" + "_all"

        self.PLOT_FRAME[str_ind] = ttk.Frame(self.notebook_obs[obsName], name="frame_plot_{}".format(str_ind))
        self.PLOT_FRAME[str_ind].configure(height=200, width=200)
        self.PLOT_FRAME[str_ind].pack(side="top")
        self.notebook_obs[obsName].add(self.PLOT_FRAME[str_ind], text=text_str)

        self.plot_data_summed(self.PLOT_FRAME[str_ind], obsName, 0)
        return None

    def load_data(self):
        self.data_path = filedialog.askdirectory()  # file should be an observation npz archive
        dat = load(self.data_path)
        self.dat = dat

        for study in dat.acoustic_study:
            self.frequency = dat.acoustic_study[study].freq_array
            for obs in range(len(dat.observation[study].observationName)):
                if dat.observation[study].observationType[obs] == 'polar':
                    obsName = dat.observation[study].observationName[obs]
                    updateResults(dat, study, bypass_xover=True)
                    self.observation[obsName] = obsName
                    self.pMic[obsName] = dat.results[study].pMicArray[obs]
                    self.angle[obsName] = dat.observation[study].theta[obs]

        self.study = study
        self.entry_minFreq.insert(0, str(int(self.frequency[0])))
        self.entry_maxFreq.insert(0, str(int(self.frequency[-1])))
        iteration = 0
        for obs in self.observation:
            # initialize observation + plotting frames
            self.notebook_obs[obs] = ttk.Notebook(self.notebook_plot, name="notebook_obs_{}".format(obs))
            self.notebook_obs[obs].configure(width=200, height=200)

            shape_pMic = self.pMic[obs].shape
            nDriver = shape_pMic[-1]
            for i in range(nDriver):
                if iteration == 0:
                    self.xovers["drv_" + str(i + 1)] = xovers(self.frequency)
                    self.wasFlipped["drv_" + str(i+1)] = False
                    self.add_driver_tab(i + 1)
                self.add_observation_tab(obs, i + 1)
            iteration += 1
            self.add_summed_pressure_tab(obs)
            self.notebook_obs[obs].pack(side="top")
            self.notebook_plot.add(self.notebook_obs[obs], text="{}".format(obs))
        self.notebook_plot.pack(side="top")
        return None

    def export_data(self):
        for i, xo in enumerate(self.xovers):
            self.dat.crossover[xo] = self.xovers[xo]
            self.dat.crossover[xo].referenceStudy = self.study
            # print("here {}".format(xo))
            # for i in range(len(self.dat.acoustic_study[self.study].radiatingSurface)):
            self.dat.crossover[xo].ref2bem = int(self.dat.acoustic_study[self.study].radiatingSurface[i])
                # print("filter {} exported as ref to {}".format(xo, int(self.dat.acoustic_study[self.study].radiatingSurface[i])))
        save(self.data_path+"_xoverDesigner", self.dat)
        return None


    def plot_data(self, masterWindow, obs, index_driver, filters):
        # import matplotlib.pyplot as plt
        ## process dataset and store initial data
        theta = self.angle[obs] * 180 / np.pi
        idx_z, _ = gtb.findInArray(self.angle[obs], 0)
        idx_30, _ = gtb.findInArray(theta, 30)
        N_angle = len(self.angle[obs])
        h = self.xovers["drv_" + str(index_driver + 1)].h

        pMic_XO = self.pMic[obs][:, :, index_driver]

        # on-axis SPL
        SPL_on_axis = gtb.gain.SPL(pMic_XO[:, idx_z])
        SPL_30_deg = gtb.gain.SPL(pMic_XO[:, idx_30])

        # SPL directivity + get default min-max SPL values
        SPL_directivity = gtb.gain.SPL(pMic_XO)
        MAX_SPL = int(np.max(SPL_directivity) + 3)
        MIN_SPL = MAX_SPL - 40
        try:
            if int(self.entry_minSPL.get()) > MIN_SPL:
                self.entry_minSPL.delete(0, "end")
                self.entry_minSPL.insert(0, str(MIN_SPL))
            if int(self.entry_maxSPL.get()) < MAX_SPL:
                self.entry_maxSPL.delete(0, "end")
                self.entry_maxSPL.insert(0, str(MAX_SPL))
        except:
            self.entry_minSPL.insert(0, str(MIN_SPL))
            self.entry_maxSPL.insert(0, str(MAX_SPL))

        # absolute directivity
        maxAngle = np.array([np.max(abs(pMic_XO), 1)])
        maxMatrix = np.repeat(maxAngle, len(self.angle[obs]), 0).T
        ABS_directivity = gtb.gain.dB(np.abs(pMic_XO) / maxMatrix)

        # min / max frequency bounds
        minf = float(self.entry_minFreq.get())
        maxf = float(self.entry_maxFreq.get())

        # plot fig
        indexFig = obs+str(index_driver + 1)
        self.fig[indexFig] = Figure(figsize=(10.417, 7.5))
        # plot 1
        plot1 = self.fig[indexFig].add_subplot(221)
        gca = plot1.contourf(self.frequency, theta, SPL_directivity.T, np.arange(MIN_SPL, MAX_SPL, 3),
                             cmap="turbo")
        plot1.set_xscale('log')
        cbar = self.fig[indexFig].colorbar(gca)
        cbar.set_label('SPL [dB]')
        plot1.set_xlabel("Frequency [Hz]")
        plot1.set_ylabel("Angle [degree]")
        plot1.set_xlim([minf, maxf])

        # plot 2
        plot2 = self.fig[indexFig].add_subplot(222)
        plot2.semilogx(self.frequency, SPL_on_axis, 'b', label='on-axis')
        plot2.semilogx(self.frequency, SPL_30_deg, 'g', label='30 deg', alpha=0.75, linewidth=0.7)
        plot2.set_xlabel('Frequency [Hz]')
        plot2.set_ylabel('SPL [dB]')
        plot2.grid(True, linestyle='dotted', which='both')
        plot2.legend(loc='best')
        plot2.set_xlim([minf, maxf])
        plot2.set_ylim([MIN_SPL, MAX_SPL])

        # plot3
        plot3 = self.fig[indexFig].add_subplot(223)
        gca3 = plot3.contourf(self.frequency, theta, ABS_directivity.T,
                              np.arange(-21, 3, 3), cmap="turbo")
        plot3.set_xscale('log')
        cbar3 = self.fig[indexFig].colorbar(gca3)
        cbar3.set_label('Gain [dB]')
        plot3.set_xlabel("Frequency [Hz]")
        plot3.set_ylabel("Angle [degree]")
        plot3.set_xlim([minf, maxf])

        # plot 4
        plot4 = self.fig[indexFig].add_subplot(224)
        plot4.semilogx(self.frequency, gtb.gain.dB(h), 'k')
        plot4.set_xlabel('Frequency [Hz]')
        plot4.set_ylabel('Gain [dB]')
        # plot4.set_ylim([-32, 3])
        plot4.grid(True, linestyle='dotted', which='both')
        plot4.set_ylim([-32, +32])

        plot4_p = plot4.twinx()
        plot4_p.semilogx(self.frequency, np.angle(h), ':r')
        plot4_p.set_ylim([-3.2, 3.2])
        plot4_p.set_ylabel("Phase [rad]", color='red')
        plot4_p.tick_params(axis='y', labelcolor='red')

        self.fig[indexFig].tight_layout()
        self.canvas[indexFig] = FigureCanvasTkAgg(self.fig[indexFig], master=masterWindow)
        self.canvas[indexFig].draw()

        # placing the canvas on the Tkinter window
        self.canvas[indexFig].get_tk_widget().pack()
        # print("index Fig - plot data: ", indexFig)
        # self.canvas_widget[indexFig] = canvas.get_tk_widget()
        # self.canvas_widget[indexFig].pack()
        # toolbar = NavigationToolbar2Tk(self.canvas[indexFig], masterWindow)
        # toolbar.update()

        # placing the toolbar on the Tkinter window
        # self.canvas_widget[indexFig].get_tk_widget().pack()
        return None

    def update_plot_data(self, masterWindow, obs, index_driver, filters):
        ## process dataset
        theta = self.angle[obs] * 180 / np.pi
        idx_z, _ = gtb.findInArray(self.angle[obs], 0)
        idx_30, _ = gtb.findInArray(theta, 30)
        N_angle = len(self.angle[obs])
        h = self.xovers["drv_" + str(index_driver + 1)].h

        XO = np.zeros([len(self.frequency), N_angle], dtype=complex)
        for i in range(N_angle):
            XO[:, i] = h
        pMic_XO = self.pMic[obs][:, :, index_driver] * XO

        # on-axis SPL
        SPL_on_axis = gtb.gain.SPL(pMic_XO[:, idx_z])
        SPL_30_deg = gtb.gain.SPL(pMic_XO[:, idx_30])

        # SPL directivity + set min/max SPL
        SPL_directivity = gtb.gain.SPL(pMic_XO)
        MIN_SPL = int(self.entry_minSPL.get())
        MAX_SPL = int(self.entry_maxSPL.get())

        # absolute directivity
        maxAngle = np.array([np.max(abs(pMic_XO), 1)])
        maxMatrix = np.repeat(maxAngle, len(self.angle[obs]), 0).T
        ABS_directivity = gtb.gain.dB(np.abs(pMic_XO) / maxMatrix)

        # min / max frequency bounds
        minf = float(self.entry_minFreq.get())
        maxf = float(self.entry_maxFreq.get())

        # plot fig
        indexFig = obs + str(index_driver + 1)
        self.fig[indexFig].clear()
        # plot 1
        plot1 = self.fig[indexFig].add_subplot(221)
        gca = plot1.contourf(self.frequency, theta, SPL_directivity.T, np.arange(MIN_SPL, MAX_SPL, 3),
                                            cmap="turbo")
        plot1.set_xscale('log')
        cbar = self.fig[indexFig].colorbar(gca)
        cbar.set_label('SPL [dB]')
        plot1.set_xlabel("Frequency [Hz]")
        plot1.set_ylabel("Angle [degree]")
        plot1.set_xlim([minf, maxf])

        # plot 2
        plot2 = self.fig[indexFig].add_subplot(222)
        plot2.semilogx(self.frequency, SPL_on_axis, 'b', label='on-axis')
        plot2.semilogx(self.frequency, SPL_30_deg, 'g', label='30 deg', alpha=0.75, linewidth=0.7)
        plot2.set_xlabel('Frequency [Hz]')
        plot2.set_ylabel('SPL [dB]')
        plot2.grid(True, linestyle='dotted', which='both')
        plot2.set_xlim([minf, maxf])
        plot2.set_ylim([MIN_SPL, MAX_SPL])
        plot2.legend(loc='best')

        # plot3
        plot3 = self.fig[indexFig].add_subplot(223)
        gca3 = plot3.contourf(self.frequency, theta, ABS_directivity.T,
                                             np.arange(-21, 3, 3), cmap="turbo")
        plot3.set_xscale('log')
        cbar3 = self.fig[indexFig].colorbar(gca3)
        cbar3.set_label('Gain [dB]')
        plot3.set_xlabel("Frequency [Hz]")
        plot3.set_ylabel("Angle [degree]")
        plot3.set_xlim([minf, maxf])

        # plot 4
        plot4 = self.fig[indexFig].add_subplot(224)
        plot4.semilogx(self.frequency, gtb.gain.dB(h), 'k')
        plot4.set_xlabel('Frequency [Hz]')
        plot4.set_ylabel('Gain [dB]')
        plot4.set_ylim([-32, +32])
        plot4.grid(True, linestyle='dotted', which='both')

        plot4_p = plot4.twinx()
        plot4_p.semilogx(self.frequency, np.angle(h), ':r')
        plot4_p.set_ylim([-3.2, 3.2])
        plot4_p.set_ylabel("Phase [rad]", color='red')
        plot4_p.tick_params(axis='y', labelcolor='red')

        self.canvas[indexFig].draw()
        return None

    def plot_data_summed(self, masterWindow, obs, filter):
        ## process dataset
        theta = self.angle[obs] * 180 / np.pi
        idx_z, _ = gtb.findInArray(self.angle[obs], 0)
        idx_30, _ = gtb.findInArray(theta, 30)
        N_angle = len(self.angle[obs])
        N_driver = np.shape(self.pMic[obs])[-1]
        pMic_XO = np.zeros([len(self.frequency), N_angle], dtype=complex)
        for i in range(N_driver):
            h = self.xovers["drv_" + str(i + 1)].h
            XO = np.tile(h, (N_angle, 1))
            pMic_XO += self.pMic[obs][:, :, i] * XO.T

        # on-axis SPL
        SPL_on_axis = gtb.gain.SPL(pMic_XO[:, idx_z])
        SPL_30_deg = gtb.gain.SPL(pMic_XO[:, idx_30])

        # SPL directivity + set min/max SPL
        SPL_directivity = gtb.gain.SPL(pMic_XO)
        MAX_SPL = int(np.max(SPL_directivity) + 3)
        MIN_SPL = MAX_SPL - 40

        try:
            if int(self.entry_minSPL.get()) > MIN_SPL:
                self.entry_minSPL.delete(0, "end")
                self.entry_minSPL.insert(0, str(MIN_SPL))
            if int(self.entry_maxSPL.get()) < MAX_SPL:
                self.entry_maxSPL.delete(0, "end")
                self.entry_maxSPL.insert(0, str(MAX_SPL))
        except:
            self.entry_minSPL.insert(0, str(MIN_SPL))
            self.entry_maxSPL.insert(0, str(MAX_SPL))

        # absolute directivity
        maxAngle = np.array([np.max(abs(pMic_XO), 1)])
        maxMatrix = np.repeat(maxAngle, len(self.angle[obs]), 0).T
        ABS_directivity = gtb.gain.dB(np.abs(pMic_XO) / maxMatrix)

        # min / max frequency bounds
        minf = float(self.entry_minFreq.get())
        maxf = float(self.entry_maxFreq.get())

        #
        index_fig = obs + "_all"
        self.fig_sum[index_fig] = Figure(figsize=(10.417, 7.5))
        # plot 1
        plot1 = self.fig_sum[index_fig].add_subplot(221)
        gca = plot1.contourf(self.frequency, theta, SPL_directivity.T, np.arange(MIN_SPL, MAX_SPL, 3),
                             cmap="turbo")
        plot1.set_xscale('log')
        cbar = self.fig_sum[index_fig].colorbar(gca)
        cbar.set_label('SPL [dB]')
        plot1.set_xlabel("Frequency [Hz]")
        plot1.set_ylabel("Angle [degree]")
        plot1.set_xlim([minf, maxf])

        # plot 2
        plot2 = self.fig_sum[index_fig].add_subplot(222)
        plot2.semilogx(self.frequency, SPL_on_axis, 'b', label="on-axis")
        plot2.semilogx(self.frequency, SPL_30_deg, 'g', label="30 deg", alpha=0.75, linewidth=0.7)
        plot2.set_xlabel('Frequency [Hz]')
        plot2.set_ylabel('SPL [dB]')
        plot2.set_xlim([minf, maxf])
        plot2.set_ylim([MIN_SPL, MAX_SPL])
        plot2.grid(True, linestyle='dotted', which='both')
        plot2.legend(loc='best')

        # plot3
        plot3 = self.fig_sum[index_fig].add_subplot(223)
        gca3 = plot3.contourf(self.frequency, theta, ABS_directivity.T,
                              np.arange(-21, 3, 3), cmap="turbo")
        plot3.set_xscale('log')
        cbar3 = self.fig_sum[index_fig].colorbar(gca3)
        cbar3.set_label('Gain [dB]')
        plot3.set_xlabel("Frequency [Hz]")
        plot3.set_ylabel("Angle [degree]")
        plot3.set_xlim([minf, maxf])

        # plot 4
        plot4 = self.fig_sum[index_fig].add_subplot(224)
        plot4.semilogx(self.frequency, gtb.gain.dB(h), 'k')
        plot4.set_xlabel('Frequency [Hz]')
        plot4.set_ylabel('Gain [dB]')
        plot4.set_ylim([-32, 32])
        plot4.grid(True, linestyle='dotted', which='both')

        plot4_p = plot4.twinx()
        plot4_p.semilogx(self.frequency, np.angle(h), ':r')
        plot4_p.set_ylim([-3.2, 3.2])
        plot4_p.set_ylabel("Phase [rad]", color='red')
        plot4_p.tick_params(axis='y', labelcolor='red')

        self.fig_sum[index_fig].tight_layout()
        self.canvas_sum[index_fig] = FigureCanvasTkAgg(self.fig_sum[index_fig], master=masterWindow)
        self.canvas_sum[index_fig].draw()

        # placing the canvas on the Tkinter window
        self.canvas_sum[index_fig].get_tk_widget().pack()
        return None


    def update_plot_data_summed(self, masterWindow, obs, filter):
        ## process dataset
        theta = self.angle[obs] * 180 / np.pi
        idx_z, _ = gtb.findInArray(self.angle[obs], 0)
        idx_30, _ = gtb.findInArray(theta, 30)
        N_angle = len(self.angle[obs])
        N_driver = np.shape(self.pMic[obs])[-1]
        pMic_XO = np.zeros([len(self.frequency), N_angle], dtype=complex)
        h_to_plot = np.zeros(len(self.frequency), dtype=complex)
        for i in range(N_driver):
            h = self.xovers["drv_" + str(i + 1)].h
            h_to_plot += h
            XO = np.tile(h, (N_angle, 1))
            pMic_XO += self.pMic[obs][:, :, i] * XO.T

        # on-axis SPL
        SPL_on_axis = gtb.gain.SPL(pMic_XO[:, idx_z])
        SPL_30_deg = gtb.gain.SPL(pMic_XO[:, idx_30])

        # SPL directivity
        SPL_directivity = gtb.gain.SPL(pMic_XO)
        MIN_SPL = int(self.entry_minSPL.get())
        MAX_SPL = int(self.entry_maxSPL.get())
        # MAX_SPL = int(np.max(SPL_directivity) + 3)
        # MIN_SPL = MAX_SPL - 40


        # absolute directivity
        maxAngle = np.array([np.max(abs(pMic_XO), 1)])
        maxMatrix = np.repeat(maxAngle, len(self.angle[obs]), 0).T
        ABS_directivity = gtb.gain.dB(np.abs(pMic_XO) / maxMatrix)

        # min / max frequency bounds
        minf = float(self.entry_minFreq.get())
        maxf = float(self.entry_maxFreq.get())

        index_fig = obs + "_all"
        self.fig_sum[index_fig].clear()
        # plot 1
        plot1 = self.fig_sum[index_fig].add_subplot(221)
        gca = plot1.contourf(self.frequency, theta, SPL_directivity.T, np.arange(MIN_SPL, MAX_SPL, 3),
                             cmap="turbo")
        plot1.set_xscale('log')
        cbar = self.fig_sum[index_fig].colorbar(gca)
        cbar.set_label('SPL [dB]')
        plot1.set_xlabel("Frequency [Hz]")
        plot1.set_ylabel("Angle [degree]")
        plot1.set_xlim([minf, maxf])

        # plot 2
        plot2 = self.fig_sum[index_fig].add_subplot(222)
        plot2.semilogx(self.frequency, SPL_on_axis, 'b', label='on-axis')
        plot2.semilogx(self.frequency, SPL_30_deg, 'g', label='30 deg', alpha=0.75, linewidth=0.7)
        plot2.set_xlabel('Frequency [Hz]')
        plot2.set_ylabel('SPL [dB]')
        plot2.grid(True, linestyle='dotted', which='both')
        plot2.set_xlim([minf, maxf])
        plot2.set_ylim([MIN_SPL, MAX_SPL])
        plot2.legend(loc='best')

        # plot3
        plot3 = self.fig_sum[index_fig].add_subplot(223)
        gca3 = plot3.contourf(self.frequency, theta, ABS_directivity.T,
                              np.arange(-21, 3, 3), cmap="turbo")
        plot3.set_xscale('log')
        cbar3 = self.fig_sum[index_fig].colorbar(gca3)
        cbar3.set_label('Gain [dB]')
        plot3.set_xlabel("Frequency [Hz]")
        plot3.set_ylabel("Angle [degree]")
        plot3.set_xlim([minf, maxf])

        # plot 4
        plot4 = self.fig_sum[index_fig].add_subplot(224)
        plot4.semilogx(self.frequency, gtb.gain.dB(h_to_plot), 'k')
        plot4.set_xlabel('Frequency [Hz]')
        plot4.set_ylabel('Gain [dB]')
        plot4.set_ylim([-32, 32])
        plot4.grid(True, linestyle='dotted', which='both')

        plot4_p = plot4.twinx()
        plot4_p.semilogx(self.frequency, np.angle(h_to_plot), ':r')
        plot4_p.set_ylim([-3.2, 3.2])
        plot4_p.set_ylabel("Phase [rad]", color='red')
        plot4_p.tick_params(axis='y', labelcolor='red')

        self.canvas_sum[index_fig].draw()
        return None


    def update_xovers(self):
        iteration = 0
        for obs in self.observation:
            N_driver = np.shape(self.pMic[obs])[-1]
            for n in range(N_driver):
                for i in range(6):
                    filt_index = str(n+1) + self.filt_index[i]
                    # print('xo: ', filt_index)
                    filtertype = self.cb_entry[filt_index].get()
                    fc = float(self.freq_entry[filt_index].get())
                    Q = float(self.q_entry[filt_index].get())
                    dB = float(self.gain_entry[filt_index].get())
                    self.assign_xover(filtertype, fc, Q, dB, str(n + 1), filt_index)
                delay = float(self.delay_entry[str(n+1)].get())
                flipPhase = self.flipPhase_box[str(n+1)].state()
                gdB = float(self.ggain_entry[str(n+1)].get()) # general gain
                self.assign_gain_delay_flip(gdB, delay, flipPhase, str(n+1))
                str_ind = obs + "_" + str(n+1)
                self.update_plot_data(self.PLOT_FRAME[str_ind], obs, n, 0)
            self.update_plot_data_summed(self.PLOT_FRAME[obs+"_all"], obs, 0)
        return None

    def assign_xover(self, filtertype, fc, Q, dB, spk_index, filt_index):
        str_drv = "drv_" + spk_index
        if filtertype == 'Lowpass':
            # print('lowpass')
            self.xovers[str_drv].addLowPassBQ(filt_index, fc, Q, dB)
        elif filtertype == 'Highpass':
            # print('highpass')
            self.xovers[str_drv].addHighPassBQ(filt_index, fc, Q, dB)
        elif filtertype == 'Bandpass':
            # print('bandpass')
            self.xovers[str_drv].addBandPassBQ(filt_index, fc, Q, dB)
        elif filtertype == 'Peaking':
            # print('peaking')
            self.xovers[str_drv].addPeakEQ(filt_index, fc, Q, dB)
        elif filtertype == 'Lowshelf':
            # print('lowshelf')
            self.xovers[str_drv].addLowShelf(filt_index, fc, Q, dB)
        elif filtertype == 'Highshelf':
            # print('highshelf')
            self.xovers[str_drv].addHighShelf(filt_index, fc, Q, dB)
        else:
            print("filter not understood")
        return None

    def assign_gain_delay_flip(self, gain, delay, isFlipped, spk_index): # could be better
        str_drv = "drv_" + spk_index
        self.xovers[str_drv].addGain('dB', gain)
        self.xovers[str_drv].addDelay('dt', delay*1e-3)
        try:
            isFlipped = isFlipped[0]
        except:
            pass
        if isFlipped == 'selected':
            if self.wasFlipped[str_drv] is False:
                self.xovers[str_drv].addPhaseFlip('pi')
                self.wasFlipped[str_drv] = True
            elif self.wasFlipped[str_drv] is True:
                pass
        elif len(isFlipped) == 0:
            if self.wasFlipped[str_drv] is True:
                self.xovers[str_drv].deleteFilter('pi')
                self.wasFlipped[str_drv] = False
            elif self.wasFlipped[str_drv] is False:
                pass
        return None

if __name__ == "__main__":
    app = XoverAppApp()
    app.run()
