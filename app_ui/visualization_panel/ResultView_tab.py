from PyQt6.QtWidgets import QWidget
from app_ui.visualization_panel.plotting_tools.directivity_plotter import MplDirectivityCanvas
from app_ui.visualization_panel.plotting_tools.field_plotter import pyvistaFieldPlotter
from pyvistaqt import QtInteractor
import pyvista as pv
import numpy as np
import generalToolbox as gtb
from generalToolbox.geometry import createPlaneArray


class ResultView:
    def __init__(self, control_panel, name, evaluations, filters, geo, frequency, system):
        self.mesh_actor = None
        self.control_panel = control_panel
        self.name = name
        self.parameters = evaluations[name]  # parameters of current observation (minAngle, maxAngle, plane, etc.)
        self.plot_param = {}
        self.results = evaluations["results"][name]  # single computed observation from the self.result of MainWindow
        self.filters = filters
        self.geo = geo
        self.frequency = frequency
        self.system = system
        self.type = None

        self.tab = QWidget()
        self.plotter = None

        if evaluations[name]["type"] == "polar":
            self.createPolarRadiationPlotter()
            self.type = "polar"
        elif evaluations[name]["type"] == "screen":
            self.createPressureFieldPlotter()
            self.plotter.plotter.resize(810, 600)
            self.type = "screen"

        control_panel.addTab(self.tab, name)

    ## CREATE POLAR PLOT
    def createPolarRadiationPlotter(self):
        """
        Create a polar radiation plotter if the current observation is type=="polar"

        :return:
        """
        minAngle = self.parameters["minAngle"]
        maxAngle = self.parameters["maxAngle"]
        step = self.parameters["step"]
        theta = np.arange(minAngle, maxAngle + step, step)
        self.plot_param["theta"] = theta

        self.plotter = MplDirectivityCanvas(self.tab, 8.1, 6.1)
        self.plotter.ax_SPL_directivity.set_xscale('log')
        self.plotter.ax_dB_directivity.set_xscale('log')
        self.plotter.ax_SPL_response.set(xlabel="Frequency [Hz]", ylabel="SPL [dB]")
        self.plotter.ax_SPL_response.grid(linestyle="dotted", which="both")
        self.plotter.ax_SPL_directivity.set(xlabel="Frequency [Hz]")

        pMic = self.system.get_pMic("GUI", self.name)
        self.plotter.initialize_plot(self.system.frequency, theta, pMic)

    ## WILL DISTRIBUTE WHERE TO SEND UPDATED PLOTS
    def updatePlot(self, surface2plot, param_update):
        if self.type == "polar":
            self.updatePolarRadiationPlotter(surface2plot, param_update["polar"])
        if self.type == "screen":
            self.updatePressureFieldPlotter(surface2plot, param_update["fields"])

    ## UPDATE POLAR PLOT
    def updatePolarRadiationPlotter(self, surface2sum, param_update):
        """
        Update polar plotter
        :return:
        """

        pMic = self.system.get_pMic("GUI", self.name, radiatingSurface=surface2sum)  # pressure to plot
        self.plotter.update_plot(self.system.frequency, self.plot_param["theta"], pMic,
                                 parameters=param_update)

    def createPressureFieldPlotter(self):
        self.plotter = pyvistaFieldPlotter(self.geo, self.tab)

        freq_idx, _ = gtb.findInArray(self.frequency, 250)

        # pressure over field
        xmic_tmp, length, width = createPlaneArray(self.parameters["L1"], self.parameters["L2"],
                                                   self.parameters["step"], self.parameters["plane"],
                                                   self.parameters["offset"])
        # print("observation name::: ", self.name)
        pMic_results = self.system.get_pMic("GUI", self.name)
        # print("SHAPE PMIC::::: ", pMic_results.shape)
        pMic_results2plot = np.reshape(pMic_results[freq_idx, :], [len(length), len(width)])

        # pressure over mesh
        pMic_spk = sumPressureArray(self.system.acoustic_study["GUI"])[freq_idx, :]  #self.system.acoustic_study["GUI"].getMeshPressure(freq_idx=freq_idx)
        pMic_spk = pMic_spk[:self.system.acoustic_study["GUI"].spaceP.grid_dof_count // self.system.acoustic_study[
            "GUI"].sizeFactor]
        toAdd = self.system.acoustic_study["GUI"].vertices - self.system.acoustic_study["GUI"].spaceP.grid_dof_count
        pMic_spk2plot = np.concatenate((pMic_spk, np.zeros(toAdd)))

        self.plotter.initialize_plot(self.name, pMic_spk2plot, pMic_results2plot, length, width, xmic_tmp)

    def addField(self, name, evaluation_setup):
        from generalToolbox.geometry import createPlaneArray
        parameters = evaluation_setup[name]
        xmic_tmp, length, width = createPlaneArray(parameters["L1"], parameters["L2"],
                                                   parameters["step"], parameters["plane"], parameters["offset"])

        freq_idx, _ = gtb.findInArray(self.frequency, 250)
        pMic_results = self.system.get_pMic("GUI", name)
        pMic_results2plot = np.reshape(pMic_results[freq_idx, :], [len(length), len(width)])
        self.plotter.addField(name, pMic_results2plot, length, width, xmic_tmp)

    def updatePressureFieldPlotter(self, surface2plot, param_update):
        freq_idx, _ = gtb.findInArray(self.system.frequency, param_update["freq2plot"])

        for field in self.plotter.fields:  # update fields value to selected surfaces and frequency
            pMic_results = self.system.get_pMic("GUI", field, radiatingSurface=surface2plot)
            pMic_results2plot = np.reshape(pMic_results[freq_idx, :],
                                           [self.plotter.length[field], self.plotter.width[field]])
            self.plotter.fieldPressure[field] = pMic_results2plot
            self.plotter.getMinMaxPressure(field, pMic_results2plot)

        # pressure over system mesh (simulation loudspeaker)
        pMic_spk = sumPressureArray(self.system.acoustic_study["GUI"], surface2plot)[freq_idx, :]
        pMic_spk = pMic_spk[:self.system.acoustic_study["GUI"].spaceP.grid_dof_count // self.system.acoustic_study[
            "GUI"].sizeFactor]
        toAdd = self.system.acoustic_study["GUI"].vertices - self.system.acoustic_study["GUI"].spaceP.grid_dof_count
        pMic_spk2plot = np.concatenate((pMic_spk, np.zeros(toAdd)))

        self.plotter.update_plot(param_update, pMic_spk2plot)



def sumPressureArray(bemObj, radiatingSurface='all'):
    # print(radiatingSurface)
    p_mesh = bemObj.p_mesh
    radSurf_system = bemObj.radiatingSurface
    if radiatingSurface == 'all':
        radiatingSurface = radSurf_system

    pressureCoeff = np.zeros([len(bemObj.freq_array), bemObj.spaceP.grid_dof_count], dtype='complex')
    for f in range(len(bemObj.freq_array)):
        for j in range(len(radiatingSurface)):
            ind_surface = np.argwhere(radSurf_system == radiatingSurface[j])[0][0]
            pressureCoeff[f, :] += p_mesh[f][ind_surface].coefficients
    return pressureCoeff



