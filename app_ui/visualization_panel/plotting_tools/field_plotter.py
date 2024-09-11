from PyQt6.QtWidgets import QWidget
from app_ui.visualization_panel.plotting_tools.directivity_plotter import MplDirectivityCanvas
from pyvistaqt import QtInteractor
import pyvista as pv
import numpy as np
import generalToolbox as gtb


class pyvistaFieldPlotter:
    def __init__(self, geometry, parent=None, width=810, height=600):
        self.geometry = geometry
        self.plotter = QtInteractor(parent)
        self.plotter.setMinimumSize(width, height)
        self.plotter.interactor.setMinimumSize(width, height)
        # self.plotter.interactor.setSizePolicy(1, 1)

        # prepare mesh
        # print("MESH PATH:: ", geometry["meshPath"])
        self.mesh = None

        # plot mesh
        self.plotter.add_axes()

        # store the fields to plot
        self.fields = {}
        self.fieldPressure = {}
        self.maxPressure = {}
        self.clim_max = None
        self.clim_min = None
        self.length = {}
        self.width = {}

    def initialize_plot(self, fieldName, meshPressure, fieldPressure, length, width, xMic):
        # reload mesh
        self.mesh = pv.read(self.geometry["meshPath"])

        mesh2plot = pv.StructuredGrid()
        mesh2plot.points = xMic
        mesh2plot.dimensions = [len(length), len(width), 1]
        self.fields[fieldName] = mesh2plot
        self.fieldPressure[fieldName] = fieldPressure
        self.length[fieldName] = len(length)
        self.width[fieldName] = len(width)

        self.getMinMaxPressure(fieldName, fieldPressure)
        self.clim_max = self.maxPressure[fieldName]["spl"] + 6
        self.clim_min = self.maxPressure[fieldName]["spl"] - 75

        # SPL FIELD
        self.plotter.add_mesh(self.fields[fieldName], scalars=gtb.gain.SPL(fieldPressure), show_edges=False,
                              cmap="turbo", n_colors=21,
                              show_scalar_bar=True, clim=[self.clim_min, self.clim_max], name=fieldName)

        # SPEAKER FIELD: REMINDER TO NORMALISE IT COMPARED TO MAX SPL FIELD!!!
        maxField = 10 ** (self.maxPressure[fieldName][
                              "spl"] / 20) * 2e-5  # pressure max radiated on field -> used to normalise pressure on mesh (otherwise it is too high)
        speakerPressure = gtb.normMinMax(np.abs(meshPressure), 0, maxField)
        speakerPressure = gtb.gain.SPL(speakerPressure)
        self.plotter.add_mesh(self.mesh, scalars=speakerPressure, show_edges=False, cmap="turbo", n_colors=21,
                              show_scalar_bar=False,
                              clim=[self.clim_min, self.clim_max], name="speaker")

    def addField(self, fieldName, fieldPressure, length, width, xMic):
        mesh2plot = pv.StructuredGrid()
        mesh2plot.points = xMic
        mesh2plot.dimensions = [len(length), len(width), 1]
        self.fields[fieldName] = mesh2plot
        self.fieldPressure[fieldName] = fieldPressure
        self.length[fieldName] = len(length)
        self.width[fieldName] = len(width)

        self.plotter.add_mesh(self.fields[fieldName], scalars=gtb.gain.SPL(fieldPressure), show_edges=False,
                              cmap="turbo", n_colors=21,
                              show_scalar_bar=True, clim=[self.clim_min, self.clim_max], name=fieldName)

    def update_plot(self, parameters, pressure_system):
        """
        parameters {transformation, cmap, show_edges, ...}
        :param fieldName:
        :param frequency:
        :param meshPressure:
        :param fieldPressure:
        :param length:
        :param width:
        :param xMic:
        :param parameters:
        :return:
        """
        self.mesh = pv.read(self.geometry["meshPath"])

        dBmin = parameters["dBmin"]
        dBmax = parameters["dBmax"]
        # dBstep = parameters["dBstep"]
        ncolor = parameters["n_colors"]
        cmap = parameters["cmap"]

        if parameters["transformation"] == "SPL":
            transformation = gtb.gain.SPL
            MAXABS = self.updateMinMaxScalars("abs")  # for the loudspeaker mesh normalization
            # print("MAXABS: ", MAXABS)
            # print("MAXSPL: ", gtb.gain.SPL(MAXABS))
            cmax, cmin = dBmax, dBmin
            normalized_system_pressure = gtb.normMinMax(np.abs(pressure_system), 0, MAXABS)
            normalized_system_pressure = gtb.gain.SPL(normalized_system_pressure)

        elif parameters["transformation"] == "Real":
            transformation = np.real
            MAXREAL = self.updateMinMaxScalars("real")  # for the loudspeaker mesh normalization
            cmax, cmin = dBmax, dBmin
            normalized_system_pressure = gtb.normMinMax(np.real(pressure_system), -MAXREAL, MAXREAL)
            normalized_system_pressure = np.real(normalized_system_pressure)

        elif parameters["transformation"] == "Phase":
            transformation = np.angle
            MAXREAL = self.updateMinMaxScalars("real")  # for the loudspeaker mesh normalization
            cmax, cmin = np.pi, -np.pi
            # normalized_system_pressure = gtb.normMinMax(np.abs(pressure_system), -np.pi, np.pi)
            normalized_system_pressure = np.angle(pressure_system)



        scalars_bar_view = [False] * len(self.fields)
        scalars_bar_view[-1] = True

        i = 0
        for fieldName in self.fields:
            self.plotter.remove_actor(fieldName)
            self.plotter.add_mesh(self.fields[fieldName], scalars=transformation(self.fieldPressure[fieldName]),
                                  cmap=cmap, n_colors=ncolor, show_scalar_bar=scalars_bar_view[i],
                                  clim=[cmin, cmax], name=fieldName)
            i += 1

        self.plotter.remove_actor("speaker")
        self.plotter.add_mesh(self.mesh, scalars=normalized_system_pressure,
                              cmap=cmap, n_colors=ncolor, show_scalar_bar=False, clim=[cmin, cmax], name="speaker")

    def getMinMaxPressure(self, fieldName, fieldPressure):
        maxPressure_ABS = np.max(np.abs(fieldPressure))
        maxPressure_SPL = np.max(gtb.gain.SPL(fieldPressure))
        maxPressure_REAL = np.max(np.real(fieldPressure))
        maxPressure_IMAG = np.max(np.imag(fieldPressure))
        max_dict = {"abs": maxPressure_ABS, "spl": maxPressure_SPL,
                    "real": maxPressure_REAL, "imag": maxPressure_IMAG}
        self.maxPressure[fieldName] = max_dict

    def updateMinMaxScalars(self, typeScalar):
        if typeScalar == "spl":
            maxSPL = []
            for fieldName in self.fields:
                maxSPL.append(self.maxPressure[fieldName]["spl"])
            MAX = np.max(maxSPL)
        elif typeScalar == "abs":
            maxSPL = []
            for fieldName in self.fields:
                maxSPL.append(self.maxPressure[fieldName]["abs"])
            MAX = np.max(maxSPL)
        elif typeScalar == "real":
            maxSPL = []
            for fieldName in self.fields:
                maxSPL.append(self.maxPressure[fieldName]["real"])
            MAX = np.max(maxSPL)
        return MAX
