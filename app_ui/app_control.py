# Control modules
import electroacPy as ep
from PyQt6.QtGui import QAction, QIcon

from control_panel.driver_tab import DriverTab
from control_panel.enclosure_tab import EnclosureTab
from control_panel.CAD_tab import CADTab
from control_panel.study_tab import StudyTab

# Visualization modules
from visualization_panel.CADView_tab import CADView
from visualization_panel.MeshView_tab import MeshView
from visualization_panel.ResultView_tab import ResultView

# Info modules
from info_panel.driver_info_tab import InfoTabDriver
from info_panel.enclosure_info_tab import InfoTabEnclosure
from info_panel.result_info_tab import InfoResultTab

# UI
import logging
import sys
from os.path import join
from pathlib import Path
from PyQt6.QtCore import (QRect, QTimer)
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QTabWidget,
                             QPushButton, QListWidget, QVBoxLayout, QLabel,
                             QLineEdit, QMenuBar, QMenu, QFileDialog, QComboBox)

# helpful tool
from numpy import savez, load
from shutil import copy2

logging.basicConfig(level=logging.WARNING)
numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)
bempp_logger = logging.getLogger('bempp')
bempp_logger.setLevel(logging.WARNING)
opencl_logger = logging.getLogger('pyopencl')
opencl_logger.setLevel(logging.WARNING)
matplotlib_logger = logging.getLogger("matplotlib")
matplotlib_logger.setLevel(logging.WARNING)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.system = ep.loudspeakerSystem()
        self.geo = {'cad_changed': False,
                    'cad_view_path': None,
                    'meshCAD': None,
                    'meshPath': None,
                    'tmp_faces': None,
                    'minSize': None,
                    'maxSize': None,
                    'ppw': None,
                    'surfaceList': {},
                    'surfaceSize': {},
                    'surfaceGroup': {},
                    'update_mesh_display': None}

        self.evaluation_setup = {"boundaries": [],
                                 "radiators": None,
                                 "offsets": [],
                                 "results": {}}

        self.filters = {}  # store filters as individual dictionaries

        self.setWindowTitle("ElectroacPy - GUI")
        self.setGeometry(100, 100, 1100, 600 + 55)

        # Create the central widget
        self.central_widget = QWidget(self)
        self.central_widget.setGeometry(QRect(0, 0, 1100, 600 + 50))

        # Create control tab
        self.control_tab = ControlPanel(self.central_widget, self.system, self.geo,
                                        self.evaluation_setup)

        # Create view tab
        self.view_tab = ViewPanel(self.central_widget, self.system, self.geo, self.evaluation_setup)


        # Create Info tab
        self.info_tab = InfoPanel(self.central_widget, self.system, self.geo, self.evaluation_setup,
                                  self.control_tab, self.view_tab)

        # create menubar
        self.menu_bar = MenuBar(self.central_widget, self.system, self.geo, self.evaluation_setup,
                                self.control_tab, self.view_tab, self.info_tab)


class MenuBar:
    def __init__(self, central_widget, system, geo, evaluation, control_tab, view_tab, info_tab):
        self.central_widget = central_widget
        self.system = system
        self.geo = geo
        self.evaluation = evaluation
        self.control_tab = control_tab
        self.view_tab = view_tab
        self.info_tab = info_tab
        self.menubar = QMenuBar(central_widget)
        self.menubar.setGeometry(QRect(0, 0, 810 + 290, 22))
        self.simulation_parameter_window = None

        # DEFINE ACTIONS
        button_open = QAction(QIcon(join("icons", "folder_50.png")), "&Open", central_widget)
        button_open.triggered.connect(self.actionLoad)
        button_save = QAction(QIcon(join("icons", "save_50.png")), "&Save", central_widget)
        button_save.triggered.connect(self.actionSave)
        button_export = QAction(QIcon(join("icons", "export_50.png")), "&Export", central_widget)
        button_export.triggered.connect(self.actionExport)

        button_frequencies = QAction("&Frequencies", central_widget)
        button_frequencies.triggered.connect(self.actionFrequencies)

        # LOAD / SAVE / EXPORT FILES
        self.menuFile = QMenu("&File", self.menubar)
        self.menuFile.addAction(button_open)
        self.menuFile.addAction(button_save)
        self.menuFile.addAction(button_export)
        self.menubar.addMenu(self.menuFile)

        # SIMULATION PARAMETERS
        self.menuSimulation_Parameters = QMenu("&Simulation Parameters", self.menubar)
        self.menuSimulation_Parameters.addAction(button_frequencies)
        self.menubar.addMenu(self.menuSimulation_Parameters)

    def actionLoad(self):
        dir_path = QFileDialog.getExistingDirectory(parent=self.central_widget,
                                                    caption="Select directory",
                                                    directory=str(Path.home()),
                                                    options=QFileDialog.Option.DontUseNativeDialog)
        if dir_path:
            from io_ui.load_system import load_system
            load_system(dir_path, self.control_tab, self.geo, self.system, self.evaluation)

    def actionSave(self):
        dir_path = QFileDialog.getExistingDirectory(parent=self.central_widget,
                                                    caption="Select directory",
                                                    directory=str(Path.home()),
                                                    options=QFileDialog.Option.DontUseNativeDialog)
        # save mesh
        if self.geo["meshPath"] is not None:
            copy2(self.geo["meshPath"], dir_path)
            self.geo["meshPath"] = dir_path + self.geo['meshPath']

        geo = self.geo
        savez(join(dir_path, "geometry.npz"), geo=geo)

        evaluation = self.evaluation
        savez(join(dir_path, "evaluation_setup.npz"), evaluation=evaluation)
        ep.save(dir_path, self.system)

    def actionExport(self):
        pass

    def actionFrequencies(self):
        if self.simulation_parameter_window is None:
            self.simulation_parameter_window = SimulationParameterWindow(self.system, self.view_tab, self.info_tab, self.evaluation)
        self.simulation_parameter_window.show()


class SimulationParameterWindow(QWidget):
    def __init__(self, system, view_tab, info_tab, evaluation_setup):
        super().__init__()
        self.setWindowTitle("Simulation parameters")
        self.system = system
        self.view_tab = view_tab
        self.info_tab = info_tab
        self.evaluation_setup = evaluation_setup

        self.labels = {}
        text = ["Frequency range",
                "start",
                "stop",
                "Nfft",
                "scale",
                "Propagation Medium",
                "c (m/s)",
                "rho (kg/m^3)"]

        dx = [20, 20, 20, 20, 20, 160, 160, 160]
        dy = [20, 50, 80, 110, 140, 20, 50, 80]
        lx = [101, 31, 31, 31, 31, 121, 41, 71]
        for i in range(len(dx)):
            self.labels[text[i]] = QLabel(text[i], self)
            self.labels[text[i]].setGeometry(QRect(dx[i], dy[i], lx[i], 16))

        self.label_warning = QLabel("Updating parameters will reset computed studies!", self)
        self.label_warning.setGeometry(QRect(20, 170, 311, 16))

        self.lineEdit_start = QLineEdit("20", self)
        self.lineEdit_start.setGeometry(QRect(60, 50, 61, 21))

        self.lineEdit_stop = QLineEdit("2500", self)
        self.lineEdit_stop.setGeometry(QRect(60, 80, 61, 21))

        self.lineEdit_nfft = QLineEdit("50", self)
        self.lineEdit_nfft.setGeometry(QRect(60, 110, 61, 21))

        self.lineEdit_c = QLineEdit("343", self)
        self.lineEdit_c.setGeometry(QRect(240, 50, 61, 21))

        self.lineEdit_rho = QLineEdit("1.22", self)
        self.lineEdit_rho.setGeometry(QRect(240, 80, 61, 21))

        # pushbutton
        self.pushbutton_update = QPushButton("Update", self)
        self.pushbutton_update.setGeometry(QRect(220, 140, 75, 24))
        self.pushbutton_update.clicked.connect(self.update_system)

        # scale
        self.combobox_fscale = QComboBox(self)
        self.combobox_fscale.setGeometry(QRect(60, 140, 61, 22))
        self.combobox_fscale.addItem("log")
        self.combobox_fscale.addItem("linear")
        self.combobox_fscale.setCurrentIndex(0)

    def update_system(self):
        import numpy as np
        from generalToolbox.freqop import freq_log10
        from app_ui.update_tools.update_driver import update_driver
        from app_ui.update_tools.update_enclosure import update_enclosure
        from app_ui.update_tools.update_study import update_study

        fmin = int(self.lineEdit_start.text())
        fmax = int(self.lineEdit_stop.text())
        nfft = int(self.lineEdit_nfft.text())
        scale = self.combobox_fscale.currentText()

        rho = float(self.lineEdit_rho.text())
        c = float(self.lineEdit_c.text())

        if scale == "linear":
            frequency = np.arange(fmin, fmax, nfft)
        elif scale == "log":
            frequency = freq_log10(fmin, fmax, nfft)

        # update global parameters
        self.system.frequency = frequency
        self.system.c = c
        self.system.rho = rho

        # update objects
        update_driver(self.system)
        update_enclosure(self.system)
        update_study(self.system, self.view_tab, self.info_tab, self.evaluation_setup)


class ControlPanel:
    def __init__(self, central_widget, system, geo, evaluation_setup):
        self.tabs = QTabWidget(central_widget)
        self.tabs.setGeometry(QRect(0, 23, 290, 360))

        # control tabs -> should follow building order: mesh > driver > enclosure > study > filters
        self.cad_tab = CADTab(self.tabs, system, geo)
        self.driver_tab = DriverTab(self.tabs, system)
        self.enclosure_tab = EnclosureTab(self.tabs, system, self.driver_tab)
        self.study_tab = StudyTab(self.tabs, system, geo, evaluation_setup,
                                  self.driver_tab, self.enclosure_tab)


class InfoPanel:
    def __init__(self, central_widget, system, geo, evaluation_setup, control_panel, view_panel):
        self.tabs = QTabWidget(central_widget)
        self.tabs.setGeometry(QRect(0, 360 + 23, 291, 261))
        self.tabs.TabPosition(QTabWidget.TabPosition.South)

        self.driver_tab_info = InfoTabDriver(self.tabs, control_panel.driver_tab, system)
        self.enclosure_tab_info = InfoTabEnclosure(self.tabs, control_panel.enclosure_tab, system)

        self.result_tab_info = InfoResultTab(self.tabs, view_panel, control_panel, system, geo, evaluation_setup)


class ViewPanel:
    def __init__(self, central_widget, system, geo, evaluation_setup):
        self.system = system
        self.geo = geo
        self.evaluation_setup = evaluation_setup

        self.tabs = QTabWidget(central_widget)
        self.tabs.setGeometry(QRect(290, 23, 810, 600 + 50))
        self.cadView_tab = CADView(self.tabs, system, geo)
        self.tabs.count()
        # check for mesh creation / update
        self.timer_mesh_detect = QTimer(self.tabs)
        self.timer_mesh_detect.timeout.connect(self.check_for_mesh_update)
        self.timer_mesh_detect.start(1000)  # Check every second
        self.meshView_tab = None  # tab to display mesh

        # check for computed observations
        self.timer_observation_detect = QTimer(self.tabs)
        self.timer_observation_detect.timeout.connect(self.check_for_computed_observation)
        self.timer_observation_detect.start(1000)
        self.obsView_tabs = {"FIELDS_VIEWER": {"name": [],
                                               "plotter": None}}

    def check_for_mesh_update(self):
        update_mesh = self.geo["update_mesh_display"]
        # print(self.geo)
        if update_mesh is True:
            self.update_mesh_view()
            self.geo["update_mesh_display"] = False
        return None

    def update_mesh_view(self):
        if self.meshView_tab is None:
            self.meshView_tab = MeshView(self.tabs, self.system, self.geo)
        else:
            self.meshView_tab.updateMesh()
        return None

    def check_for_computed_observation(self):
        # if self.system.observation:
        #     if self.system.observation["GUI"].observationName is not False:  # check is study is initialized
        #         for obs in self.system.observation["GUI"].computedObservations:
        #             self.obsView_tabs[obs] = ResultView(self.tabs, self.evaluation_setup, None, self.geo)
        if self.meshView_tab is not None and bool(self.evaluation_setup["results"]) is not False:
            # plot observations after mesh
            for obs in self.evaluation_setup["results"]:
                if obs in self.obsView_tabs or obs in self.obsView_tabs["FIELDS_VIEWER"]["name"]:
                    pass
                else:
                    if self.evaluation_setup[obs]["type"] == "polar" or self.evaluation_setup[obs][
                        "type"] == "point":
                        self.obsView_tabs[obs] = ResultView(self.tabs, obs, self.evaluation_setup,
                                                            False, self.geo, self.system.frequency,
                                                            self.system)
                    elif self.evaluation_setup[obs]["type"] == "screen":
                        if self.obsView_tabs["FIELDS_VIEWER"]["plotter"] is None:
                            self.obsView_tabs["FIELDS_VIEWER"]["plotter"] = ResultView(self.tabs, obs,
                                                                                       self.evaluation_setup,
                                                                                       False, self.geo,
                                                                                       self.system.frequency,
                                                                                       self.system)
                            self.obsView_tabs["FIELDS_VIEWER"]["name"].append(obs)
                        else:
                            self.obsView_tabs["FIELDS_VIEWER"]["plotter"].addField(obs, self.evaluation_setup)
                            self.obsView_tabs["FIELDS_VIEWER"]["name"].append(obs)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
