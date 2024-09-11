from PyQt6.QtCore import QRect, QTimer
from PyQt6.QtWidgets import QWidget, QLabel, QLineEdit, QPushButton, QListWidget, QComboBox, QFrame, QFileDialog, \
    QAbstractItemView, QPlainTextEdit, QRadioButton

import numpy as np


class InfoResultTab:
    def __init__(self, info_panel, view_panel, control_panel, system, geo, evaluation_setup):
        self.info_panel = info_panel
        self.view_panel = view_panel
        self.control_panel = control_panel
        self.system = system
        self.geo = geo
        self.evaluation_setup = evaluation_setup
        self.plot_parameters = {"polar": {},
                                "fields": {}}

        # additional windows
        self.polar_parameter_window = None
        self.field_parameter_window = None

        self.tab = QWidget()

        # labels
        self.label_radGroupSelect = QLabel("Radiating groups", self.tab)
        self.label_radGroupSelect.setGeometry(QRect(10, 20, 101, 16))

        self.label_fieldFreqSelect = QLabel("Fields frequency", self.tab)
        self.label_fieldFreqSelect.setGeometry(QRect(10, 140, 131, 16))

        # self.label_fieldFreqUnit = QLabel("Hz", self.tab)
        # self.label_fieldFreqUnit.setGeometry(QRect(130, 160, 31, 16))

        # pushbutton
        self.pushbutton_updatePlots = QPushButton("Update plots", self.tab)
        self.pushbutton_updatePlots.setGeometry(QRect(200, 200, 81, 24))
        self.pushbutton_updatePlots.clicked.connect(self.updatePlot)

        self.pushbutton_polarParameters = QPushButton("Polar parameters", self.tab)
        self.pushbutton_polarParameters.setGeometry(QRect(10, 160, 101, 24))
        self.pushbutton_polarParameters.clicked.connect(self.setPolarParameters)

        self.pushbutton_fieldParameters = QPushButton("Field parameters", self.tab)
        self.pushbutton_fieldParameters.setGeometry(QRect(10, 190, 101, 24))
        self.pushbutton_fieldParameters.clicked.connect(self.setFieldParameters)

        # listWidget
        self.listWidget = QListWidget(self.tab)
        self.listWidget.setGeometry(QRect(10, 50, 141, 61))
        self.listWidget.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)

        # lineEdit
        # self.lineEdit_freq = QLineEdit(self.tab)
        # self.lineEdit_freq.setGeometry(QRect(10, 160, 113, 21))

        # update eval tab
        self.isTabInitialized = False
        self.timer_eval_detect = QTimer(self.tab)
        self.timer_eval_detect.timeout.connect(self.buildInfoTab)
        self.timer_eval_detect.start(1000)  # Check every second

        info_panel.addTab(self.tab, "Results")

    def buildInfoTab(self):
        # check if studies are computed and if the tab is not already built
        if bool(self.evaluation_setup["results"]) is True and self.isTabInitialized is False:
            for rad in self.evaluation_setup['radiators']:
                item2add = ""
                if self.system.radiator_id[rad] == "SPKBOX":
                    ref2bem = self.system.enclosure[rad].ref2bem
                    for i in range(len(ref2bem)):
                        item2add = "{} - ".format(ref2bem[i]) + rad
                        self.listWidget.addItem(item2add)
                elif self.system.radiator_id[rad] == "EAD":
                    ref2bem = self.system.driver[rad].ref2bem
                    for i in range(len(ref2bem)):
                        item2add = "{} - ".format(ref2bem[i]) + rad
                        self.listWidget.addItem(item2add)
                elif self.system.radiator_id[rad] == "PLV":
                    ref2bem = self.system.laser_acc[rad].ref2bem
                    for i in range(len(ref2bem)):
                        item2add = "{} - ".format(ref2bem[i]) + rad
                        self.listWidget.addItem(item2add)
            self.isTabInitialized = True

    def updatePlot(self, type=None):
        if type is None:  # UPDATE ALL PLOTS
            for name in self.evaluation_setup: # UPDATE POLAR RADIATIONS
                if (name == "boundaries" or name == "radiators" or name == "offsets"
                        or name == "results"):
                    pass
                else:
                    selectedSurfaces = []
                    for i in range(len(self.listWidget.selectedItems())):
                        selectedSurfaces.append(extract_number(self.listWidget.selectedItems()[i].text()))
                    # print("SELECTED SURFACES: ", selectedSurfaces)
                    if name not in self.view_panel.obsView_tabs["FIELDS_VIEWER"]["name"]:
                        self.view_panel.obsView_tabs[name].updatePlot(selectedSurfaces, self.plot_parameters)

            if bool(self.view_panel.obsView_tabs["FIELDS_VIEWER"]["name"]) is True:  # UPDATE FIELDS VIEWER
                self.view_panel.obsView_tabs["FIELDS_VIEWER"]["plotter"].updatePlot(selectedSurfaces,
                                                                                    self.plot_parameters)



    def setPolarParameters(self):
        if self.polar_parameter_window is None:
            self.polar_parameter_window = PolarParameterWindow(self, self.plot_parameters)
        self.polar_parameter_window.show()
        return None

    def setFieldParameters(self):
        if self.field_parameter_window is None:
            self.field_parameter_window = FieldParameterWindow(self, self.plot_parameters)
        self.field_parameter_window.show()
        return None


def extract_number(input_string):
    # Split the string at the first occurrence of " - "
    parts = input_string.split(" - ", 1)
    # Return the part before the "-", converted to an integer
    return int(parts[0])


class PolarParameterWindow(QWidget):
    def __init__(self, parent, plot_parameters):
        super().__init__()
        self.param_update = plot_parameters
        self.setWindowTitle("Polar observation parameters")
        self.parent = parent

        # SET ALL LABELS
        labels = ["Polar observation parameters",
                  "Frequency axis",
                  "f min", "f max",
                  "Level axis",
                  "max dB", "min dB", "step", "colormap",
                  "Frequency response at angle",
                  "Polar response at frequency"]
        self.labels = {}
        dx = [10, 30, 130, 130, 30, 30, 30, 30, 140, 30, 30]
        dy = [10, 40, 70, 100, 150, 180, 210, 240, 150, 280, 410]
        lx = [161, 91, 31, 31, 61, 49, 49, 49, 61, 161, 161]

        for i, label in enumerate(labels):
            self.labels[label] = QLabel(label, self)
            self.labels[label].setGeometry(QRect(dx[i], dy[i], lx[i], 16))

        # radio button
        self.radiobutton_logx = QRadioButton("log", self)
        self.radiobutton_logx.setGeometry(QRect(30, 70, 89, 20))
        self.radiobutton_logx.toggle()

        self.radiobutton_linear = QRadioButton("linear", self)
        self.radiobutton_linear.setGeometry(QRect(30, 100, 89, 20))

        # lineEdit
        self.lineEdit_f_min = QLineEdit(str(int(parent.system.frequency[0])), self)
        self.lineEdit_f_min.setGeometry(QRect(170, 70, 41, 21))

        self.lineEdit_f_max = QLineEdit(str(int(parent.system.frequency[-1])), self)
        self.lineEdit_f_max.setGeometry(QRect(170, 100, 41, 21))

        self.lineEdit_dB_max = QLineEdit("100", self)
        self.lineEdit_dB_max.setGeometry(QRect(80, 180, 31, 21))

        self.lineEdit_dB_min = QLineEdit("10", self)
        self.lineEdit_dB_min.setGeometry(QRect(80, 210, 31, 21))

        self.lineEdit_step = QLineEdit("6", self)
        self.lineEdit_step.setGeometry(QRect(80, 240, 31, 21))

        # combobox
        self.combobox_cmap = QComboBox(self)
        self.combobox_cmap.setGeometry(QRect(140, 170, 71, 22))
        self.combobox_cmap.addItem("turbo")
        self.combobox_cmap.addItem("viridis")
        self.combobox_cmap.addItem("jet")
        self.combobox_cmap.addItem("seismic")
        self.combobox_cmap.addItem("magma")

        # pushbutton
        self.pushbutton_plot = QPushButton("Update plots", self)
        self.pushbutton_plot.setGeometry(QRect(30, 540, 81, 24))
        self.pushbutton_plot.clicked.connect(self.sendNewParam)

        # listWidget
        self.listWidget_angle = QListWidget(self)
        self.listWidget_angle.setGeometry(QRect(30, 310, 141, 91))
        self.listWidget_angle.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        theta = np.arange(-180, 181, 1)
        for idx, valt in enumerate(theta):  # implementation is a bit stupid
            self.listWidget_angle.addItem(str(valt))

        self.listWidget_frequency = QListWidget(self)
        self.listWidget_frequency.setGeometry(QRect(30, 440, 141, 91))
        for idx, valf in enumerate(parent.system.frequency):
            self.listWidget_frequency.addItem(str(round(valf, 3)))
        self.listWidget_frequency.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)

    def sendNewParam(self):
        self.param_update["polar"] = {"dBmin": float(self.lineEdit_dB_min.text()),
                                      "dBmax": float(self.lineEdit_dB_max.text()),
                                      "dBstep": float(self.lineEdit_step.text()),
                                      "cmap": self.combobox_cmap.currentText(),
                                      "fmin": float(self.lineEdit_f_min.text()),
                                      "fmax": float(self.lineEdit_f_max.text())}
        if self.radiobutton_linear.isChecked():
            self.param_update["polar"]["xscale"] = "linear"
        elif self.radiobutton_logx.isChecked():
            self.param_update["polar"]["xscale"] = "log"

        selectedFrequencies = []
        for i in range(len(self.listWidget_frequency.selectedItems())):
            selectedFrequencies.append(float(self.listWidget_frequency.selectedItems()[i].text()))
        self.param_update["polar"]["freq2plot"] = selectedFrequencies

        selectedAngles = []
        for i in range(len(self.listWidget_angle.selectedItems())):
            selectedAngles.append(float(self.listWidget_angle.selectedItems()[i].text()))
        self.param_update["polar"]["angle2plot"] = selectedAngles

        self.parent.updatePlot()


class FieldParameterWindow(QWidget):
    def __init__(self, parent, plot_parameters):
        super().__init__()
        self.param_update = plot_parameters
        self.setWindowTitle("Fields parameters")
        self.parent = parent

        # SET ALL LABELS
        labels = ["Pressure field parameters",
                  "Level axis",
                  "max Lvl", "min Lvl", "n_colors", "colormap",
                  "Frequency to plot"]
        self.labels = {}
        dx = [10, 20, 20, 20, 20, 130, 10]
        dy = [10, 40, 70, 100, 130, 40, 170]
        lx = [161, 91, 49, 49, 61, 61, 101]

        for i, label in enumerate(labels):
            self.labels[label] = QLabel(label, self)
            self.labels[label].setGeometry(QRect(dx[i], dy[i], lx[i], 16))

        # lineEdit
        self.lineEdit_dB_max = QLineEdit("120", self)
        self.lineEdit_dB_max.setGeometry(QRect(70, 70, 31, 21))

        self.lineEdit_dB_min = QLineEdit("30", self)
        self.lineEdit_dB_min.setGeometry(QRect(70, 100, 31, 21))

        self.lineEdit_step = QLineEdit("21", self)
        self.lineEdit_step.setGeometry(QRect(70, 130, 31, 21))

        # combobox - colormap
        self.combobox_cmap = QComboBox(self)
        self.combobox_cmap.setGeometry(QRect(130, 60, 71, 22))
        self.combobox_cmap.addItem("turbo")
        self.combobox_cmap.addItem("viridis")
        self.combobox_cmap.addItem("jet")
        self.combobox_cmap.addItem("seismic")
        self.combobox_cmap.addItem("magma")


        # combobox - transformation
        self.combobox_transform = QComboBox(self)
        self.combobox_transform.setGeometry(QRect(130, 120, 71, 22))
        self.combobox_transform.addItem("SPL")
        self.combobox_transform.addItem("Real")
        self.combobox_transform.addItem("Phase")

        # listWidget
        self.listWidget_frequency = QListWidget(self)
        self.listWidget_frequency.setGeometry(QRect(20, 190, 121, 351))

        for idx, valf in enumerate(parent.system.frequency):
            self.listWidget_frequency.addItem(str(round(valf, 3)))

        # pushbutton
        self.pushbutton_plot = QPushButton("Update plots", self)
        self.pushbutton_plot.setGeometry(QRect(154, 540, 81, 24))
        self.pushbutton_plot.clicked.connect(self.sendNewParam)

    def sendNewParam(self):
        self.param_update["fields"] = {"dBmin": float(self.lineEdit_dB_min.text()),
                                       "dBmax": float(self.lineEdit_dB_max.text()),
                                       "n_colors": int(self.lineEdit_step.text()),
                                       "cmap": self.combobox_cmap.currentText(),
                                       "transformation": self.combobox_transform.currentText(),
                                       "ncolor": 21}

        self.param_update["fields"]["freq2plot"] = float(self.listWidget_frequency.currentItem().text())
        self.parent.updatePlot()
