from PyQt6.QtCore import QRect
from PyQt6.QtWidgets import (QWidget, QLabel, QLineEdit, QPushButton, QListWidget, QComboBox,
                             QCheckBox, QRadioButton, QGraphicsView, QTableWidget, QTableWidgetItem)

# for visualization
from OCC.Display.backend import load_backend

load_backend("pyqt6")

from OCC.Display import qtDisplay
from OCC.Core.Quantity import (Quantity_NOC_ALICEBLUE,
                               Quantity_NOC_WHITE,
                               Quantity_NOC_WHITE,
                               Quantity_NOC_GRAY40,
                               Quantity_NOC_BLACK,
                               Quantity_NOC_RED,
                               Quantity_NOC_BLUE1)

from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Ax2
from OCC.Core.Geom import Geom_CartesianPoint, Geom_Line
from OCC.Core.Quantity import Quantity_Color
from OCC.Core.Aspect import Aspect_TOM_O
from OCC.Core.AIS import AIS_Point, AIS_Axis
from OCC.Core.Prs3d import Prs3d_PointAspect
from OCC.Extend.DataExchange import read_step_file
from OCC.Extend.TopologyUtils import TopologyExplorer

import numpy as np


class StudyTab:
    def __init__(self, control_panel, system, geo,
                 evaluation_setup, driver_tab, enclosure_tab):
        self.system = system
        self.geo = geo
        self.driverList = driver_tab.listWidget
        self.enclosureList = enclosure_tab.listWidget
        # self.results = results

        self.tab = QWidget()
        self.evaluation_setup = evaluation_setup

        # extra setup windows
        self.setup_window = None

        # Labels
        self.label_setup = QLabel(self.tab)
        self.label_setup.setGeometry(QRect(10, 20, 111, 22))
        self.label_setup.setText("Evaluation setup")

        self.label_boundaries = QLabel(self.tab)
        self.label_boundaries.setGeometry(QRect(10, 190, 121, 16))
        self.label_boundaries.setText("Infinite boundaries")

        self.label_plane = QLabel(self.tab)
        self.label_plane.setGeometry(QRect(10, 210, 49, 16))
        self.label_plane.setText("Plane")

        self.label_offset = QLabel(self.tab)
        self.label_offset.setGeometry(QRect(70, 210, 51, 16))
        self.label_offset.setText("Offset")

        planes = ["xy", "yz", "zx"]
        planes_z = [230, 260, 290]


        for i in range(3):
            self.label_offset_xyz = QLabel(self.tab)
            self.label_offset_xyz.setGeometry(QRect(110, planes_z[i], 21, 21))
            self.label_offset_xyz.setText("m")

        # checkbox for boundary position
        self.check_plane_xy = QCheckBox(self.tab)
        self.check_plane_xy.setGeometry(QRect(10, 230, 75, 20))
        self.check_plane_xy.setText(planes[0])

        self.check_plane_yz = QCheckBox(self.tab)
        self.check_plane_yz.setGeometry(QRect(10, 260, 75, 20))
        self.check_plane_yz.setText(planes[1])

        self.check_plane_zx = QCheckBox(self.tab)
        self.check_plane_zx.setGeometry(QRect(10, 290, 75, 20))
        self.check_plane_zx.setText(planes[2])

        self.lineEdit_offset_xy = QLineEdit(self.tab)
        self.lineEdit_offset_xy.setGeometry(QRect(70, planes_z[0], 31, 20))
        self.lineEdit_offset_xy.setText("0")

        self.lineEdit_offset_yz = QLineEdit(self.tab)
        self.lineEdit_offset_yz.setGeometry(QRect(70, planes_z[1], 31, 20))
        self.lineEdit_offset_yz.setText("0")

        self.lineEdit_offset_zx = QLineEdit(self.tab)
        self.lineEdit_offset_zx.setGeometry(QRect(70, planes_z[2], 31, 20))
        self.lineEdit_offset_zx.setText("0")

        # pushbutton
        # ADD EVAL TO SYSTEM
        self.pushbutton_addEvaluation = QPushButton(self.tab)
        self.pushbutton_addEvaluation.setGeometry(QRect(180, 40, 91, 24))
        self.pushbutton_addEvaluation.setText("Add evaluation")
        self.pushbutton_addEvaluation.clicked.connect(self.add_observation)

        # REMOVE EVAL FROM SYSTEM
        self.pushbutton_removeSelection = QPushButton(self.tab)
        self.pushbutton_removeSelection.setGeometry(QRect(10, 160, 111, 24))
        self.pushbutton_removeSelection.setText("Remove selection")
        self.pushbutton_removeSelection.clicked.connect(self.remove_selection)

        # MODIFY SELECTED EVAL
        self.pushbutton_modifySelection = QPushButton(self.tab)
        self.pushbutton_modifySelection.setGeometry(QRect(160, 160, 111, 24))
        self.pushbutton_modifySelection.setText("Modify selection")
        self.pushbutton_modifySelection.clicked.connect(self.modify_selection)

        # RUN FULL STUDY
        self.pushbutton_run = QPushButton(self.tab)
        self.pushbutton_run.setGeometry(QRect(160, 270, 111, 24))
        self.pushbutton_run.setText("Run full study")
        self.pushbutton_run.clicked.connect(self.run_study)

        # RUN OBSERVATION ONLY
        self.pushbutton_update_eval = QPushButton(self.tab)
        self.pushbutton_update_eval.setGeometry(QRect(160, 300, 111, 24))
        self.pushbutton_update_eval.setText("Update evaluations")
        self.pushbutton_update_eval.clicked.connect(self.run_observation)

        # QlistWidget
        self.list_eval = QListWidget(self.tab)
        self.list_eval.setGeometry(QRect(10, 70, 256, 81))

        # ComboBox for evaluation selection
        self.combo_eval = QComboBox(self.tab)
        self.combo_eval.setGeometry(QRect(10, 40, 111, 22))
        items = ["Polar", "Field", "Point"] #, "Sphere"]
        for i in range(len(items)):
            self.combo_eval.addItem(items[i])

        control_panel.addTab(self.tab, "Study")

    def add_observation(self):
        currentItem = self.combo_eval.currentText()

        # POLAR RADIATION
        if currentItem == "Polar":
            if self.setup_window is None:
                self.setup_window = addPolarArrayWindow(self.system, self.geo,
                                                        self.evaluation_setup, self.list_eval)
            else:
                self.setup_window.close()
                # self.polar_window_setup = None
                self.setup_window = addPolarArrayWindow(self.system, self.geo,
                                                        self.evaluation_setup, self.list_eval)
            self.setup_window.show()

        # PRESSURE FIELD
        elif currentItem == "Field":
            if self.setup_window is None:
                self.setup_window = addPressureFieldWindow(self.system, self.geo,
                                                           self.evaluation_setup, self.list_eval)
            else:
                self.setup_window.close()
                self.setup_window = addPressureFieldWindow(self.system, self.geo,
                                                           self.evaluation_setup, self.list_eval)
            self.setup_window.show()
        # PRESSURE POINT
        elif currentItem == "Point":
            if self.setup_window is None:
                self.setup_window = addPressurePointWindow(self.system, self.geo,
                                                           self.evaluation_setup, self.list_eval)
            else:
                self.setup_window.close()
                self.setup_window = addPressurePointWindow(self.system, self.geo,
                                                           self.evaluation_setup, self.list_eval)
            self.setup_window.show()

        else:
            pass
        return None

    def modify_selection(self):
        # TODO -> open window setup and modify observation parameters
        # TODO -> THIS PART COULD ADD SOME ISSUES TO THE SELF.EVALUATION_SETUP VARIABLE
        if self.list_eval.currentItem():
            if self.setup_window is not None:
                self.setup_window.close()
            if self.evaluation_setup[self.list_eval.currentItem().text()]["type"] == "polar":
                current_eval = self.list_eval.currentItem().text()
                self.setup_window = addPolarArrayWindow(self.system, self.geo,
                                                        self.evaluation_setup, self.list_eval)
                fill_polar_window(self.evaluation_setup, current_eval, self.setup_window)

                # remove previous eval from evaluation_list and QList widget
                self.evaluation_setup.pop(current_eval)
                self.list_eval.takeItem(self.list_eval.row(self.list_eval.currentItem()))
                self.setup_window.show()

            elif self.evaluation_setup[self.list_eval.currentItem().text()]["type"] == "screen":
                current_eval = self.list_eval.currentItem().text()
                self.setup_window = addPressureFieldWindow(self.system, self.geo,
                                                        self.evaluation_setup, self.list_eval)
                fill_pressure_field_window(self.evaluation_setup, current_eval, self.setup_window)

                # remove previous eval from evaluation_list and QList widget
                self.evaluation_setup.pop(current_eval)
                self.list_eval.takeItem(self.list_eval.row(self.list_eval.currentItem()))
                self.setup_window.show()

            elif self.evaluation_setup[self.list_eval.currentItem().text()]["type"] == "point":
                current_eval = self.list_eval.currentItem().text()
                self.setup_window = addPressurePointWindow(self.system, self.geo,
                                                        self.evaluation_setup, self.list_eval)
                fill_pressure_point_window(self.evaluation_setup, current_eval, self.setup_window)

                # remove previous eval from evaluation_list and QList widget
                self.evaluation_setup.pop(current_eval)
                self.list_eval.takeItem(self.list_eval.row(self.list_eval.currentItem()))
                self.setup_window.show()
        pass

    def remove_selection(self):
        if self.list_eval.currentItem():
            self.evaluation_setup.pop(self.list_eval.currentItem().text())
            self.list_eval.takeItem(self.list_eval.row(self.list_eval.currentItem()))
        return None

    def run_study(self):
        self.build_study()
        self.system.run()
        self.store_results()


    def run_observation(self):
        pass

    def store_results(self):
        """
        Stores the computed evaluations into the self.evaluation_setup["results"] dictionary.
        self.evaluation_setup["results"].__dict__ should return as many sub-dictionaries as computed observations.
        :return:
        """

        for i, name in enumerate(self.evaluation_setup):
            if name == "boundaries" or name == "radiators" or name == "offsets" or name == "results":
                pass
            else:
                self.evaluation_setup["results"][name] = self.system.observation["GUI"].pMicArray[i-4]  # -> ndarray
                print(self.evaluation_setup["results"][name].shape)

    def build_study(self):
        acoustic_radiator = []
        boundary = []
        offsets = []

        # Add acoustic radiators to evaluation setup
        for drv in self.system.driver:
            if self.system.driver[drv].ref2bem is not None:
                acoustic_radiator.append(drv)

        for box in self.system.enclosure:
            if self.system.enclosure[box] is not None:
                acoustic_radiator.append(box)

        self.evaluation_setup["radiators"] = acoustic_radiator
        # print("ACOUSTIC RADIATORS: ", acoustic_radiator)


        # Add boundaries to evaluation setup
        if self.check_plane_xy.isChecked():
            boundary.append("z")
            offsets.append(float(self.lineEdit_offset_xy.text()))
        if self.check_plane_yz.isChecked():
            boundary.append("x")
            offsets.append(float(self.lineEdit_offset_yz.text()))
        if self.check_plane_zx.isChecked():
            boundary.append("y")
            offsets.append(float(self.lineEdit_offset_zx.text()))
        if bool(boundary) is False:
            boundary = False
            offsets = False

        self.evaluation_setup["boundaries"] = boundary
        self.evaluation_setup["offsets"] = offsets

        # print("BOUNDARIES: ", boundary)
        # print("OFFSETS: ", offsets)
        # print("EVALUATION SETUP: ", self.evaluation_setup)
        # BUIlD STUDY
        if self.geo["meshPath"]:
            self.system.study_exteriorBEM("GUI", self.geo["meshPath"], acoustic_radiator,
                                          boundary=boundary, offset=offsets)

            # ADD OBSERVATIONS TO STUDY IF AT LEAST 1 EVAL IS DEFINED
            if len(self.evaluation_setup) > 4:
                for eval in self.evaluation_setup:
                    print("EVALUATION: ", eval)
                    if eval == "boundaries" or eval == "radiators" or eval == "offsets" or eval == "results":
                        var = None
                    else:
                        # ADD POLAR RADIATION
                        if self.evaluation_setup[eval]["type"] == "polar":
                            minAngle = self.evaluation_setup[eval]["minAngle"]
                            maxAngle = self.evaluation_setup[eval]["maxAngle"]
                            step = self.evaluation_setup[eval]["step"]
                            # plane = self.evaluation_setup[eval]["plane"]
                            on_axis = self.evaluation_setup[eval]["on_axis"]
                            direction = self.evaluation_setup[eval]["direction"]
                            radius = self.evaluation_setup[eval]["radius"]
                            offset = self.evaluation_setup[eval]['offset']
                            self.system.observation_polarRadiation("GUI", eval,
                                                                   minAngle, maxAngle, step,
                                                                   on_axis, direction,
                                                                   radius, offset)
                        elif self.evaluation_setup[eval]["type"] == "screen":
                            L1 = self.evaluation_setup[eval]["L1"]
                            L2 = self.evaluation_setup[eval]["L2"]
                            step = self.evaluation_setup[eval]["step"]
                            offset = self.evaluation_setup[eval]["offset"]
                            plane = self.evaluation_setup[eval]["plane"]
                            self.system.observation_pressureField("GUI", eval, L1, L2, step, plane, offset)
                        elif self.evaluation_setup[eval]["type"] == "point":
                            xMic = self.evaluation_setup[eval]["xMic"]
                            labels = []
                            for i in range(len(xMic)):
                                labels.append("x: {}, y: {}, z: {}".format(xMic[i, 0], xMic[i, 1], xMic[i, 2]))
                            self.system.observation_pressureResponse("GUI", eval, xMic, labels=labels)
                # print(self.system.observation["GUI"].observationName)
                # print(self.system.observation["GUI"].xMic[0].shape)
                # print(self.system.observation["GUI"].xMic[1].shape)
                # self.system.plot_system("GUI")
        else:
            print("no surface mesh defined")
        return None

class addPolarArrayWindow(QWidget):
    def __init__(self, system, geo, evaluation_setup, evaluation_list):
        super().__init__()
        self.system = system
        self.geo = geo
        self.evaluation_setup = evaluation_setup
        self.evaluation_list = evaluation_list
        self.setWindowTitle("Polar Evaluation Setup")

        # plotter tab
        self.eval_view = QGraphicsView(self)
        self.eval_view.setGeometry(QRect(210, 20, 461, 371))

        # OpenCascade plotter
        self.plotter = qtDisplay.qtViewer3d(self.eval_view)  #QtInteractor(self.tab) #for pyvistaqt
        self.plotter.resize(461, 371)
        # self.plotter.InitDriver()
        self.display = self.plotter._display
        self.display.View.SetBgGradientColors(
            Quantity_Color(Quantity_NOC_ALICEBLUE),
            Quantity_Color(Quantity_NOC_WHITE),
            2,
            True,
        )
        self.display.Repaint()
        self.is_eval_displayed = False

        # param
        self.parameters = {"name": None,
                           "minAngle": None,
                           "maxAngle": None,
                           "step": None,
                           "radius": None,
                           "offset": [0, 0, 0],
                           "plane": None,
                           "on_axis": "x",
                           "direction": "y",
                           "type": "polar",
                           "toUpdate": False}

        # labels
        labels = ["Evaluation name", "Angle", "min", "max",
                  "step", "Radius (m)", "On-axis", "Offset", "x", "y", "z"]
        dimensions = [[20, 20, 91, 16], [20, 80, 49, 16], [20, 100, 31, 16],
                      [70, 100, 31, 16], [120, 100, 31, 16],
                      [20, 150, 61, 16], [20, 210, 49, 16], [20, 270, 49, 16],
                      [20, 290, 31, 16], [70, 290, 31, 16], [120, 290, 31, 16]]

        self.label = {}
        for index, label in enumerate(labels):
            self.label[label] = QLabel(label, self)
            x, y, dx, dy = (dimensions[index][0], dimensions[index][1],
                            dimensions[index][2], dimensions[index][3])
            self.label[label].setGeometry(QRect(x, y, dx, dy))

        # lineEdit
        self.lineEdit_name = QLineEdit(self)
        self.lineEdit_name.setGeometry(QRect(20, 40, 131, 21))

        self.lineEdit_minAngle = QLineEdit(self)
        self.lineEdit_minAngle.setGeometry(QRect(20, 120, 41, 21))
        self.lineEdit_minAngle.setText("-180")

        self.lineEdit_maxAngle = QLineEdit(self)
        self.lineEdit_maxAngle.setGeometry(QRect(70, 120, 41, 21))
        self.lineEdit_maxAngle.setText("180")

        self.lineEdit_step = QLineEdit(self)
        self.lineEdit_step.setGeometry(QRect(120, 120, 41, 21))
        self.lineEdit_step.setText("5")

        self.lineEdit_radius = QLineEdit(self)
        self.lineEdit_radius.setGeometry(QRect(20, 170, 41, 21))
        self.lineEdit_radius.setText("1.8")

        self.lineEdit_xOffset = QLineEdit(self)
        self.lineEdit_xOffset.setGeometry(QRect(20, 310, 41, 21))
        self.lineEdit_xOffset.setText("0")

        self.lineEdit_yOffset = QLineEdit(self)
        self.lineEdit_yOffset.setGeometry(QRect(70, 310, 41, 21))
        self.lineEdit_yOffset.setText("0")

        self.lineEdit_zOffset = QLineEdit(self)
        self.lineEdit_zOffset.setGeometry(QRect(120, 310, 41, 21))
        self.lineEdit_zOffset.setText("0")

        # radiobutton
        # self.radio_xy = QRadioButton("xy", self)
        # self.radio_xy.setGeometry(QRect(20, 220, 40, 20))
        # self.radio_xy.toggle()
        #
        # self.radio_yz = QRadioButton("yz", self)
        # self.radio_yz.setGeometry(QRect(70, 220, 40, 20))
        #
        # self.radio_zx = QRadioButton("zx", self)
        # self.radio_zx.setGeometry(QRect(120, 220, 40, 20))

        # combobox
        self.combo_onAxis = QComboBox(self)
        self.combo_onAxis.setGeometry(QRect(20, 230, 51, 22))
        self.combo_onAxis.addItem("x")
        self.combo_onAxis.addItem("-x")
        self.combo_onAxis.addItem("y")
        self.combo_onAxis.addItem("-y")
        self.combo_onAxis.addItem("z")
        self.combo_onAxis.addItem("-z")

        self.combo_direction = QComboBox(self)
        self.combo_direction.setGeometry(QRect(80, 230, 51, 22))
        self.combo_direction.addItem("x")
        self.combo_direction.addItem("-x")
        self.combo_direction.addItem("y")
        self.combo_direction.addItem("-y")
        self.combo_direction.addItem("z")
        self.combo_direction.addItem("-z")
        self.combo_direction.setCurrentIndex(2)
        # TODO: UPDATE THE CORRESPONDING ARRAY CREATION IN STUDY TAB AND INFO / PLOTTING WINDOW

        # checkbox
        # self.check_reverse = QCheckBox("Reverse on-axis", self)
        # self.check_reverse.setGeometry(QRect(20, 320, 101, 20))

        # pushbutton
        self.pushbutton_confirm = QPushButton("Confirm", self)
        self.pushbutton_confirm.setGeometry(QRect(110, 360, 75, 24))
        self.pushbutton_confirm.clicked.connect(self.confirm_eval)

        self.pushbutton_preview = QPushButton("Preview", self)
        self.pushbutton_preview.setGeometry(QRect(20, 360, 75, 24))
        self.pushbutton_preview.clicked.connect(self.update_evaluation)

    def update_evaluation(self):
        if self.lineEdit_name.text():
            self.parameters["name"] = self.lineEdit_name.text()
            self.parameters["minAngle"] = int(self.lineEdit_minAngle.text())
            self.parameters["maxAngle"] = int(self.lineEdit_maxAngle.text())
            self.parameters["step"] = float(self.lineEdit_step.text())
            self.parameters["radius"] = float(self.lineEdit_radius.text())
            xoff = float(self.lineEdit_xOffset.text())
            yoff = float(self.lineEdit_yOffset.text())
            zoff = float(self.lineEdit_zOffset.text())
            self.parameters["offset"] = [xoff, yoff, zoff]

            self.parameters["on_axis"] = self.combo_onAxis.currentText()
            self.parameters["direction"] = self.combo_direction.currentText()

            # if self.radio_xy.isChecked() is True:
            #     self.parameters["plane"] = "xy"
            # elif self.radio_yz.isChecked() is True:
            #     self.parameters["plane"] = "yz"
            # elif self.radio_zx.isChecked() is True:
            #     self.parameters["plane"] = "zx"
            # if self.check_reverse.isChecked() is True:
            #     self.parameters["plane"] = self.parameters["plane"][::-1]

            self.plot_evaluation()
        return None

    def plot_evaluation(self):
        self.display.EraseAll()
        from generalToolbox.geometry import create_circular_array_2
        from generalToolbox import findInArray

        theta = np.arange(self.parameters['minAngle'],
                          self.parameters["maxAngle"] + self.parameters["step"],
                          self.parameters["step"])

        # xMic = createCircArray(theta, self.parameters['plane'],
        #                        self.parameters['radius'],
        #                        self.parameters['offset'])
        xMic = create_circular_array_2(theta, self.parameters['on_axis'], self.parameters['direction'],
                                       self.parameters['radius'], self.parameters["offset"])
        idx_0, _ = findInArray(theta, 0)

        for i in range(len(theta)):
            x = xMic[i, 0] * 1e3
            y = xMic[i, 1] * 1e3
            z = xMic[i, 2] * 1e3  # system size is in meters so need to scale up
            p = Geom_CartesianPoint(gp_Pnt(x, y, z))
            ais_point = AIS_Point(p)
            if i == idx_0:
                color = Quantity_Color(Quantity_NOC_RED)
            else:
                color = Quantity_Color(Quantity_NOC_BLACK)

            drawer = ais_point.Attributes()
            asp = Prs3d_PointAspect(Aspect_TOM_O, color, 1)
            drawer.SetPointAspect(asp)
            ais_point.SetAttributes(drawer)
            self.display.Context.Display(ais_point, False)

        if self.geo["cad_view_path"]:
            shp = read_step_file(self.geo["cad_view_path"])
            t = TopologyExplorer(shp)
            system = t.my_shape

            self.display.DisplayShape(system, color=Quantity_Color(Quantity_NOC_GRAY40))

        self.display.display_triedron()
        self.display.FitAll()
        self.is_eval_displayed = True
        return None

    def confirm_eval(self):
        # Update eval otherwise it won't register
        if self.is_eval_displayed is False:
            self.update_evaluation()

        name = self.parameters["name"]
        self.evaluation_setup[name] = self.parameters
        self.evaluation_list.addItem(name)

        self.display.EraseAll()
        self.close()


class addPressureFieldWindow(QWidget):
    def __init__(self, system, geo, evaluation_setup, evaluation_list):
        super().__init__()
        self.system = system
        self.geo = geo
        self.evaluation_setup = evaluation_setup
        self.evaluation_list = evaluation_list
        self.setWindowTitle("Pressure Field Evaluation Setup")

        # plotter tab
        self.eval_view = QGraphicsView(self)
        self.eval_view.setGeometry(QRect(210, 20, 461, 371))

        # OpenCascade plotter
        self.plotter = qtDisplay.qtViewer3d(self.eval_view)  # QtInteractor(self.tab) #for pyvistaqt
        self.plotter.resize(461, 371)
        # self.plotter.InitDriver()
        self.display = self.plotter._display
        self.display.View.SetBgGradientColors(
            Quantity_Color(Quantity_NOC_ALICEBLUE),
            Quantity_Color(Quantity_NOC_WHITE),
            2,
            True,
        )
        self.display.Repaint()
        self.is_eval_displayed = False

        # param
        self.parameters = {"name": None,
                           "L1": None,
                           "L2": None,
                           "step": None,
                           "offset": [0, 0, 0],
                           "plane": None,
                           "type": "screen",
                           "toUpdate": False}

        # labels
        labels = ["Evaluation name", "Dimensions (m)", "L1", "L2",
                  "step", "Plane", "Offset", "x", "y", "z"]
        dimensions = [[20, 20, 91, 16], [20, 80, 91, 16], [20, 110, 41, 16],
                      [20, 140, 41, 16], [20, 170, 31, 16],
                      [20, 210, 49, 16], [20, 260, 49, 16], [20, 280, 31, 16],
                      [70, 280, 31, 16], [120, 280, 31, 16]]

        self.label = {}
        for index, label in enumerate(labels):
            self.label[label] = QLabel(label, self)
            x, y, dx, dy = (dimensions[index][0], dimensions[index][1],
                            dimensions[index][2], dimensions[index][3])
            self.label[label].setGeometry(QRect(x, y, dx, dy))

        # lineEdit
        self.lineEdit_name = QLineEdit(self)
        self.lineEdit_name.setGeometry(QRect(20, 40, 131, 21))

        self.lineEdit_L1 = QLineEdit(self)
        self.lineEdit_L1.setGeometry(QRect(50, 110, 41, 21))
        self.lineEdit_L1.setText("1")

        self.lineEdit_L2 = QLineEdit(self)
        self.lineEdit_L2.setGeometry(QRect(50, 140, 41, 21))
        self.lineEdit_L2.setText("1")

        self.lineEdit_step = QLineEdit(self)
        self.lineEdit_step.setGeometry(QRect(50, 170, 41, 21))
        self.lineEdit_step.setText("5e-2")

        self.lineEdit_xOffset = QLineEdit(self)
        self.lineEdit_xOffset.setGeometry(QRect(20, 300, 41, 21))
        self.lineEdit_xOffset.setText("0")

        self.lineEdit_yOffset = QLineEdit(self)
        self.lineEdit_yOffset.setGeometry(QRect(70, 300, 41, 21))
        self.lineEdit_yOffset.setText("0")

        self.lineEdit_zOffset = QLineEdit(self)
        self.lineEdit_zOffset.setGeometry(QRect(120, 300, 41, 21))
        self.lineEdit_zOffset.setText("0")

        # radiobutton
        self.radio_xy = QRadioButton("xy", self)
        self.radio_xy.setGeometry(QRect(20, 230, 40, 20))
        self.radio_xy.toggle()

        self.radio_yz = QRadioButton("yz", self)
        self.radio_yz.setGeometry(QRect(70, 230, 40, 20))

        self.radio_zx = QRadioButton("zx", self)
        self.radio_zx.setGeometry(QRect(120, 230, 40, 20))

        # pushbutton
        self.pushbutton_confirm = QPushButton("Confirm", self)
        self.pushbutton_confirm.setGeometry(QRect(110, 360, 75, 24))
        self.pushbutton_confirm.clicked.connect(self.confirm_eval)

        self.pushbutton_preview = QPushButton("Preview", self)
        self.pushbutton_preview.setGeometry(QRect(20, 360, 75, 24))
        self.pushbutton_preview.clicked.connect(self.update_evaluation)

    def update_evaluation(self):
        if self.lineEdit_name.text():
            self.parameters["name"] = self.lineEdit_name.text()
            self.parameters["L1"] = int(self.lineEdit_L1.text())
            self.parameters["L2"] = int(self.lineEdit_L2.text())
            self.parameters["step"] = float(self.lineEdit_step.text())
            xoff = float(self.lineEdit_xOffset.text())
            yoff = float(self.lineEdit_yOffset.text())
            zoff = float(self.lineEdit_zOffset.text())
            self.parameters["offset"] = [xoff, yoff, zoff]

            if self.radio_xy.isChecked() is True:
                self.parameters["plane"] = "xy"
            elif self.radio_yz.isChecked() is True:
                self.parameters["plane"] = "yz"
            elif self.radio_zx.isChecked() is True:
                self.parameters["plane"] = "zx"

            # print(self.parameters)
            self.plot_evaluation()


    def plot_evaluation(self):
        self.display.EraseAll()
        from generalToolbox.geometry import createPlaneArray
        xMic, Ln, Wn = createPlaneArray(self.parameters["L1"], self.parameters["L2"],
                                self.parameters["step"], self.parameters["plane"],
                                self.parameters["offset"])

        for i in range(len(xMic)):
            x = xMic[i, 0] * 1e3
            y = xMic[i, 1] * 1e3
            z = xMic[i, 2] * 1e3  # system size is in meters so need to scale up
            p = Geom_CartesianPoint(gp_Pnt(x, y, z))
            ais_point = AIS_Point(p)
            color = Quantity_Color(Quantity_NOC_BLACK)

            drawer = ais_point.Attributes()
            asp = Prs3d_PointAspect(Aspect_TOM_O, color, 1)
            drawer.SetPointAspect(asp)
            ais_point.SetAttributes(drawer)
            self.display.Context.Display(ais_point, False)

        if self.geo["cad_view_path"]:
            shp = read_step_file(self.geo["cad_view_path"])
            t = TopologyExplorer(shp)
            system = t.my_shape

            self.display.DisplayShape(system, color=Quantity_Color(Quantity_NOC_GRAY40))

        self.display.FitAll()
        self.is_eval_displayed = True


    def confirm_eval(self):
        # Update eval otherwise it won't register
        if self.is_eval_displayed is False:
            self.update_evaluation()

        name = self.parameters["name"]
        self.evaluation_setup[name] = self.parameters
        self.evaluation_list.addItem(name)

        self.display.EraseAll()
        self.close()



class addPressurePointWindow(QWidget):
    def __init__(self, system, geo, evaluation_setup, evaluation_list):
        super().__init__()
        self.system = system
        self.geo = geo
        self.evaluation_setup = evaluation_setup
        self.evaluation_list = evaluation_list
        self.setWindowTitle("Pressure Point Evaluation Setup")

        # plotter tab
        self.eval_view = QGraphicsView(self)
        self.eval_view.setGeometry(QRect(210, 20, 461, 371))

        # OpenCascade plotter
        self.plotter = qtDisplay.qtViewer3d(self.eval_view)  # QtInteractor(self.tab) #for pyvistaqt
        self.plotter.resize(461, 371)
        # self.plotter.InitDriver()
        self.display = self.plotter._display
        self.display.View.SetBgGradientColors(
            Quantity_Color(Quantity_NOC_ALICEBLUE),
            Quantity_Color(Quantity_NOC_WHITE),
            2,
            True,
        )
        self.display.Repaint()
        self.is_eval_displayed = False

        # param
        self.parameters = {"name": None,
                           "xMic": None,
                           "toUpdate": False,
                           "type": "point"}

        # NAME
        self.label_name = QLabel("Evaluation Name", self)
        self.label_name.setGeometry(QRect(20, 20, 91, 16))

        self.lineEdit_name = QLineEdit(self)
        self.lineEdit_name.setGeometry(QRect(20, 40, 131, 21))

        # table widget
        self.tableWidget = QTableWidget(self)
        self.tableWidget.setGeometry(QRect(10, 80, 191, 231))
        self.tableWidget.setAlternatingRowColors(True)
        self.tableWidget.setColumnCount(3)
        self.tableWidget.setRowCount(3)
        self.tableWidget.setHorizontalHeaderLabels(["x (m)", "y (m)", "z (m)"])
        self.tableWidget.setColumnWidth(0, 58)
        self.tableWidget.setColumnWidth(1, 58)
        self.tableWidget.setColumnWidth(2, 58)

        # pushbuttons
        self.pushbutton_addRow = QPushButton("Add row", self)
        self.pushbutton_addRow.setGeometry(QRect(20, 320, 75, 24))
        self.pushbutton_addRow.clicked.connect(self.add_row_2_table)

        self.pushbutton_confirm = QPushButton("Confirm", self)
        self.pushbutton_confirm.setGeometry(QRect(110, 360, 75, 24))
        self.pushbutton_confirm.clicked.connect(self.confirm_eval)

        self.pushbutton_preview = QPushButton("Preview", self)
        self.pushbutton_preview.setGeometry(QRect(20, 360, 75, 24))
        self.pushbutton_preview.clicked.connect(self.update_evaluation)

    def add_row_2_table(self):
        row_count = self.tableWidget.rowCount()
        self.tableWidget.insertRow(row_count)
        return None

    def update_evaluation(self):
        if self.lineEdit_name.text():
            self.parameters["name"] = self.lineEdit_name.text()

            row_count = self.tableWidget.rowCount()
            col_count = self.tableWidget.columnCount()
            xMic = []
            for row in range(row_count):
                mic_tmp = []
                for col in range(col_count):
                    if self.tableWidget.item(row, col) is not None:
                        mic_tmp.append(float(self.tableWidget.item(row, col).text()))
                if bool(mic_tmp) is True:
                    xMic.append(mic_tmp)
            xMic = np.array(xMic)
            self.parameters["xMic"] = xMic

            # print(self.parameters)
            self.plot_evaluation()


    def plot_evaluation(self):
        self.display.EraseAll()
        xMic = self.parameters["xMic"]

        for i in range(len(xMic)):
            x = xMic[i, 0] * 1e3
            y = xMic[i, 1] * 1e3
            z = xMic[i, 2] * 1e3  # system size is in meters so need to scale up
            p = Geom_CartesianPoint(gp_Pnt(x, y, z))
            ais_point = AIS_Point(p)
            color = Quantity_Color(Quantity_NOC_BLACK)

            drawer = ais_point.Attributes()
            asp = Prs3d_PointAspect(Aspect_TOM_O, color, 1)
            drawer.SetPointAspect(asp)
            ais_point.SetAttributes(drawer)
            self.display.Context.Display(ais_point, False)

        if self.geo["cad_view_path"]:
            shp = read_step_file(self.geo["cad_view_path"])
            t = TopologyExplorer(shp)
            system = t.my_shape

            self.display.DisplayShape(system, color=Quantity_Color(Quantity_NOC_GRAY40))

        self.display.FitAll()
        self.is_eval_displayed = True


    def confirm_eval(self):
        # Update eval otherwise it won't register
        if self.is_eval_displayed is False:
            self.update_evaluation()

        name = self.parameters["name"]
        self.evaluation_setup[name] = self.parameters
        self.evaluation_list.addItem(name)

        self.display.EraseAll()
        self.close()




def fill_polar_window(evaluation, current_eval, polarWindow):
    """
    When user modifies polar window -> read from evaluation and writes in widgets
    :param evaluation:
    :param current_eval:
    :param pressureFieldWindow:
    :return:
    """
    name = evaluation[current_eval]["name"]
    minAngle = str(evaluation[current_eval]["minAngle"])
    maxAngle = str(evaluation[current_eval]["maxAngle"])
    step = str(evaluation[current_eval]["step"])
    radius = str(evaluation[current_eval]["radius"])
    # plane = evaluation[current_eval]["plane"]
    on_axis = evaluation[current_eval]["on_axis"]
    direction = evaluation[current_eval]["direction"]
    axis_list = ["x", "-x", "y", "-y", "z", "-z"]
    offset = evaluation[current_eval]["offset"]
    xOffset = str(offset[0])
    yOffset = str(offset[1])
    zOffset = str(offset[2])

    # fill lineEdit
    polarWindow.lineEdit_name.setText(name)
    polarWindow.lineEdit_minAngle.setText(minAngle)
    polarWindow.lineEdit_maxAngle.setText(maxAngle)
    polarWindow.lineEdit_step.setText(step)
    polarWindow.lineEdit_radius.setText(radius)
    polarWindow.lineEdit_xOffset.setText(xOffset)
    polarWindow.lineEdit_yOffset.setText(yOffset)
    polarWindow.lineEdit_zOffset.setText(zOffset)

    axis_index = axis_list.index(on_axis)
    axis_direction = axis_list.index(direction)

    polarWindow.combo_onAxis.setCurrentIndex(axis_index)
    polarWindow.combo_direction.setCurrentIndex(axis_direction)

    # check plane
    # if plane == "xy":
    #     polarWindow.radio_xy.toggle()
    # elif plane == "yz":
    #     polarWindow.radio_yz.toggle()
    # elif plane == "zx":
    #     polarWindow.radio_zx.toggle()
    # elif plane == "yx":
    #     polarWindow.radio_xy.toggle()
    #     polarWindow.check_reverse.toggle()
    # elif plane == "zy":
    #     polarWindow.radio_yz.toggle()
    #     polarWindow.check_reverse.toggle()
    # elif plane == "xz":
    #     polarWindow.radio_zx.toggle()
    #     polarWindow.check_reverse.toggle()  #.isChecked()

    evaluation[current_eval]["toUpdate"] = True


def fill_pressure_field_window(evaluation, current_eval, pressureFieldWindow):
    """
    When user modifies pressure window -> read from evaluation and writes in widgets
    :param evaluation:
    :param current_eval:
    :param pressureFieldWindow:
    :return:
    """
    name = evaluation[current_eval]["name"]
    L1 = str(evaluation[current_eval]["L1"])
    L2 = str(evaluation[current_eval]["L2"])
    step = str(evaluation[current_eval]["step"])
    plane = evaluation[current_eval]["plane"]
    offset = evaluation[current_eval]["offset"]
    xOffset = str(offset[0])
    yOffset = str(offset[1])
    zOffset = str(offset[2])

    pressureFieldWindow.lineEdit_name.setText(name)
    pressureFieldWindow.lineEdit_L1.setText(L1)
    pressureFieldWindow.lineEdit_L2.setText(L2)
    pressureFieldWindow.lineEdit_step.setText(step)
    pressureFieldWindow.lineEdit_xOffset.setText(xOffset)
    pressureFieldWindow.lineEdit_yOffset.setText(yOffset)
    pressureFieldWindow.lineEdit_zOffset.setText(zOffset)

    if plane == "xy":
        pressureFieldWindow.radio_xy.toggle()
    elif plane == "yz":
        pressureFieldWindow.radio_yz.toggle()
    elif plane == "zx":
        pressureFieldWindow.radio_zx.toggle()

    evaluation[current_eval]["toUpdate"] = True

def fill_pressure_point_window(evaluation, current_eval, pressurePointWindow):
    """
    When user modifies pressure point window -> read from evaluation and writes in widgets
    :param evaluation:
    :param current_eval:
    :param pressureFieldWindow:
    :return:
    """
    name = evaluation[current_eval]["name"]
    xMic = evaluation[current_eval]["xMic"]
    row_count = len(xMic)
    col_count = 3
    pressurePointWindow.tableWidget.setRowCount(row_count)
    pressurePointWindow.label_name.setText(name)

    for row in range(row_count):
        for col in range(col_count):
            pressurePointWindow.tableWidget.setItem(row, col, QTableWidgetItem(str(xMic[row, col])))

    evaluation[current_eval]["toUpdate"] = True





