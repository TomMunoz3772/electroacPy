from PyQt6.QtCore import QRect, QTimer
from PyQt6.QtWidgets import QWidget, QLabel, QLineEdit, QPushButton, QListWidget, QComboBox, QFrame, QVBoxLayout
from pyvistaqt import QtInteractor
import pyvista as pv

from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Extend.DataExchange import read_step_file, read_step_file_with_names_colors
from OCC.Display.backend import load_backend
from OCC.Core.Quantity import (Quantity_Color, Quantity_TOC_RGB, Quantity_NOC_ALICEBLUE,
                               Quantity_NOC_ANTIQUEWHITE, Quantity_NOC_IVORY,
                               Quantity_NOC_WHITE, Quantity_NOC_GRAY40)

load_backend("pyqt6")
import OCC.Display.qtDisplay as qtDisplay


class CADView:
    def __init__(self, control_panel, system, geo):
        self.system = system
        self.geo = geo
        self.tab = QWidget()
        self.geometryRender = None

        # Create a layout for the tab widget
        # layout = QVBoxLayout()
        # self.tab.setLayout(layout)

        # Create the OpenCascade plotter
        self.plotter = qtDisplay.qtViewer3d(self.tab)  #QtInteractor(self.tab) #for pyvistaqt
        self.plotter.resize(810, 600)
        self.plotter.InitDriver()
        self.display = self.plotter._display
        self.display.View.SetBgGradientColors(
            Quantity_Color(Quantity_NOC_ALICEBLUE),
            Quantity_Color(Quantity_NOC_WHITE),
            2,
            True,
        )
        self.display.display_triedron()
        self.display.Repaint()

        # Create a QTimer to check for changes every second
        self.timer = QTimer(self.tab)
        self.timer.timeout.connect(self.check_for_changes)
        self.timer.start(1000)  # Check every second

        # Create a label to display surface information
        self.info_label = QLabel(self.tab)
        self.info_label.setGeometry(QRect(10, 580, 111, 16))
        self.info_label.setText("Selection ID: N/A")
        control_panel.addTab(self.tab, "CAD view")

    def check_for_changes(self):
        has_changed = self.geo["cad_changed"]
        if has_changed is True:
            self.update_cad_view()
            self.geo["cad_changed"] = False

    def update_cad_view(self):
        shp = read_step_file(self.geo["cad_view_path"])
        t = TopologyExplorer(shp)
        a_box = t.my_shape  # BRepPrimAPI_MakeBox(10.0, 20.0, 30.0).Shape()

        if self.geometryRender is None:
            self.geometryRender = self.display.DisplayShape(a_box,
                                                            color=Quantity_Color(Quantity_NOC_GRAY40))
        else:
            self.display.EraseAll()
            self.geometryRender = self.display.DisplayShape(a_box,
                                                            color=Quantity_Color(Quantity_NOC_GRAY40))
        self.display.FitAll()
        # self.display.SetSelectionModeFace()
        return None

    def getFaceID(self):

        return
