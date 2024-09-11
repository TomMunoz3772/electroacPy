# from OCC.STEPControl import STEPControl_Reader
# import os
# import sys
#
# from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
# from OCC.Extend.TopologyUtils import TopologyExplorer
# from OCC.Extend.DataExchange import read_step_file
# from OCC.Display.OCCViewer import OffscreenRenderer
# from OCC.Display.backend import load_backend
#
# load_backend("pyqt6")
# import OCC.Display.qtDisplay as qtDisplay
#
# from PyQt6.QtWidgets import QDialog, QVBoxLayout, QGroupBox, QHBoxLayout, QPushButton, QApplication
#
# shp = read_step_file("606s3_straight_edges.step")
# t = TopologyExplorer(shp)
#
# # my_renderer = OffscreenRenderer()
# # my_renderer.DisplayShape(t.my_shape)
#
#
# class App(QDialog):
#     def __init__(self):
#         super().__init__()
#         self.title = "PyQt6 / pythonOCC"
#         self.left = 300
#         self.top = 300
#         self.width = 800
#         self.height = 300
#         self.initUI()
#
#     def initUI(self):
#         self.setWindowTitle(self.title)
#         self.setGeometry(self.left, self.top, self.width, self.height)
#         self.createHorizontalLayout()
#
#         windowLayout = QVBoxLayout()
#         windowLayout.addWidget(self.horizontalGroupBox)
#         self.setLayout(windowLayout)
#         self.show()
#         self.canvas.InitDriver()
#         self.canvas.resize(200, 200)
#         self.display = self.canvas._display
#
#     def createHorizontalLayout(self):
#         self.horizontalGroupBox = QGroupBox("Display PythonOCC")
#         layout = QHBoxLayout()
#
#         disp = QPushButton("Display Box", self)
#         disp.clicked.connect(self.displayBOX)
#         layout.addWidget(disp)
#
#         eras = QPushButton("Erase Box", self)
#         eras.clicked.connect(self.eraseBOX)
#         layout.addWidget(eras)
#
#         self.canvas = qtDisplay.qtViewer3d(self)
#         layout.addWidget(self.canvas)
#         self.horizontalGroupBox.setLayout(layout)
#
#     def displayBOX(self):
#         shp = read_step_file("606s3_straight_edges.step")
#         t = TopologyExplorer(shp)
#         a_box = t.my_shape  #BRepPrimAPI_MakeBox(10.0, 20.0, 30.0).Shape()
#         self.ais_box = self.display.DisplayShape(a_box)[0]
#         self.display.FitAll()
#
#     def eraseBOX(self):
#         self.display.Context.Erase(self.ais_box, True)
#
#
# if __name__ == "__main__":
#     app = QApplication(sys.argv)
#     window = App()
#     window.show()
#     sys.exit(app.exec())

import pyvista as pv
import numpy as np

# Load the .msh file
mesh = pv.read("mesh_tmp.msh")

# Display information about the mesh
print(mesh)

# Assuming the enclosure surfaces are marked with a specific scalar value
# For this example, let's assume the "enclosure" surfaces have a scalar value of 1
enclosure_value = 1

# Create a boolean array where True indicates the enclosure cells
# enclosure_cells = mesh.cell_data['ScalarValues'] == enclosure_value
max_scalar = np.max(mesh.active_scalars)
# Create a color array where enclosure cells are white and others are red
colors = []

for scalar in mesh.active_scalars:
    if scalar != max_scalar:
        colors.append([255, 177, 0])  # White color for enclosure
    else:
        colors.append([87, 212, 82])  # Red color for non-enclosure

# Add the colors to the mesh
mesh.cell_data['colors'] = colors

# Create a plotter
plotter = pv.Plotter()

# Add the mesh to the plotter with the specified colors
plotter.add_mesh(mesh, scalars='colors', rgb=True, show_edges=True)

# Show the plot
plotter.show()
