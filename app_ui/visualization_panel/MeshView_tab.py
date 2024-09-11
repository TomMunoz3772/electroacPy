from PyQt6.QtCore import QRect, QTimer
from PyQt6.QtWidgets import QWidget, QLabel, QLineEdit, QPushButton, QListWidget, QComboBox, QFrame, QVBoxLayout
from pyvistaqt import QtInteractor
import pyvista as pv
import numpy as np

class MeshView:
    def __init__(self, control_panel, system, geo):
        self.system = system
        self.geo = geo
        self.control_panel = control_panel

        self.tab = QWidget()

        # Create PyVista plotter
        self.plotter = QtInteractor(self.tab)
        self.plotter.setMinimumSize(810, 600)
        self.plotter.interactor.setMinimumSize(810, 600)

        # prepare mesh
        mesh = pv.read(self.geo["meshPath"])
        max_scalar = np.max(mesh.active_scalars)
        colors = []
        for scalar in mesh.active_scalars:
            if scalar != max_scalar:
                colors.append([255, 177, 0])  # Green color for enclosure
            else:
                colors.append([87, 212, 82])  # Red color for non-enclosure
        mesh.cell_data['colors'] = colors

        # plot mesh
        self.mesh_actor = "MESH_DISPLAY_ACTOR"
        self.plotter.add_mesh(mesh, show_edges=True, show_scalar_bar=False,
                              name=self.mesh_actor, scalars='colors', rgb=True)
        self.plotter.add_axes()
        self.plotter.resize(810, 600)

        control_panel.addTab(self.tab, "Mesh View")

    def updateMesh(self):
        self.plotter.remove_actor(self.mesh_actor)

        # reload mesh
        mesh = pv.read(self.geo["meshPath"])
        max_scalar = np.max(mesh.active_scalars)
        colors = []
        for scalar in mesh.active_scalars:
            if scalar != max_scalar:
                colors.append([255, 177, 0])  # Green color for enclosure
            else:
                colors.append([87, 212, 82])  # Red color for non-enclosure
        mesh.cell_data['colors'] = colors

        # plot
        self.plotter.add_mesh(mesh, show_edges=True, show_scalar_bar=False,
                              name=self.mesh_actor, scalars='colors', rgb=True)
        return None
