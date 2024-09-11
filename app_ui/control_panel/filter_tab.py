from PyQt6.QtCore import QRect
from PyQt6.QtWidgets import (QWidget, QLabel, QLineEdit, QPushButton, QListWidget, QComboBox,
                             QCheckBox, QRadioButton, QGraphicsView, QTableWidget, QTableWidgetItem)


import numpy as np


class FilterTab:
    def __init__(self, control_panel, system, evaluation_setup, driver_tab, enclosure_tab):
        # need to get which surface radiates from the driver and enclosure tab (similar to what is done in study_tab)
        self.control_panel = control_panel
        self.system = system
        self.evaluation_setup = evaluation_setup
        self.driver_tab = driver_tab
        self.enclosure_tab = enclosure_tab
        


