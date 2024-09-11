from PyQt6.QtCore import QRect, QTimer
from PyQt6.QtWidgets import QWidget, QLabel, QLineEdit, QPushButton, QListWidget, QComboBox, QFrame, QFileDialog, \
    QAbstractItemView, QPlainTextEdit

import numpy as np


class InfoTabEnclosure:
    def __init__(self, info_panel, enclosure_tab, system):
        self.enclosure_tab = enclosure_tab
        self.info_panel = info_panel
        self.system = system
        self.tab = QWidget()

        self.label_enclosureSelect = QLabel("Select enclosure", self.tab)
        self.label_enclosureSelect.setGeometry(QRect(10, 20, 91, 16))

        # combobox
        self.comboSelection = QComboBox(self.tab)
        self.comboSelection.setGeometry(QRect(110, 20, 91, 22))
        self.comboSelection.currentIndexChanged.connect(self.updateBoxInfoText)
        self.nEnclosure = 0

        # Text display
        self.textEnclosure_display = QPlainTextEdit(self.tab)
        self.textEnclosure_display.setGeometry(QRect(10, 60, 261, 161))
        self.textEnclosure_display.setReadOnly(True)

        # update info tab
        self.timer_enclosure_detect = QTimer(self.tab)
        self.timer_enclosure_detect.timeout.connect(self.buildInfoTab)
        self.timer_enclosure_detect.start(1000)  # Check every second

        info_panel.addTab(self.tab, "Enclosure")

    def buildInfoTab(self):
        nEnclosure_tmp = self.enclosure_tab.listWidget.count()

        if nEnclosure_tmp != self.nEnclosure:
            self.comboSelection.clear()
            for index in range(nEnclosure_tmp):
                item = self.enclosure_tab.listWidget.item(index)
                self.comboSelection.addItem(item.text())
            self.nEnclosure = np.copy(nEnclosure_tmp)

    def updateBoxInfoText(self):
        self.textEnclosure_display.clear()
        box = self.comboSelection.currentText()
        boxType = self.system.enclosure[box].config

        if boxType == "sealed":
            Vb = self.system.enclosure[box].Vb
            ref2bem = self.system.enclosure[box].ref2bem
            if self.system.enclosure[box].whichDriver is not None:
                driver = self.system.enclosure[box].whichDriver
            else:
                driver = "None"
            text = ("--- Enclosure Parameters ---" + "\n" +
                    "Config: sealed" + "\n" +
                    "Vb = " + str(round(Vb, 3)*1e3) + " L" + "\n" +
                    "--- BEM Setup ---" + "\n" +
                    "BEM reference: " + str(ref2bem) + "\n" +
                    "Driver: " + driver)
            self.textEnclosure_display.setPlainText(text)
        elif boxType == "vented":
            Vb = self.system.enclosure[box].Vb
            Lp = self.system.enclosure[box].Lp
            Sp = self.system.enclosure[box].Sp
            rp = self.system.enclosure[box].rp

            ref2bem = self.system.enclosure[box].ref2bem
            if self.system.enclosure[box].whichDriver is not None:
                driver = self.system.enclosure[box].whichDriver
            else:
                driver = "None"

            text =  ("--- Enclosure Parameters ---" + "\n" +
                    "Config: vented" + "\n" +
                    "Vb = " + str(round(Vb, 3)*1e3) + " L" + "\n" +
                    "Lp = " + str(round(Lp, 3)*1e2) + " cm" + "\n" +
                    "rp = " + str(round(rp, 3) * 1e2) + " cm" + "\n" +
                    "Sp = " + str(round(Sp, 3)*1e4) + " cm^2" + "\n" +
                    "--- BEM Setup ---" + "\n" +
                    "BEM reference: " + str(ref2bem) + "\n" +
                    "Driver: " + driver)
            self.textEnclosure_display.setPlainText(text)
