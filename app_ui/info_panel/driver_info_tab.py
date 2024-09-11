from PyQt6.QtCore import QRect, QTimer
from PyQt6.QtWidgets import QWidget, QLabel, QLineEdit, QPushButton, QListWidget, QComboBox, QFrame, QFileDialog, \
    QAbstractItemView, QPlainTextEdit

import numpy as np

class InfoTabDriver:
    def __init__(self, info_panel, driver_tab, system):
        self.driver_tab = driver_tab
        self.info_panel = info_panel
        self.system = system

        self.tab = QWidget()
        # self.buildInfoTab(control_panel.driver_tab)

        # labels
        self.label_driverSelect = QLabel("Select driver", self.tab)
        self.label_driverSelect.setGeometry(QRect(10, 20, 71, 16))

        # combobox
        self.comboSelection = QComboBox(self.tab)
        self.comboSelection.setGeometry(QRect(80, 20, 91, 22))
        self.comboSelection.currentIndexChanged.connect(self.updateLPMInfoText)
        self.nDriver = 0

        # Text display
        self.textLPM_display = QPlainTextEdit(self.tab)
        self.textLPM_display.setGeometry(QRect(10, 60, 261, 161))
        self.textLPM_display.setReadOnly(True)

        # update info tab
        self.timer_driver_detect = QTimer(self.tab)
        self.timer_driver_detect.timeout.connect(self.buildInfoTab)
        self.timer_driver_detect.start(1000)  # Check every second


        info_panel.addTab(self.tab, "Driver")


    def buildInfoTab(self):
        nDriver_tmp = self.driver_tab.listWidget.count()

        if nDriver_tmp != self.nDriver:
            self.comboSelection.clear()
            for index in range(nDriver_tmp):
                item = self.driver_tab.listWidget.item(index)
                self.comboSelection.addItem(item.text())
            self.nDriver = np.copy(nDriver_tmp)

    def updateLPMInfoText(self):
        self.textLPM_display.clear()

        # load drv object
        drv = self.comboSelection.currentText()

        Electrical = "--- Electrical ---"
        Re = "Re = " + str(round(self.system.driver[drv].Re, 3)) + " Ohm"
        Le = "Le = " + str(round(self.system.driver[drv].Le*1e3, 3)) + " mH"
        Bl = "Bl = " + str(round(self.system.driver[drv].Bl, 3)) + " T.m"

        Mechanical = "--- Mechanical ---"
        Mms = "Mms = " + str(round(self.system.driver[drv].Mms*1e3, 3)) + " g"
        Cms = "Cms = " + str(round(self.system.driver[drv].Cms*1e3, 3)) + " mm/N"
        Kms = "Kms = " + str(round(1 / self.system.driver[drv].Cms, 3)) + " N/m"
        Rms = "Rms = " + str(round(self.system.driver[drv].Rms, 3)) + " N.s/m"
        Sd = "Sd = " + str(round(self.system.driver[drv].Sd*1e4, 3)) + " cm^2"

        QualityFactors = "--- Quality Factors ---"
        Qes = "Qes = " + str(round(self.system.driver[drv].Qes, 3))
        Qms = "Qms = " + str(round(self.system.driver[drv].Qms, 3))
        Qts = "Qts = " + str(round(self.system.driver[drv].Qts, 3))

        Other = "--- Other ---"
        Fs = "Fs = " + str(round(self.system.driver[drv].Fs, 3)) + "  Hz"
        Vas = "Vas = " + str(round(self.system.driver[drv].Vas*1e3, 3)) + " L"
        if self.system.driver[drv].ref2bem is False:
            Ref2bem = "BEM reference: None"
        else:
            Ref2bem = "BEM reference: " + str(self.system.driver[drv].ref2bem)

        self.textLPM_display.setPlainText(
            Electrical + "\n" +
            Re + "\n" +
            Le + "\n" +
            Bl + "\n" +
            Mechanical + "\n" +
            Mms + "\n" +
            Cms + "\n" +
            Kms + "\n" +
            Rms + "\n" +
            Sd + "\n" +
            QualityFactors + "\n" +
            Qes + "\n" +
            Qms + "\n" +
            Qts + "\n" +
            Other + "\n" +
            Fs + "\n" +
            Vas + "\n" +
            Ref2bem
            )
