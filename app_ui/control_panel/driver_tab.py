from PyQt6.QtCore import QRect
from PyQt6.QtWidgets import QWidget, QLabel, QLineEdit, QPushButton, QListWidget, QFileDialog
from electroacPy.speakerSim.electroAcousticDriver import loadLPM

class DriverTab:
    def __init__(self, control_panel, system):
        self.system = system
        self.tab = QWidget()
        self.driver_lpm_window = None
        self.driver_vib_window = None

        # add from lpm data
        self.pushButton_lpm = QPushButton(self.tab)
        self.pushButton_lpm.setGeometry(QRect(10, 40, 91, 24))
        self.pushButton_lpm.setText("Add from LPM")
        self.pushButton_lpm.clicked.connect(self.addDriverFromLPM)

        # add from acceleration data
        self.pushButton_vib = QPushButton(self.tab)
        self.pushButton_vib.setGeometry(QRect(10, 80, 141, 24))
        self.pushButton_vib.setText("Add acceleration data")
        self.pushButton_vib.clicked.connect(self.addDriverFromAcceleration)

        # list with currently defined drivers
        self.listWidget = QListWidget(self.tab)
        self.listWidget.setGeometry(QRect(10, 120, 256, 192))

        # push button to plot selected driver info
        self.pushButton_plot = QPushButton(self.tab)
        self.pushButton_plot.setGeometry(QRect(180, 80, 75, 24))
        self.pushButton_plot.setText("Plot data")
        self.pushButton_plot.clicked.connect(self.plotDriverData)

        control_panel.addTab(self.tab, "Driver")

    def addDriverFromLPM(self):
        if self.driver_lpm_window is None:
            self.driver_lpm_window = DriverLpmWindow(self.system, self.listWidget)
        self.driver_lpm_window.show()

    def addDriverFromAcceleration(self):
        print("Acceleration - not implemented yet")

    def plotDriverData(self):
        selected_driver = self.listWidget.currentItem().text()
        if selected_driver:
            self.system.driver[selected_driver].plotXVA()
            self.system.driver[selected_driver].plotZe()
        else:
            pass


class DriverLpmWindow(QWidget):
    def __init__(self, loudspeakerSystem, driverList):
        super().__init__()
        self.system = loudspeakerSystem
        self.driverList = driverList
        self.setWindowTitle("Linear Parameters")
        labels = ["Re (Ohm)", "Le (mH)", "Mms (g)", "Cms (µm/N)",
                  "Rms ()", "Bl (T.m)", "Sd (cm^2)"]
        # self.ll = labels
        py = [50, 80, 120, 150, 180, 220, 250]
        self.label = {}
        self.lineEdit = {}
        for i in range(len(labels)):
            self.label[labels[i]] = QLabel(self)
            self.label[labels[i]].setText(labels[i])
            self.label[labels[i]].setGeometry(QRect(20, py[i], 66, 20))
            self.lineEdit[labels[i]] = QLineEdit(self)
            self.lineEdit[labels[i]].setGeometry(QRect(110, py[i], 40, 20))

        self.pushButton_load = QPushButton(self)
        self.pushButton_load.setGeometry(QRect(10, 290, 81, 24))
        self.pushButton_load.setText("Load LPM file")
        self.pushButton_load.clicked.connect(self.loadLPM_param)

        self.pushButton_confirm = QPushButton(self)
        self.pushButton_confirm.setGeometry(QRect(100, 290, 75, 24))
        self.pushButton_confirm.setText("Add driver")
        self.pushButton_confirm.clicked.connect(self.addDriverPushed)

        # Window name label
        self.label_windowName = QLabel(self)
        self.label_windowName.setGeometry(QRect(20, 20, 141, 16))
        self.label_windowName.setText("Linear parameters")

        # 2nd column
        self.label_driverName = QLabel(self)
        self.label_driverName.setGeometry(QRect(170, 20, 100, 16))
        self.label_driverName.setText("Driver name")
        self.lineEdit["driver name"] = QLineEdit(self)
        self.lineEdit["driver name"].setGeometry(QRect(170, 50, 151, 21))

        # input voltage and BEM reference
        self.label_inputVoltage = QLabel(self)
        self.label_inputVoltage.setGeometry(QRect(170, 80, 71, 16))
        self.label_inputVoltage.setText("Input voltage")
        self.lineEdit["input voltage"] = QLineEdit(self)
        self.lineEdit["input voltage"].setGeometry(QRect(270, 80, 51, 21))

        self.label_bemref = QLabel(self)
        self.label_bemref.setGeometry(QRect(170, 120, 81, 16))
        self.label_bemref.setText("BEM reference")
        self.lineEdit["bem ref"] = QLineEdit(self)
        self.lineEdit["bem ref"].setGeometry(QRect(270, 120, 51, 21))

    def addDriverPushed(self):
        # print(self.lineEdit["Re (Ohm)"].text())
        Name = self.lineEdit["driver name"].text()
        Re = float(self.lineEdit["Re (Ohm)"].text())
        Le = float(self.lineEdit["Le (mH)"].text()) * 1e-3
        Mms = float(self.lineEdit["Mms (g)"].text()) * 1e-3
        Cms = float(self.lineEdit["Cms (µm/N)"].text()) * 1e-6
        Rms = float(self.lineEdit["Rms ()"].text())
        Bl = float(self.lineEdit["Bl (T.m)"].text())
        Sd = float(self.lineEdit["Sd (cm^2)"].text()) * 1e-4
        U = float(self.lineEdit["input voltage"].text())
        try:  # check if user inputs a reference for BEM
            bem = float(self.lineEdit["bem ref"].text())
        except:
            bem = None
        if self.lineEdit["driver name"].text():
            self.system.lem_driver(Name, U, Le, Re, Cms, Mms, Rms, Bl, Sd, ref2bem=bem)
            self.driverList.addItem(Name)
        else:
            pass

    def loadLPM_param(self):
        """
        Open a file dialog to load LPM data from a Klippel .txt file
        """
        # options = QFileDialog()
        file_name, _ = QFileDialog.getOpenFileName(self, "Search File", "",
                                                   "All Files (*);;Text Files (*.txt)")
        ead = loadLPM(file_name, self.system.frequency)
        self.lineEdit["Re (Ohm)"].setText(str(ead.Re))
        self.lineEdit["Le (mH)"].setText(str(ead.Le*1e3))
        self.lineEdit["Mms (g)"].setText(str(ead.Mms*1e3))
        self.lineEdit["Cms (µm/N)"].setText(str(ead.Cms*1e6))
        self.lineEdit["Rms ()"].setText(str(ead.Rms))
        self.lineEdit["Bl (T.m)"].setText(str(ead.Bl))
        self.lineEdit["Sd (cm^2)"].setText(str(ead.Sd*1e4))
        self.lineEdit["input voltage"].setText(str(1))



