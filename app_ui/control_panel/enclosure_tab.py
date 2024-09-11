from PyQt6.QtCore import QRect
from PyQt6.QtWidgets import QWidget, QLabel, QLineEdit, QPushButton, QListWidget, QComboBox
import ast


class EnclosureTab:
    def __init__(self, control_panel, loudspeakerSystem, driverTab):
        self.system = loudspeakerSystem
        self.control_panel = control_panel
        self.tab = QWidget()

        # setup window initialization
        self.sealed_enclosure_window = None
        self.ported_enclosure_window = None
        self.pr_enclosure_window = None

        # tab label
        self.label_tab = QLabel(self.tab)
        self.label_tab.setGeometry(QRect(10, 20, 91, 16))
        self.label_tab.setText("Enclosure type")

        # combo box for enclosure type selection
        self.comboBox_enclosure = QComboBox(self.tab)
        self.comboBox_enclosure.setGeometry(QRect(10, 40, 191, 22))
        self.comboBox_enclosure.addItem("Sealed")
        self.comboBox_enclosure.addItem("Ported")
        self.comboBox_enclosure.addItem("Passive radiator")
        self.comboBox_enclosure.addItem("Bandpass (4th order ported)")
        self.comboBox_enclosure.addItem("Bandpass (4th order PR)")

        # pushbutton to add enclosure
        self.pushbutton_add_enclosure = QPushButton(self.tab)
        self.pushbutton_add_enclosure.setGeometry(QRect(10, 80, 91, 24))
        self.pushbutton_add_enclosure.setText("Add enclosure")
        self.pushbutton_add_enclosure.clicked.connect(lambda: self.enclosure_setup(driverTab))
        control_panel.addTab(self.tab, "Enclosure")

        # pushbutton to plot acceleration and impedance
        self.pushbutton_plot_data = QPushButton(self.tab)
        self.pushbutton_plot_data.setGeometry(QRect(180, 80, 75, 24))
        self.pushbutton_plot_data.setText("Plot Data")
        self.pushbutton_plot_data.clicked.connect(self.plot_enclosure_data)

        # listWidget of all enclosure defined
        self.listWidget = QListWidget(self.tab)
        self.listWidget.setGeometry(QRect(10, 120, 256, 192))

    def enclosure_setup(self, driverTab):
        if self.comboBox_enclosure.currentText() == "Sealed":
            if self.sealed_enclosure_window is None:
                self.sealed_enclosure_window = SealedEnclosureWindow(self.system, self.listWidget, driverTab.listWidget)
            else:
                self.sealed_enclosure_window = SealedEnclosureWindow(self.system, self.listWidget, driverTab.listWidget)
            self.sealed_enclosure_window.show()
        elif self.comboBox_enclosure.currentText() == "Ported":
            if self.ported_enclosure_window is None:
                self.ported_enclosure_window = PortedEnclosureWindow(self.system, self.listWidget, driverTab.listWidget)
            self.ported_enclosure_window.show()

    def plot_enclosure_data(self):
        selected_box = self.listWidget.currentItem().text()
        if selected_box:
            self.system.enclosure[selected_box].plotXVA()
            self.system.enclosure[selected_box].plotZe()
            print(self.system.enclosure[selected_box].config)
        else:
            pass


class SealedEnclosureWindow(QWidget):
    def __init__(self, loudspeakerSystem, enclosureListWidget, driverList):
        super().__init__()
        self.system = loudspeakerSystem
        self.driverList = driverList
        self.enclosureListWidget = enclosureListWidget
        self.setWindowTitle("Sealed Enclosure Setup")

        # Volume parameter
        self.label_Vb = QLabel(self)
        self.label_Vb.setGeometry(QRect(20, 40, 49, 16))
        self.label_Vb.setText("Vb (L)")
        self.lineEdit_Vb = QLineEdit(self)
        self.lineEdit_Vb.setGeometry(QRect(100, 40, 51, 21))

        # Loss factor
        self.label_eta = QLabel(self)
        self.label_eta.setGeometry(QRect(20, 70, 61, 16))
        self.label_eta.setText("Loss factor")
        self.lineEdit_eta = QLineEdit(self)
        self.lineEdit_eta.setGeometry(QRect(100, 70, 51, 21))
        self.lineEdit_eta.setText(str(1e-5))

        # wiring and speaker number
        self.label_wiring = QLabel(self)
        self.label_wiring.setGeometry(QRect(290, 170, 49, 16))
        self.label_wiring.setText("Wiring")
        self.comboBox_w = QComboBox(self)
        self.comboBox_w.setGeometry(QRect(290, 190, 68, 22))
        self.comboBox_w.addItem("parallel")
        self.comboBox_w.addItem("series")

        self.label_Nd = QLabel(self)
        self.label_Nd.setGeometry(QRect(290, 240, 81, 16))
        self.label_Nd.setText("Driver number")
        self.lineEdit_Nd = QLineEdit(self)
        self.lineEdit_Nd.setGeometry(QRect(290, 260, 71, 21))
        self.lineEdit_Nd.setText("1")

        # pushbutton to get tuning help
        self.pushbutton_tuning = QPushButton(self)
        self.pushbutton_tuning.setGeometry(QRect(290, 130, 81, 24))
        self.pushbutton_tuning.setText("Tuning helper")
        self.pushbutton_tuning.clicked.connect(self.tuning_helper)

        # pushbutton to add enclosure to system
        self.pushbutton_add_box = QPushButton(self)
        self.pushbutton_add_box.setGeometry(QRect(30, 390, 91, 24))
        self.pushbutton_add_box.setText("Add enclosure")
        self.pushbutton_add_box.clicked.connect(self.add_sealed_enclosure)

        # BEM reference
        self.label_bem = QLabel(self)
        self.label_bem.setGeometry(QRect(160, 340, 81, 16))
        self.label_bem.setText("BEM reference")
        self.lineEdit_bem = QLineEdit(self)
        self.lineEdit_bem.setGeometry(QRect(160, 360, 41, 21))

        # Enclosure Name
        self.label_box_name = QLabel(self)
        self.label_box_name.setGeometry(QRect(20, 340, 101, 16))
        self.label_box_name.setText("Enclosure name")
        self.lineEdit_box_name = QLineEdit(self)
        self.lineEdit_box_name.setGeometry(QRect(20, 360, 113, 21))

        # ListWidget -> add the driver list from the driver tab
        self.label_driver_select = QLabel(self)
        self.label_driver_select.setGeometry(QRect(20, 110, 121, 16))
        self.label_driver_select.setText("Driver selection")
        self.listWidget = QListWidget(self)
        self.listWidget.setGeometry(QRect(20, 130, 256, 192))
        for index in range(driverList.count()):
            item = driverList.item(index)
            self.listWidget.addItem(item.text())

    def tuning_helper(self):
        selected_driver = self.listWidget.currentItem().text()
        if selected_driver:
            align = self.system.driver[selected_driver].sealedAlignment()
        else:
            align = None
        return align

    def add_sealed_enclosure(self):
        selected_driver = self.listWidget.currentItem().text()
        bem = self.lineEdit_bem.text()
        Vb = float(self.lineEdit_Vb.text()) * 1e-3
        eta = float(self.lineEdit_eta.text())
        name = self.lineEdit_box_name.text()
        wiring = self.comboBox_w.currentText()
        Nd = int(self.lineEdit_Nd.text())

        try:
            bem = string_to_list(bem)
        except:
            bem = string_to_list("[" + bem + "]")

        if selected_driver:
            if Vb:
                if eta:
                    if name:
                        if Nd:
                            self.system.lem_enclosure_2(name, Vb,
                                                        eta=eta,
                                                        setDriver=selected_driver,
                                                        wiring=wiring,
                                                        Nd=Nd, ref2bem=bem)
                            self.enclosureListWidget.addItem(name)

        return None


class PortedEnclosureWindow(QWidget):
    def __init__(self, loudspeakerSystem, enclosureListWidget, driverList):
        super().__init__()
        self.system = loudspeakerSystem
        self.driverList = driverList
        self.enclosureListWidget = enclosureListWidget
        self.setWindowTitle("Ported Enclosure Setup")

        # Volume parameter
        self.label_Vb = QLabel(self)
        self.label_Vb.setGeometry(QRect(20, 40, 49, 16))
        self.label_Vb.setText("Vb (L)")
        self.lineEdit_Vb = QLineEdit(self)
        self.lineEdit_Vb.setGeometry(QRect(100, 40, 51, 21))

        # Loss factor
        self.label_eta = QLabel(self)
        self.label_eta.setGeometry(QRect(20, 70, 61, 16))
        self.label_eta.setText("Loss factor")
        self.lineEdit_eta = QLineEdit(self)
        self.lineEdit_eta.setGeometry(QRect(100, 70, 51, 21))
        self.lineEdit_eta.setText(str(1e-5))

        # Port length
        self.label_Lp = QLabel(self)
        self.label_Lp.setGeometry(QRect(210, 40, 49, 16))
        self.label_Lp.setText("Lp (cm)")
        self.lineEdit_Lp = QLineEdit(self)
        self.lineEdit_Lp.setGeometry(QRect(260, 40, 51, 21))

        # Port Radius
        self.label_rp = QLabel(self)
        self.label_rp.setGeometry(QRect(210, 70, 49, 16))
        self.label_rp.setText("rp (cm)")
        self.lineEdit_rp = QLineEdit(self)
        self.lineEdit_rp.setGeometry(QRect(260, 70, 51, 21))

        # wiring and speaker number
        self.label_wiring = QLabel(self)
        self.label_wiring.setGeometry(QRect(290, 170, 49, 16))
        self.label_wiring.setText("Wiring")
        self.comboBox_w = QComboBox(self)
        self.comboBox_w.setGeometry(QRect(290, 190, 68, 22))
        self.comboBox_w.addItem("parallel")
        self.comboBox_w.addItem("series")

        self.label_Nd = QLabel(self)
        self.label_Nd.setGeometry(QRect(290, 240, 81, 16))
        self.label_Nd.setText("Driver number")
        self.lineEdit_Nd = QLineEdit(self)
        self.lineEdit_Nd.setGeometry(QRect(290, 260, 71, 21))
        self.lineEdit_Nd.setText("1")

        # pushbutton to get tuning help
        self.pushbutton_tuning = QPushButton(self)
        self.pushbutton_tuning.setGeometry(QRect(290, 130, 81, 24))
        self.pushbutton_tuning.setText("Tuning helper")
        self.pushbutton_tuning.clicked.connect(self.tuning_helper)

        # pushbutton to add enclosure to system
        self.pushbutton_add_box = QPushButton(self)
        self.pushbutton_add_box.setGeometry(QRect(30, 390, 91, 24))
        self.pushbutton_add_box.setText("Add enclosure")
        self.pushbutton_add_box.clicked.connect(self.add_sealed_enclosure)

        # BEM reference
        self.label_bem = QLabel(self)
        self.label_bem.setGeometry(QRect(160, 340, 81, 16))
        self.label_bem.setText("BEM reference")

        self.sublabel_driver = QLabel(self)
        self.sublabel_port = QLabel(self)
        self.sublabel_driver.setGeometry(QRect(170, 380, 49, 16))
        self.sublabel_port.setGeometry(QRect(220, 380, 49, 16))
        self.sublabel_driver.setText("driver")
        self.sublabel_port.setText("port")

        self.lineEdit_bem = QLineEdit(self)
        self.lineEdit_bem.setGeometry(QRect(170, 360, 41, 21))
        self.lineEdit_bemp = QLineEdit(self)
        self.lineEdit_bemp.setGeometry(QRect(220, 360, 41, 21))

        # Enclosure Name
        self.label_box_name = QLabel(self)
        self.label_box_name.setGeometry(QRect(20, 340, 101, 16))
        self.label_box_name.setText("Enclosure name")
        self.lineEdit_box_name = QLineEdit(self)
        self.lineEdit_box_name.setGeometry(QRect(20, 360, 113, 21))

        # ListWidget -> add the driver list from the driver tab
        self.label_driver_select = QLabel(self)
        self.label_driver_select.setGeometry(QRect(20, 110, 121, 16))
        self.label_driver_select.setText("Driver selection")
        self.listWidget = QListWidget(self)
        self.listWidget.setGeometry(QRect(20, 130, 256, 192))
        for index in range(driverList.count()):
            item = driverList.item(index)
            self.listWidget.addItem(item.text())

    def tuning_helper(self):
        selected_driver = self.listWidget.currentItem().text()
        if selected_driver:
            align = self.system.driver[selected_driver].portedAlignment()
        else:
            align = None
        return align

    def add_sealed_enclosure(self):
        selected_driver = self.listWidget.currentItem().text()
        bem_drv = self.lineEdit_bem.text()
        bem_port = self.lineEdit_bemp.text()
        Vb = float(self.lineEdit_Vb.text()) * 1e-3
        eta = float(self.lineEdit_eta.text())
        name = self.lineEdit_box_name.text()
        wiring = self.comboBox_w.currentText()
        Nd = int(self.lineEdit_Nd.text())
        Lp = float(self.lineEdit_Lp.text()) * 1e-2
        rp = float(self.lineEdit_rp.text()) * 1e-2

        try:
            bem_drv = string_to_list(bem_drv)
        except:
            bem_drv = string_to_list("[" + bem_drv + "]")
        bem_port = string_to_list("[" + bem_port + "]")

        # parameters = ['selected_driver',
        #               'Vb', 'rp', 'Lp', 'eta',
        #               'name', 'Nd', 'bem_drv',
        #               'bem_port']

        if None in (selected_driver, bem_drv, bem_port, Vb, eta, name, wiring, Nd, Lp, rp):
            pass
        else:
            self.system.lem_enclosure_2(name, Vb, Lp=Lp, rp=rp,
                                        eta=eta,
                                        setDriver=selected_driver,
                                        wiring=wiring,
                                        Nd=Nd, ref2bem=bem_drv + bem_port)
            self.enclosureListWidget.addItem(name)
        return None


def string_to_list(input_str):
    try:
        # Use ast.literal_eval to safely evaluate the string
        result = ast.literal_eval(input_str)
        # Ensure the result is a list of integers
        if isinstance(result, list) and all(isinstance(i, int) for i in result):
            return result
        else:
            raise ValueError("Input string does not represent a list of integers")
    except (ValueError, SyntaxError):
        raise ValueError("Invalid input string")
