from PyQt6.QtCore import QRect
from PyQt6.QtWidgets import QWidget, QLabel, QLineEdit, QPushButton, QListWidget, QComboBox, QFrame, QFileDialog, \
    QAbstractItemView
from generalToolbox import meshCAD

class CADTab:
    def __init__(self, control_panel, system, geo):
        self.system = system
        self.tab = QWidget()
        self.geo = geo
        self.group_selection_window = None
        self.file_name = None

        # LOAD STEP FILE
        self.label_stp = QLabel(self.tab)
        self.label_stp.setGeometry(QRect(10, 10, 81, 16))
        self.label_stp.setText("*.step location")

        self.pushbutton_search_step = QPushButton(self.tab)
        self.pushbutton_search_step.setGeometry(QRect(200, 30, 75, 24))
        self.pushbutton_search_step.setText("Search")
        self.pushbutton_search_step.clicked.connect(self.search_step_file)

        self.lineEdit_location = QLineEdit(self.tab)
        self.lineEdit_location.setGeometry(QRect(10, 30, 181, 21))

        self.pushbutton_load_step = QPushButton(self.tab)
        self.pushbutton_load_step.setGeometry(QRect(10, 60, 91, 24))
        self.pushbutton_load_step.setText("Load")
        self.pushbutton_load_step.clicked.connect(self.load_step_file)

        # self.line = QFrame(self.tab)
        # self.line.setGeometry(QRect(10, 90, 261, 16))
        # self.line.setFrameShape()

        # MESH SIZE DEFINITION
        self.label_maxFreq = QLabel(self.tab)
        self.label_maxFreq.setGeometry(QRect(10, 120, 111, 16))
        self.label_maxFreq.setText("Max frequency (Hz)")

        self.lineEdit_maxFreq = QLineEdit(self.tab)
        self.lineEdit_maxFreq.setGeometry(QRect(150, 120, 91, 21))
        self.lineEdit_maxFreq.setText(str(int(self.system.frequency[-1])))

        self.label_ppw = QLabel(self.tab)
        self.label_ppw.setGeometry(QRect(10, 150, 121, 16))
        self.label_ppw.setText("Points per wavelength")

        self.lineEdit_ppw = QLineEdit(self.tab)
        self.lineEdit_ppw.setGeometry(QRect(150, 150, 31, 21))
        self.lineEdit_ppw.setText("6")

        self.label_minFactor = QLabel(self.tab)
        self.label_minFactor.setGeometry(QRect(10, 180, 81, 16))
        self.label_minFactor.setText("Min size factor")

        self.lineEdit_minFactor = QLineEdit(self.tab)
        self.lineEdit_minFactor.setGeometry(QRect(150, 180, 31, 21))
        self.lineEdit_minFactor.setText("10")

        # GROUP SURFACES
        self.pushbutton_def_rad = QPushButton(self.tab)
        self.pushbutton_def_rad.setGeometry(QRect(10, 210, 141, 24))
        self.pushbutton_def_rad.setText("Define radiating surface")
        self.pushbutton_def_rad.clicked.connect(self.group_surfaces)
        self.listWidget_radSurf = QListWidget(self.tab)
        self.listWidget_radSurf.setGeometry(QRect(10, 240, 181, 71))

        # MESH
        self.pushbutton_mesh = QPushButton(self.tab)
        self.pushbutton_mesh.setGeometry(QRect(200, 290, 75, 24))
        self.pushbutton_mesh.setText("Mesh")
        self.pushbutton_mesh.clicked.connect(self.mesh_geo)

        control_panel.addTab(self.tab, "CAD")

    def search_step_file(self):
        file_name, _ = QFileDialog.getOpenFileName(self.tab, "Search File", "", "All Files (*);;Text Files (*.txt)")
        self.lineEdit_location.setText(file_name)

    def load_step_file(self):
        if self.lineEdit_location.text():
            self.geo["cad_view_path"] = self.lineEdit_location.text()
            self.geo["cad_changed"] = True
        return None

    def group_surfaces(self):
        # check mesh size
        maxFreq = float(self.lineEdit_maxFreq.text())
        ppw = int(self.lineEdit_ppw.text())
        minSizeFactor = float(self.lineEdit_minFactor.text())
        self.geo["maxSize"] = self.system.c / maxFreq / ppw
        self.geo['minSize'] = self.geo["maxSize"] / minSizeFactor
        self.geo['ppw'] = ppw

        # group surfaces
        if self.group_selection_window is None:
            # print(self.listWidget_radSurf.item())
            self.group_selection_window = GroupSelectionWindow(self.system, self.geo, self.listWidget_radSurf)
        else:
            if self.geo["meshCAD"]:  # otherwise it doesn't properly "clean" the 3D viewer
                import gmsh
                gmsh.initialize()
                gmsh.finalize()
            self.group_selection_window.close()
            self.group_selection_window = GroupSelectionWindow(self.system, self.geo, self.listWidget_radSurf)
        self.group_selection_window.show()

    def mesh_geo(self):
        """
        Mesh CAD geometry using surface group defined using GroupSelectionWindow
        :return:
        """

        # recheck mesh size incase it changes
        maxFreq = float(self.lineEdit_maxFreq.text())
        ppw = int(self.lineEdit_ppw.text())
        minSizeFactor = float(self.lineEdit_minFactor.text())
        self.geo["maxSize"] = self.system.c / maxFreq / ppw
        self.geo['minSize'] = self.geo["maxSize"] / minSizeFactor

        # Create GEO and mesh
        self.geo["meshCAD"] = None
        self.geo["meshCAD"] = meshCAD(self.geo["cad_view_path"],
                                      minSize=self.geo['minSize'],
                                      maxSize=self.geo['maxSize'])

        for name in self.geo['surfaceGroup']:
            surfaceGroup_tmp = self.geo['surfaceGroup'][name]
            surfaceList_tmp = self.geo['surfaceList'][name]
            surfaceSize_tmp = self.geo['surfaceSize'][name]
            self.geo["meshCAD"].addSurfaceGroup(name, surfaceList_tmp, surfaceGroup_tmp, surfaceSize_tmp)

        self.geo["meshCAD"].mesh("mesh_tmp")

        import os
        current_directory = os.getcwd()
        self.geo["meshPath"] = os.path.join(current_directory, "mesh_tmp.msh")
        self.geo["update_mesh_display"] = True
        return None


class GroupSelectionWindow(QWidget):
    def __init__(self, loudspeakerSystem, geometry, listRadSurf):
        super().__init__()
        self.system = loudspeakerSystem
        self.geo = geometry
        self.list_radSurf = listRadSurf
        self.setWindowTitle("Radiating surface selection")
        self.group_attribution_window = None

        # labels
        self.label_CAD = QLabel(self)
        self.label_CAD.setGeometry(QRect(80, 20, 31, 16))
        self.label_CAD.setText("CAD")

        self.label_radGroup = QLabel(self)
        self.label_radGroup.setGeometry(QRect(270, 20, 101, 16))
        self.label_radGroup.setText("Radiating groups")

        # push buttons
        self.pushbutton_add = QPushButton(self)
        self.pushbutton_add.setGeometry(QRect(20, 50, 75, 24))
        self.pushbutton_add.setText("Add")
        self.pushbutton_add.clicked.connect(self.addSurfaceToGroup)

        self.pushbutton_remove = QPushButton(self)
        self.pushbutton_remove.setGeometry(QRect(220, 50, 75, 24))
        self.pushbutton_remove.setText("Remove")
        self.pushbutton_remove.clicked.connect(self.removeGroup)

        self.pushbutton_validate = QPushButton(self)
        self.pushbutton_validate.setGeometry(QRect(310, 280, 75, 24))
        self.pushbutton_validate.setText("Validate")
        self.pushbutton_validate.clicked.connect(self.validate)

        # Qlistwidgets
        self.list_faces = QListWidget(self)
        self.list_faces.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.list_faces.setGeometry(QRect(20, 85, 171, 192))

        self.list_group = QListWidget(self)
        self.list_group.setGeometry(QRect(220, 85, 171, 192))


        # Fill surface list view
        if self.geo['cad_view_path']:
            self.list_faces.clear()
            self.geo["meshCAD"] = None
            self.geo["meshCAD"] = meshCAD(self.geo["cad_view_path"])  # get surface list
            for i in range(len(self.geo["meshCAD"].surface_list)):
                self.list_faces.addItem(str(self.geo["meshCAD"].surface_list[i]))

        # Fill groups list view
        if self.geo['surfaceList'] is not False:
            for name in self.geo['surfaceList']:
                self.list_group.addItem(name)

    def addSurfaceToGroup(self):
        selectedSurfaces = []
        for i in range(len(self.list_faces.selectedItems())):
            selectedSurfaces.append(int(self.list_faces.selectedItems()[i].text()))
        self.geo["tmp_faces"] = selectedSurfaces

        if self.group_attribution_window is None:
            self.group_attribution_window = GroupAttributionWindow(self.geo,
                                                                   group_listWidget=self.list_group,
                                                                   group_listWidget_CADTAB=self.list_radSurf)
        self.group_attribution_window.show()
        return None

    def removeGroup(self):
        # might need to modify that! TODO -> check in study_tap class to see how to remove items
        if self.list_group.currentItem():
            for index in range(self.list_group.count()):
                item = self.list_group.item(index)
                if item.text() == self.list_group.currentItem().text():
                    # print("ON EST ICI:: ", item.text())
                    self.geo['surfaceList'].pop(item.text())
                    self.geo['surfaceSize'].pop(item.text())
                    self.geo['surfaceGroup'].pop(item.text())
                    self.list_group.takeItem(index)
                    self.list_radSurf.takeItem(index)
        return None

    def validate(self):
        self.close()
        return None


class GroupAttributionWindow(QWidget):
    def __init__(self, geometry, group_listWidget,
                 group_listWidget_CADTAB):
        super().__init__()
        self.geo = geometry
        self.list_group = group_listWidget
        self.list_group_CAD_TAB = group_listWidget_CADTAB
        self.setWindowTitle("Set to group")

        # labels
        self.label_set2group = QLabel(self)
        self.label_set2group.setGeometry(QRect(10, 10, 81, 16))
        self.label_set2group.setText("Group number")

        self.label_meshSize = QLabel(self)
        self.label_meshSize.setGeometry(QRect(100, 10, 71, 21))
        self.label_meshSize.setText("Mesh size")

        self.label_mm = QLabel(self)
        self.label_mm.setGeometry(QRect(180, 40, 49, 16))
        self.label_mm.setText("[mm]")

        self.label_groupName = QLabel(self)
        self.label_groupName.setGeometry(QRect(10, 70, 81, 16))
        self.label_groupName.setText("Group name")

        # QlineEdit
        self.lineEdit_group = QLineEdit(self)
        self.lineEdit_group.setGeometry(QRect(10, 30, 71, 21))

        self.lineEdit_meshSize = QLineEdit(self)
        self.lineEdit_meshSize.setGeometry((QRect(100, 30, 71, 21)))
        self.lineEdit_meshSize.setText(str(round(self.geo['maxSize'] * 1e3, 2)))

        self.lineEdit_groupName = QLineEdit(self)
        self.lineEdit_groupName.setGeometry(QRect(100, 70, 111, 21))

        # pushbutton
        self.pushbutton_validate = QPushButton(self)
        self.pushbutton_validate.setGeometry(QRect(10, 100, 75, 24))
        self.pushbutton_validate.setText("Validate")
        self.pushbutton_validate.clicked.connect(self.createRadiatingGroup)

    def createRadiatingGroup(self):
        name = self.lineEdit_groupName.text()
        # faces = self.geo['tmp_faces']
        meshSize = float(self.lineEdit_meshSize.text())
        group_id = int(self.lineEdit_group.text())
        # self.geo['meshCAD'].addSurfaceGroup(name, faces, group_id, meshSize)

        # group creation in GEO file -> mesh is done later
        list_widget_id = "{} - ".format(group_id) + name
        self.geo['surfaceList'][list_widget_id] = self.geo["tmp_faces"]
        self.geo['surfaceGroup'][list_widget_id] = group_id
        self.geo['surfaceSize'][list_widget_id] = meshSize

        # add group name to list widget
        item = []
        for i in range(self.list_group.count()):
            # print(self.list_group.item(i).text())
            item.append(self.list_group.item(i).text())

        if list_widget_id not in item:
            self.list_group.addItem(list_widget_id)
            self.list_group_CAD_TAB.addItem(list_widget_id)
        else:
            for index in range(self.list_group.count()):
                item = self.list_group.item(index)
                if item.text() == list_widget_id:
                    self.list_group.takeItem(index)
                    self.list_group.addItem(list_widget_id)
                    self.list_group.takeItem(list_widget_id)
                    self.list_group_CAD_TAB.addItem(list_widget_id)

        self.close()
        return None
