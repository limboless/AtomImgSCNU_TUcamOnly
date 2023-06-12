from PyQt5.QtWidgets import *
import numpy as np
from pathlib import Path
from PyQt5.QtCore import *
from PyQt5 import QtGui
from PyQt5.QtGui import QFont, QIntValidator, QDoubleValidator
import Utilities.Helper.settings as settings


class ImgDisplaySetting(QWidget):
    img_sub = pyqtSignal(object)
    img_sub2 = pyqtSignal(object)

    def __init__(self, parent=None):
        super(ImgDisplaySetting, self).__init__(parent=parent)
        self.parent = parent

        # background image
        self.GroupBox2 = QGroupBox()
        layoutT2 = QGridLayout()#camera setting

        self.bkgStatus = QCheckBox("Deduct RefImg", self)
        self.bkgLoad = QPushButton("Load RefImg", self)

        # photon filter
        self.pfStatus = QPushButton('Image Filter', self)

        self.pfMin = QLineEdit(self)
        self.pfMin.setValidator(QIntValidator(self.pfMin))
        self.pfMin.setMaxLength(5)

        self.tolab = QLabel('to', self)
        self.pfMax = QLineEdit(self)
        self.pfMax.setValidator(QIntValidator(self.pfMax))
        self.pfMax.setMaxLength(5)

        self.roi = QCheckBox("ROI", self)
        self.roi_size = QLineEdit(str(settings.widget_params["Analyse Data Setting"]["roisize"]), self)
        self.roi_size.setValidator(QIntValidator(self.roi_size))
        self.roi_size.setMaxLength(3)
        self.rawdata = QCheckBox("CrossHair", self)
        # self.cross_axes = QCheckBox("cross axes", self)
        self.roicnt = QLabel(" ROI-cnt", self)
        # self.roicnt.setAlignment(Qt.AlignRight)
        self.xpos = QLineEdit(str(settings.widget_params["Analyse Data Setting"]["xpos"]), self)
        self.ypos = QLineEdit(str(settings.widget_params["Analyse Data Setting"]["ypos"]), self)
        self.xpos.setValidator(QDoubleValidator(self.xpos))
        self.xpos.setMaxLength(6)
        self.xpos.setToolTip('X axis')
        self.ypos.setValidator(QDoubleValidator(self.ypos))
        self.ypos.setMaxLength(6)
        self.ypos.setToolTip('Y axis')
        self.realtime = QCheckBox("realtime", self)
        self.colormap = QComboBox()
        self.colormap.addItems(settings.colorlist)
        # self.colormap.activated.connect()



        layoutT2.addWidget(self.bkgStatus,0,0,1,2)
        layoutT2.addWidget(self.bkgLoad,0,2,1,3)
        layoutT2.addWidget(self.pfStatus,1,0,1,2)
        layoutT2.addWidget(self.pfMin,1,2,1,1)
        layoutT2.addWidget(self.tolab,1,3,1,1)
        layoutT2.addWidget(self.pfMax,1,4,1,1)
        layoutT2.addWidget(self.roi,2,0,1,1)
        layoutT2.addWidget(self.roi_size,2,1,1,1)
        layoutT2.addWidget(self.rawdata,2,2,1,2)
        layoutT2.addWidget(self.realtime,2,4,1,1)
        layoutT2.addWidget(self.roicnt,3,0,1,1)
        layoutT2.addWidget(self.xpos,3,1,1,1)
        layoutT2.addWidget(self.ypos,3,2,1,1)
        layoutT2.addWidget(self.colormap,3,4,1,1)

        self.GroupBox2.setLayout(layoutT2)
        self.GroupBox2.setFont(QFont("Roman times", 12))

        layoutv = QVBoxLayout()
        layoutv.addWidget(self.GroupBox2)

        self.setLayout(layoutv)
        self.default_setting()

        self.bkgStatus.stateChanged.connect(lambda: self.ckbstate(self.bkgStatus))
        self.pfStatus.clicked.connect(lambda: self.ckbstate(self.pfStatus))

        self.pfMin.textChanged.connect(self.changePfMin)
        self.pfMax.textChanged.connect(self.changePfMax)

        self.bkgLoad.clicked.connect(self.loadbkgImg)


        screen = QtGui.QDesktopWidget().screenGeometry()#Control window size
        self.setFixedWidth(screen.width()*24/100)

        ########
    def change_roiposition(self,position):
        if position[0] != self.xpos.text():
            settings.widget_params["Analyse Data Setting"]["xpos"] = position[0]
            self.xpos.setText(str(settings.widget_params["Analyse Data Setting"]["xpos"]))
        if position[1] != self.ypos.text():
            settings.widget_params["Analyse Data Setting"]["ypos"] = position[1]
            self.ypos.setText(str(settings.widget_params["Analyse Data Setting"]["ypos"]))

    def loadbkgImg(self):
        path = QFileDialog.getOpenFileName(self, "Open File")  # name path
        strimg_path = str(path)
        img_file = strimg_path[2:len(strimg_path) - 19]
        img_path = Path(img_file)

        pathjud = str(img_path)
        pathjud = pathjud[len(pathjud) - 3:]  # Get the version of the file
        if pathjud == 'ata':
            file = open(img_path)
            linesc = file.readlines()  # Read the file as a behavior unit
            rows = len(linesc)  # get the numbers fo line
            lines = len(linesc[0].strip().split(' '))
            img_data = np.zeros((rows, lines))  # Initialization matrix
            row = 0
            for line2 in linesc:
                line2 = line2.strip().split(' ')
                img_data[row, :] = line2[:]
                row += 1
            file.close()
        else:
            try:
                img_data = np.load(img_path, encoding='bytes')
            except TypeError:
                return
            except PermissionError:
                return

        img = img_data[::-1]
        # path, _ = QFileDialog.getOpenFileName(self, 'Open Image', 'c:\\', 'Image files(*.jpg *.gif *.png)')
        # img = Image.open(path)
        settings.imgData["BkgImg"] = np.array(img)
        print('The background image has been added.')

    def default_setting(self):
        # self.bkgStatus.setChecked(settings.widget_params["Image Display Setting"]["bkgStatus"])
        # self.pfStatus.setChecked(settings.widget_params["Image Display Setting"]["pfStatus"])

        self.pfMax.setText(str(settings.widget_params["Image Display Setting"]["pfMax"]))
        self.pfMin.setText(str(settings.widget_params["Image Display Setting"]["pfMin"]))

    def changeMagValue(self):
        settings.widget_params["Image Display Setting"]["magValue"] = self.magValue.value()

    def changePfMin(self):
        settings.widget_params["Image Display Setting"]["pfMin"] = self.pfMin.text()

    def changePfMax(self):
        settings.widget_params["Image Display Setting"]["pfMax"] = self.pfMax.text()

    def ckbstate(self, b):
        if b.text() == "Deduct RefImg":
            if b.isChecked():
                # if b.isChecked() == True:
                if settings.imgData["BkgImg"] !=[]:
                    settings.widget_params["Image Display Setting"]["bkgStatus"] = True
                    print('subtract background image ： finish.')
                    self.img_sub.emit(1)
                    # print("background status", settings.widget_params["Image Display Setting"]["bkgStatus"])
                else:
                    print('No background image.')
                # else:
                #     settings.widget_params["Image Display Setting"]["bkgStatus"] = False

        if b.text() == "Image Filter":
            # if b.isChecked() == True:
            settings.widget_params["Image Display Setting"]["pfStatus"] = True
            self.img_sub2.emit(1)
                # print('photon filter ： finish.')
                # print("photon filter Status status", settings.widget_params["Image Display Setting"]["pfStatus"])
            # else:
            #     settings.widget_params["Image Display Setting"]["pfStatus"] = False


