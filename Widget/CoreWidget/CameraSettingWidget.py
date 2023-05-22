from PyQt5.QtWidgets import *
# from PIL import Image
# import numpy as np
# from pathlib import Path

from PyQt5.QtCore import *
from PyQt5 import QtGui
from PyQt5.QtGui import QFont
from Utilities.IO.IOHelper import get_camconfig_setting, CAMCONFIGT_FILE_PATH
from Model.Instruments.Camera.Chameleon import Chameleon
from TUCAM import *
import Utilities.Helper.settings as settings


class CameraSetting(QWidget):
    img_sub = pyqtSignal(object)
    img_sub2 = pyqtSignal(object)

    def __init__(self, parent=None):
        super(CameraSetting, self).__init__(parent=parent)
        self.parent = parent

        self.GroupBox1 = QGroupBox()
        layouth = QHBoxLayout()
        layoutT = QGridLayout()#camera setting

        self.AbsTriger = QCheckBox("AbsImg", self)
        self.video_mode = QCheckBox("Video", self)
        self.hardware_mode = QCheckBox("ExtTrig", self)
        self.cb = QComboBox()#camera select
        self.further_setting = QPushButton("CamSetting")
        self.auto_save = QCheckBox("AutoSave", self)
        self.prefix_label = QLabel('ImgFilePrefix', self)
        # self.prefix_label.setAlignment(Qt.AlignRight)
        self.prefix_text = QLineEdit('Default', self)
        # self.lab = QLabel('',self)
        layoutT.addWidget(self.AbsTriger,1,2,1,1)
        layoutT.addWidget(self.video_mode,1,0,1,1)
        layoutT.addWidget(self.hardware_mode,1,1,1,1)
        layoutT.addWidget(self.cb,0,0,1,3)
        layoutT.addWidget(self.further_setting,0,3,1,1)
        layoutT.addWidget(self.auto_save,1,3,1,1)
        layoutT.addWidget(self.prefix_label, 2, 0, 1, 2)
        layoutT.addWidget(self.prefix_text, 2, 2, 1, 2)
        # layoutT.addWidget(self.lab,2,0,1,1)
        self.GroupBox1.setLayout(layoutT)
        layouth.addWidget(self.GroupBox1)

        self.further_setting.clicked.connect(self.camera_setting)
        # after click the start experiment button, then can change camera setting in detail
        self.further_setting.setEnabled(True)

        layoutv = QVBoxLayout()
        layoutv.addLayout(layouth)

        self.GroupBox1.setFont(QFont("Roman times", 12))

        self.setLayout(layoutv)

        self.default_setting()

        screen = QtGui.QDesktopWidget().screenGeometry()#Control window size
        self.setFixedWidth(screen.width()*23/100)

        self.d = QDialog()  # create a dialog
        dialog_layout = QVBoxLayout()
        self.apply_button = QPushButton("Apply")
        self.camera_further_setting = CameraSettingWidget()  # set three parameters
        dialog_layout.addWidget(self.camera_further_setting)
        dialog_layout.addWidget(self.apply_button)
        self.d.setLayout(dialog_layout)
        self.detect()
        self.cb.activated.connect(self.select_camera_index)

    #######################
    def detect(self):
        existingdata = get_camconfig_setting()
        #################
        # Path = './'
        # # Path = 'AtomImgSCNU/Debug/tiffile/'
        # TUCAM_Api_Uninit = TUSDKdll.TUCAM_Api_Uninit
        # TUCAM_Api_Uninit
        
        # TUCAM_Api_Init = TUSDKdll.TUCAM_Api_Init
        # TUCAMINIT = TUCAM_INIT(0, Path.encode('utf-8'))
        # # TUCAM_Api_Init(pointer(TUCAMINIT));
       
        # print("TUCAMINIT.uiCamCount is ",TUCAMINIT.uiCamCount)
        # print("Chameleon.getPortInfo() is ",Chameleon.getPortInfo())
        # print("TUCAMINIT.pstrConfigPath is ",TUCAMINIT.pstrConfigPath)
        try:
            Path = './'
            TUCAM_Api_Init = TUSDKdll.TUCAM_Api_Init
            TUCAMINIT = TUCAM_INIT(0, Path.encode('utf-8'))
            TUCAM_Api_Init(pointer(TUCAMINIT));
            print(TUCAMINIT.uiCamCount)
            print(TUCAMINIT.pstrConfigPath)
            TUCAM_Dev_Open = TUSDKdll.TUCAM_Dev_Open
            TUCAMOPEN = TUCAM_OPEN(0, 0)
            TUCAM_Dev_Open(pointer(TUCAMOPEN));
            print(TUCAMOPEN.uiIdxOpen)
            print(TUCAMOPEN.hIdxTUCam)

            # Get Camera Info
            TUCAM_Dev_GetInfo = TUSDKdll.TUCAM_Dev_GetInfo
            # Camera name:  
            m_infoid = TUCAM_IDINFO
            TUCAMVALUEINFO = TUCAM_VALUE_INFO(m_infoid.TUIDI_CAMERA_MODEL.value, 0, 0, 0)
            TUCAM_Dev_GetInfo(c_int64(TUCAMOPEN.hIdxTUCam), pointer(TUCAMVALUEINFO))
            print(TUCAMVALUEINFO.pText)

            # Camera VID
            TUCAMVALUEINFO = TUCAM_VALUE_INFO(m_infoid.TUIDI_VENDOR.value, 0, 0, 0)
            TUCAM_Dev_GetInfo(c_int64(TUCAMOPEN.hIdxTUCam), pointer(TUCAMVALUEINFO))
            print('%#X'%TUCAMVALUEINFO.nValue)

            # Camera PID
            TUCAMVALUEINFO = TUCAM_VALUE_INFO(m_infoid.TUIDI_PRODUCT.value, 0, 0, 0)
            TUCAM_Dev_GetInfo(c_int64(TUCAMOPEN.hIdxTUCam), pointer(TUCAMVALUEINFO))
            print('%#X'%TUCAMVALUEINFO.nValue)

            # Sdk API
            TUCAMVALUEINFO = TUCAM_VALUE_INFO(m_infoid.TUIDI_VERSION_API.value, 0, 0, 0)
            TUCAM_Dev_GetInfo(c_int64(TUCAMOPEN.hIdxTUCam), pointer(TUCAMVALUEINFO))
            print(TUCAMVALUEINFO.pText)


            # CloseCamera
            TUCAM_Dev_Close = TUSDKdll.TUCAM_Dev_Close
            TUCAM_Dev_Close(c_int64(TUCAMOPEN.hIdxTUCam))
            # TUCAM_Api_Uninit = TUSDKdll.TUCAM_Api_Uninit
            # TUCAM_Api_Uninit
            print("camera closed")

            camera_infos = '%#X'%TUCAMVALUEINFO.nValue
            # camera_infos = Chameleon.getPortInfo()
            camera_infos2 = []
            m = 0
            if camera_infos is not None:
                for n in camera_infos:
                    n = n[0:-2]
                    try:
                        camindex = existingdata.index(n)
                        camera_infos2.append(existingdata[camindex+1]+'&'+str(m))
                    except ValueError:
                        #Add cameras that have not been added
                        camera_infos2.append(camera_infos[m])
                        config_path = CAMCONFIGT_FILE_PATH
                        caminfo = open(config_path, 'r+')
                        caminfo.read()
                        caminfo.write("{},".format(camera_infos[m][0:-2]))
                        caminfo.write("{},".format('newcamera'))
                        caminfo.write("{}\n".format(10))
                        caminfo.close()
                    m = m + 1
            if camera_infos is not None:
                self.cb.clear()
                self.cb.addItems(camera_infos2)  # get camera num
                # if detect cameras, default camera index is 0
                # settings.instrument_params["Camera"]["index"] = 0
                self.select_camera_index()
                self.cb.setEnabled(True)

            else:
                print("No camera detected !")
                self.cb.clear()
                self.cb.setEnabled(False)
                settings.instrument_params["Camera"]["index"] = None
                # return
        except:
            print("No camera(TUCAM) detected")

    def select_camera_index(self):
        settings.instrument_params["Camera"]["index"] = self.cb.currentText()[-1]
        existingdata = get_camconfig_setting()
        dataindex = existingdata.index(self.cb.currentText()[0:-2])
        try:
            settings.instrument_params["Camera"]["shutter time"] = existingdata[dataindex + 1]
            self.camera_further_setting.shutter_time.setValue(float(settings.instrument_params["Camera"]["shutter time"]))
        except:
            settings.instrument_params["Camera"]["shutter time"] = existingdata[dataindex + 2]
            self.camera_further_setting.shutter_time.setValue(float(settings.instrument_params["Camera"]["shutter time"]))



    def camera_setting(self):
        self.d.setWindowTitle(self.cb.currentText())
        self.d.setWindowModality(Qt.ApplicationModal)
        self.d.exec_()

    def default_setting(self):
        self.video_mode.setChecked(True)
        # self.hardware_mode.setChecked(False)
        self.video_mode.setEnabled(False)


#点击CamSetting按钮弹出的相机参数窗口
class CameraSettingWidget(QWidget):
    """
        camera setting and control widget for initialization and running,
        including basis camera settings and control.

    """
    def __init__(self, parent=None):
        super(CameraSettingWidget, self).__init__(parent)
        self.parent = parent
        self.GroupBox = QGroupBox("Camera Setting")
        layout = QVBoxLayout()
        exposure = QHBoxLayout()
        self.exposure_time_label = QLabel("Exposure time: ")
        self.exposure_time = QDoubleSpinBox()
        self.exposure_time.setRange(10, 80)
        self.exposure_time.setSingleStep(1)
        exposure.addWidget(self.exposure_time_label)
        exposure.addWidget(self.exposure_time)

        shutter = QHBoxLayout()
        self.shutter_label = QLabel("Shutter time: ")
        self.shutter_time = QDoubleSpinBox()
        self.shutter_time.setRange(10, 80)
        self.shutter_time.setSingleStep(1)
        shutter.addWidget(self.shutter_label)
        shutter.addWidget(self.shutter_time)

        gain = QHBoxLayout()
        self.gain_label = QLabel("Gain: ")
        self.gain_value = QDoubleSpinBox()
        self.gain_value.setRange(1, 10)
        self.gain_value.setSingleStep(1)
        gain.addWidget(self.gain_label)
        gain.addWidget(self.gain_value)
       
       #layout.addLayout(shutter)
        layout.addLayout(exposure)
        
        # layout.addLayout(gain)

        self.GroupBox.setLayout(layout)
        self.vertical_layout = QVBoxLayout()
        self.vertical_layout.addWidget(self.GroupBox)

        self.setLayout(self.vertical_layout)

        self.default_setting()

        # self.exposure_time.valueChanged.connect(self.change_exposure)
        # self.shutter_time.valueChanged.connect(self.change_shutter)
        # self.gain_value.valueChanged.connect(self.change_gain)

    def default_setting(self):
        self.shutter_time.setValue(settings.instrument_params["Camera"]["shutter time"])
        self.exposure_time.setValue(settings.instrument_params["Camera"]["exposure time"])
        self.gain_value.setValue(settings.instrument_params["Camera"]["gain value"])

    def change_shutter(self):
        settings.instrument_params["Camera"]["shutter time"] = self.shutter_time.value()
        # print("shutter time is ", settings.instrument_params["Camera"]["shutter time"])

    def change_exposure(self):
        settings.instrument_params["Camera"]["exposure time"] = self.exposure_time.value()
        # print("exposure time is ", settings.instrument_params["Camera"]["exposure time"])

    def change_gain(self):
        settings.instrument_params["Camera"]["gain value"] = self.gain_value.value()
        # print("gain value is ", settings.instrument_params["Camera"]["gain value"])

