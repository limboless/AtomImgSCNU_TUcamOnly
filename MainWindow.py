#from types import ClassMethodDescriptorType
#############from Model.Instruments.Camera.Chameleon import Chameleon
###from TUDefine import *
# from Widget.CoreWidget.AnalyseDataWidget import PlotWindow
from Widget.CoreWidget.ImgQueueWidget import PlotWindow
from TUCAM import *
from Utilities.Helper import settings, Helper
from Utilities.IO import IOHelper
from Model.DataAnalysis.CaculateAtoms import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtGui import QIntValidator

from Widget.CoreWidget.PlotMainWindowWidget import PlotMainWindow
from Widget.CoreWidget.ImgQueueWidget import ImgQueueWidget
from Widget.CoreWidget.ImgDisplaySettingWidget import ImgDisplaySetting
from Widget.CoreWidget.ParametersWidget import ImgParameters
# from Widget.CoreWidget.PromptWidget import PromptWidget
from Widget.CoreWidget.ResultWidget import Absorption, Flurence
from Widget.CoreWidget.PromptWidget import PromptWidget
from Widget.CoreWidget.FittingdataWidget import FittingdataWidget
from Widget.CoreWidget.CameraSettingWidget import CameraSetting

import numpy as np
from numpy import save
import sys
from PIL import Image
import time
from pathlib import Path
import datetime
from Utilities.IO.IOHelper import get_camconfig_setting
from queue import Queue
import warnings

from pathlib import Path
#from AutoLoadImg import AutoLoadImg

class MainWindow(QMainWindow):

    sig_abort_workers = pyqtSignal()
    tempnum = pyqtSignal(object)

    def __init__(self):
        super(MainWindow, self).__init__()
        self.move(82, 0)
        ### MENUS AND TOOLBARS ###
        self.fileMenu = self.menuBar().addMenu("File")
        self.windowMenu = self.menuBar().addMenu("Window")
        self.cameraMenu = self.menuBar().addMenu("Camera")
        self.optionMenu = self.menuBar().addMenu("Options")
        fpath = IOHelper.get_config_setting('DATA_PATH')
        self.path = self.menuBar().addMenu('##  Save file to: ' + str(fpath) + '  ##')
        # self.path.setEnabled(True)

        self.plotToolbar = self.addToolBar("Plot")
        self.expToolbar = self.addToolBar("Experiment")

        # experiment start/stop buttons
        self.start_exp_action = Helper.create_action(self, "Start Experiment", slot=self.start_exp, icon="start")#name and action and image
        self.stop_exp_action = Helper.create_action(self, "Stop Experiment", slot=self.stop_exp, icon="stop")
        self.detect_action = Helper.create_action(self, "Detect Camera", slot=self.detect, icon="cam")#name and action and image
        self.capture_action = Helper.create_action(self, "capture image", slot=self.capture_image, icon="capture")
        self.stop_exp_action.setEnabled(False)

        # plot buttons
        self.clear_img_stack_action = Helper.create_action(self, "clear image stack", slot=self.clear_img_stack, icon="clear_img_stack")
        self.clear_main_win_action = Helper.create_action(self, "clear main window", slot=self.clear_main_win, icon="clear_main_win")


        ### CREATE WIDGET ###
        # global parameters全局参数
        settings.inintParams()

        self.plot_main_window = PlotMainWindow()
        self.setCentralWidget(self.plot_main_window)# set central

        # Create a slider 创建滑动条
        self.img_queue = ImgQueueWidget()
        slider1 = QWidget()
        self.scroll = QScrollArea()     #为堆栈添加滑动条
        self.scroll.setWidget(self.img_queue)
        self.slidervbox = QVBoxLayout()
        self.slidervbox.addWidget(self.scroll)
        slider1.setLayout(self.slidervbox)
        screen = QDesktopWidget().screenGeometry()
        slider1.setFixedSize(screen.width()*13.5/100,screen.width()*(9/16)*(60/100)) #设置滑动条窗口大小

        # image queue dock
        # create a QDockWidget
        imgQueueDockWidget = QDockWidget("img Stack", self)
        imgQueueDockWidget.setObjectName("imgStackDockWidget")
        imgQueueDockWidget.setAllowedAreas(Qt.LeftDockWidgetArea)
        imgQueueDockWidget.setWidget(slider1)
        self.addDockWidget(Qt.LeftDockWidgetArea, imgQueueDockWidget)
        self.windowMenu.addAction(imgQueueDockWidget.toggleViewAction())
        # imgQueueDockWidget.setEnabled(False)

        # camera setting dock
        self.camera_setting = CameraSetting()
        cameraSettingDockWidget = QDockWidget("Camera Setting", self)
        cameraSettingDockWidget.setObjectName("cameraSettingDockWidget")
        cameraSettingDockWidget.setAllowedAreas(Qt.BottomDockWidgetArea)
        cameraSettingDockWidget.setWidget(self.camera_setting)
        self.addDockWidget(Qt.BottomDockWidgetArea, cameraSettingDockWidget)
        self.windowMenu.addAction(cameraSettingDockWidget.toggleViewAction())

        self.camera_setting.prefix_text.editingFinished.connect(self.editFinished)

        # image display setting dock
        self.img_display_setting = ImgDisplaySetting()
        # create a QDockWidget
        displaySettingDockWidget = QDockWidget("Display Setting", self)
        displaySettingDockWidget.setObjectName("displaySettingDockWidget")
        displaySettingDockWidget.setAllowedAreas(Qt.BottomDockWidgetArea)
        displaySettingDockWidget.setWidget(self.img_display_setting)
        self.addDockWidget(Qt.BottomDockWidgetArea, displaySettingDockWidget)
        # enable the toggle view action
        self.windowMenu.addAction(displaySettingDockWidget.toggleViewAction())
        self.img_display_setting.roi_size.editingFinished.connect(self.editFinished2)
        self.img_display_setting.xpos.editingFinished.connect(self.editFinished3)
        self.img_display_setting.ypos.editingFinished.connect(self.editFinished4)
        self.img_display_setting.img_sub.connect(self.plot_main_window.img_plot2)
        self.img_display_setting.img_sub2.connect(self.plot_main_window.img_plot3)
        self.img_display_setting.img_sub.connect(lambda: self.plot_main_window.judrcf(self.img_display_setting.roi, self.img_display_setting.rawdata, self.Fitting_dock.Fitting))
        self.img_display_setting.img_sub2.connect(lambda: self.plot_main_window.judrcf(self.img_display_setting.roi, self.img_display_setting.rawdata, self.Fitting_dock.Fitting))
        # self.img_display_setting.gotoposition.clicked.connect(self.gotopos)
        self.img_display_setting.realtime.stateChanged.connect(lambda: self.realtime(self.img_display_setting.realtime))


        # image analyse setting dock
        self.img_analyse_setting = ImgParameters()
        analyseDataDockWidget = QDockWidget("Image Setup and Laser", self)
        analyseDataDockWidget.setObjectName("analyseDataDockWidget")
        analyseDataDockWidget.setAllowedAreas(Qt.RightDockWidgetArea)
        analyseDataDockWidget.setWidget(self.img_analyse_setting)
        self.addDockWidget(Qt.RightDockWidgetArea, analyseDataDockWidget)
        self.windowMenu.addAction(analyseDataDockWidget.toggleViewAction())

        self.img_analyse_setting.magValue.editingFinished.connect(self.change_mag)
        self.img_analyse_setting.TOF.editingFinished.connect(self.change_TOF)


        # output dock
        # self.prompt_dock = PromptWidget()
        # promptDockWidget = QDockWidget("Output Console", self)
        # promptDockWidget.setObjectName("consoleDockWidget")
        # # promptDockWidget.setAllowedAreas(Qt.RightDockWidgetArea)
        # promptDockWidget.setWidget(self.prompt_dock)
        # self.addDockWidget(Qt.RightDockWidgetArea, promptDockWidget)
        # # redirect print statements to show a copy on "console"
        sys.stdout = Helper.print_redirect()
        sys.stdout.print_signal.connect(self.update_console)
        # # sys.stdout.print_signal.connect(self.update_pro)
        # self.windowMenu.addAction(promptDockWidget.toggleViewAction())

        # Fitting dock 拟合停靠
        self.Fitting_dock = FittingdataWidget()
        fittingDockWidget = QDockWidget("Fitting Data", self)
        fittingDockWidget.setObjectName("FittingDockWidget")
        fittingDockWidget.setAllowedAreas(Qt.BottomDockWidgetArea)
        fittingDockWidget.setWidget(self.Fitting_dock)
        self.addDockWidget(Qt.BottomDockWidgetArea, fittingDockWidget)
        self.windowMenu.addAction(fittingDockWidget.toggleViewAction())

        self.promptwin_dock = PromptWidget()
        promptwin_dock = QDockWidget("Prompt", self)
        promptwin_dock.setObjectName("Prompt")
        promptwin_dock.setAllowedAreas(Qt.BottomDockWidgetArea)
        promptwin_dock.setWidget(self.promptwin_dock)
        self.addDockWidget(Qt.BottomDockWidgetArea, promptwin_dock)
        self.windowMenu.addAction(promptwin_dock.toggleViewAction())

        # result dock
        self.FlurenceD = Flurence()
        self.result_dock = Absorption()
        FlurenceDockWidget = QDockWidget("Flurence Image", self)
        FlurenceDockWidget.setObjectName("resultDockWidget")
        FlurenceDockWidget.setAllowedAreas(Qt.RightDockWidgetArea)
        FlurenceDockWidget.setWidget(self.FlurenceD)
        resultDockWidget = QDockWidget("Absorption Image", self)
        resultDockWidget.setObjectName("resultDockWidget")
        resultDockWidget.setAllowedAreas(Qt.RightDockWidgetArea)
        resultDockWidget.setWidget(self.result_dock)

        self.addDockWidget(Qt.RightDockWidgetArea, FlurenceDockWidget)
        self.windowMenu.addAction(FlurenceDockWidget.toggleViewAction())
        self.addDockWidget(Qt.RightDockWidgetArea, resultDockWidget)
        self.windowMenu.addAction(resultDockWidget.toggleViewAction())
        self.tabifyDockWidget(FlurenceDockWidget, resultDockWidget)

        self.FlurenceD.Detu.editingFinished.connect(self.change_Detu)
        self.FlurenceD.Dia.editingFinished.connect(self.change_Dia)
        self.FlurenceD.ToPwr.editingFinished.connect(self.change_ToPwr)



        ### TOOLBAR MENU ###
        self.expToolbar.setObjectName("ExperimentToolbar")
        self.expToolbar.addAction(self.start_exp_action)
        self.expToolbar.addAction(self.stop_exp_action)
        self.expToolbar.addAction(self.detect_action)
        self.expToolbar.addAction(self.capture_action)

        self.plotToolbar.setObjectName("PlotToolbar")
        self.plotToolbar.addAction(self.clear_img_stack_action)
        self.plotToolbar.addAction(self.clear_main_win_action)


        self.tempque = Queue()
        self.tempdata = []
        self.caltemp = 0


        self.LoadfolderAction = Helper.create_action(self,
                                                           "Load folder",
                                                           slot=self.load_img2stack,
                                                           shortcut=None,
                                                           icon=None,
                                                           tip="Load previous images to image stack from file")

        self.LoadImgAction = Helper.create_action(self,
                                                      "Load Images",
                                                      slot=self.file_load_imgs,
                                                      shortcut=None,
                                                      icon=None,
                                                      tip="Load previous images to image stack")

        self.SaveImgAction = Helper.create_action(self,
                                                       "Save all stack images",
                                                       slot=self.file_save_imgs,
                                                       shortcut=None,
                                                       icon=None,
                                                       tip="Save image stack's images")

        self.SaveMainImgAction = Helper.create_action(self,
                                                  "Save MainWindow's images",
                                                  slot=self.Mainwindowfile_save_imgs,
                                                  shortcut=None,
                                                  icon=None,
                                                  tip="Save MainWidnow's images")

        self.SetpathAction = Helper.create_action(self,
                                                      "Set the default save path",
                                                      slot=self.Setpath,
                                                      shortcut=None,
                                                      icon=None,
                                                      tip="Set the default save path")

        # self.AbsorbImageAction = Helper.create_action(self,
        #                                           "Absorb Image",
        #                                           slot=self.img_analyse_setting.absorb_setting,
        #                                           shortcut=None,
        #                                           icon=None,
        #                                           tip="Processing absorption imaging")

        # self.PrefixsettingAction = Helper.create_action(self,
        #                                               "Prefix setting",
        #                                               slot=self.img_analyse_setting.prefix_setting,
        #                                               shortcut=None,
        #                                               icon=None,
        #                                               tip="Prefix setting")

        self.setcolorAction = Helper.create_action(self,
                                                        "colour setting",
                                                        slot=self.setcolour,
                                                        shortcut=None,
                                                        icon=None,
                                                        tip="colour setting")

        self.camsettingAction = Helper.create_action(self,
                                                   "camera information",
                                                   slot=self.caminformation,
                                                   shortcut=None,
                                                   icon=None,
                                                   tip="camera information")
        
        self.camsettingAction2 = Helper.create_action(self,
                                                   "camera information2",
                                                   slot=self.caminformation2,
                                                   shortcut=None,
                                                   icon=None,
                                                   tip="camera information2")

        self.showbkgAction = Helper.create_action(self,
                                                     "background image",
                                                     slot=self.showbkg,
                                                     shortcut=None,
                                                     icon=None,
                                                     tip="background image")

        # self.tempreature = Helper.create_action(self,
        #                                           "T1",
        #                                           slot=self.tempre,
        #                                           shortcut=None,
        #                                           icon=None,
        #                                           tip="tempreature")
        #
        # self.tempreature2 = Helper.create_action(self,
        #                                         "T2",
        #                                         slot=self.tempre2,
        #                                         shortcut=None,
        #                                         icon=None,
        #                                         tip="tempreature")

        self.openfolderAction = Helper.create_action(self,
                                                "Open folder",
                                                slot=self.Openf,
                                                shortcut=None,
                                                icon=None,
                                                tip="Open folder")

        self.fileMenu.addAction(self.LoadImgAction)
        self.fileMenu.addAction(self.LoadfolderAction)
        self.fileMenu.addAction(self.SaveMainImgAction)
        self.fileMenu.addAction(self.SaveImgAction)
        # self.optionMenu.addAction(self.SetpathAction)
        #modi
        self.cameraMenu.addAction(self.camsettingAction)
        self.cameraMenu.addAction(self.camsettingAction2)
        # self.optionMenu.addAction(self.showbkgAction)
        # self.optionMenu.addAction(self.AbsorbImageAction) #
        self.optionMenu.addAction(self.setcolorAction)
        # self.optionMenu.addAction(self.tempreature)
        # self.optionMenu.addAction(self.tempreature2)
        self.path.addAction(self.SetpathAction)
        self.path.addAction(self.openfolderAction)
        # self.optionMenu.addAction(self.PrefixsettingAction)


        # queue for update main window when camera is in video mode
        self.acquiring = False
        # thread for acquiring image from camera to queue
        self.thread = None
        self.worker = None
        self.connect_slot2signal()

        icpath = IOHelper.get_configt_setting('DATA_PATH')
        # icpath = Path(icpath)
        icpath = str(icpath).replace('\\', '/')
        icpath = icpath + '/icon/UALab.png'
        # print(icpath)
        self.setWindowIcon(QIcon(icpath))
        self.show()

    def gotopos(self):
        self.plot_main_window.goto_pos()

    def realtime(self, mode):
        if mode.isChecked() == True:
            settings.widget_params["Analyse Data Setting"]["realtime"] = True
        else:
            settings.widget_params["Analyse Data Setting"]["realtime"] = False
        self.plot_main_window.setroichange()



    def editFinished(self):
        settings.widget_params["Analyse Data Setting"]["Prefix"] = str(self.camera_setting.prefix_text.text())
        # print(settings.widget_params["Analyse Data Setting"]["Prefix"])
        # print(type(settings.widget_params["Analyse Data Setting"]["Prefix"]))

    def editFinished2(self):
        settings.widget_params["Analyse Data Setting"]["roisize"] = int(self.img_display_setting.roi_size.text())
        self.plot_main_window.change_roisize()
        # print(settings.widget_params["Analyse Data Setting"]["Prefix"])
        # print(type(settings.widget_params["Analyse Data Setting"]["Prefix"]))

    def editFinished3(self):
        settings.widget_params["Analyse Data Setting"]["xpos"] = float(self.img_display_setting.xpos.text())
        self.gotopos()
        # self.plot_main_window.change_roisize()

    def editFinished4(self):
        settings.widget_params["Analyse Data Setting"]["ypos"] = float(self.img_display_setting.ypos.text())
        self.gotopos()
        # self.plot_main_window.change_roisize()

    def change_Detu(self):
        settings.widget_params["Analyse Data Setting"]["Detu"] = float(self.FlurenceD.Detu.text())
        self.plot_main_window.calculate_roi()

    def change_Dia(self):
        settings.widget_params["Analyse Data Setting"]["Dia"] = float(self.FlurenceD.Dia.text())
        self.plot_main_window.calculate_roi()

    def change_ToPwr(self):
        settings.widget_params["Analyse Data Setting"]["ToPwr"] = float(self.FlurenceD.ToPwr.text())
        self.plot_main_window.calculate_roi()

    def change_mag(self):
        settings.widget_params["Analyse Data Setting"]["magValue"] = float(self.img_analyse_setting.magValue.text())
        self.plot_main_window.calculate_roi()

    def change_TOF(self):
        settings.widget_params["Analyse Data Setting"]["TOF"] = float(self.img_analyse_setting.TOF.text())
        self.plot_main_window.measureT()


    def change_camera_params(self):
        self.camera_setting.apply_button.setEnabled(False)
        if self.acquiring:
            self.sig_abort_workers.emit()
            self.thread.quit()  # this will quit **as soon as thread event loop unblocks**
            self.thread.wait()  # <- so you need to wait for it to *actually* quit
            print("camera thread quit")
            self.worker = Worker()
            self.thread = QThread()
            self.worker.moveToThread(self.thread)
            self.worker.sig_video_mode_img.connect(self.update_main_plot_win)
            self.worker.sig_hardware_mode_img.connect(self.update_image_queue)
            # control worker:
            self.sig_abort_workers.connect(self.worker.abort)
            self.thread.started.connect(self.worker.work)
            self.thread.start()  # this will emit 'started' and start thread's event loop
            print("camera setting is applied ")
        self.camera_setting.apply_button.setEnabled(True)

    def change_camera_mode(self, mode):
        if self.acquiring:
            self.sig_abort_workers.emit()
            self.thread.quit()  # this will quit **as soon as thread event loop unblocks**
            self.thread.wait()  # <- so you need to wait for it to *actually* quit
            print("camera thread quit")
        if mode.isChecked():
            # self.sig_abort_workers.emit()
            # self.thread.quit()  # this will quit **as soon as thread event loop unblocks**
            # self.thread.wait()  # <- so you need to wait for it to *actually* quit
            # print("camera thread quit")
            if mode.text() == 'Video':
                settings.widget_params["Image Display Setting"]["mode"] = 0
                self.camera_setting.hardware_mode.setEnabled(True)
                self.camera_setting.video_mode.setEnabled(False)
                self.camera_setting.hardware_mode.setChecked(False)
                self.camera_setting.apply_button.setEnabled(True)
                self.camera_setting.AbsTriger.setEnabled(True)
                self.camera_setting.AbsTriger.setChecked(False)
                self.camera_setting.camera_further_setting.gain_value.setEnabled(True)
                self.camera_setting.camera_further_setting.exposure_time.setEnabled(True)
                self.camera_setting.camera_further_setting.shutter_time.setEnabled(True)
                print('Video mode')
            elif mode.text() == 'ExtTrig':
                settings.widget_params["Image Display Setting"]["mode"] = 1
                self.camera_setting.hardware_mode.setEnabled(False)
                self.camera_setting.video_mode.setChecked(False)
                self.camera_setting.video_mode.setEnabled(True)
                self.camera_setting.apply_button.setEnabled(False)
                self.camera_setting.AbsTriger.setEnabled(True)
                self.camera_setting.AbsTriger.setChecked(False)
                self.camera_setting.camera_further_setting.gain_value.setEnabled(False)
                self.camera_setting.camera_further_setting.exposure_time.setEnabled(False)
                self.camera_setting.camera_further_setting.shutter_time.setEnabled(False)
                # self.img_display_setting.video_mode.setChecked(True)
                print('hardware mode')

            elif mode.text() == 'AbsImg':
                settings.widget_params["Image Display Setting"]["mode"] = 2
                self.camera_setting.AbsTriger.setEnabled(False)
                self.camera_setting.hardware_mode.setEnabled(True)
                self.camera_setting.video_mode.setEnabled(True)
                self.camera_setting.hardware_mode.setChecked(False)
                self.camera_setting.video_mode.setChecked(False)
                print('Absorption imaging mode')

        if self.acquiring:
            self.worker = Worker()
            self.thread = QThread()
            self.worker.moveToThread(self.thread)
            self.worker.sig_video_mode_img.connect(self.update_main_plot_win)
            self.worker.sig_hardware_mode_img.connect(self.update_image_queue)
            # control worker:
            self.sig_abort_workers.connect(self.worker.abort)
            self.thread.started.connect(self.worker.work)
            self.thread.start()  # this will emit 'started' and start thread's event loop
        # # print("camera is in new mode")

    def start_exp(self):
        """
        start basis experiment include capturing images, more operations can be
        added here or use a script file to control instrument accurately.
        :return:
        """
        if settings.instrument_params["Camera"]["index"] is not None:

            self.start_exp_action.setEnabled(False)#CANNOT START

            self.LoadImgAction.setEnabled(False)
            self.LoadfolderAction.setEnabled(False)#CANNOT LOAD AND SAVE
            self.SaveImgAction.setEnabled(False)
            self.SaveMainImgAction.setEnabled(False)
            self.detect_action.setEnabled(False)
            self.camsettingAction.setEnabled(False)


            self.camera_setting.AbsTriger.setEnabled(False)
            self.camera_setting.video_mode.setEnabled(False)
            self.camera_setting.hardware_mode.setEnabled(False)

            self.clear_img_stack_action.setEnabled(False)
            self.clear_main_win_action.setEnabled(False)
            try:
            # if 1:
                print("try1")
                self.worker = Worker()
                print("try1.1")
                self.thread = QThread()
                print("try1.2")
                self.worker.moveToThread(self.thread)
                self.worker.sig_video_mode_img.connect(self.update_main_plot_win)
                self.worker.sig_hardware_mode_img.connect(self.update_image_queue)
                # control worker:
                print("try2")
                self.sig_abort_workers.connect(self.worker.abort)
                self.thread.started.connect(self.worker.work)
                self.thread.start()  # this will emit 'started' and start thread's event loop
                #modi
                print("try3")
                # finish camera index setting, then can't change camera index during experiment,
                # if want to change camera index, then stop experiment
                self.camera_setting.cb.setEnabled(False)
                self.camera_setting.further_setting.setEnabled(True)
                self.camera_setting.apply_button.setEnabled(True)
                settings.widget_params["Image Display Setting"]["imgSource"] = "camera"
                self.acquiring = True
                settings.cameraON = True
                if settings.widget_params["Image Display Setting"]["mode"] == 0:
                    self.img_display_setting.roi.setCheckState(0)
                else:
                    self.plot_main_window.ClearPlotwin()

                self.stop_exp_action.setEnabled(True)
            except:
            # else:
                self.stop_exp()
                print('Wrong camera selection')
        else:
            print("select a camera for further experiment")

    def stop_exp(self):
        """
        stop basis experiment include capturing images when image source is camera.
        :return:
        """
        self.stop_exp_action.setEnabled(False)  #already stop, so connot stop
        if self.acquiring:
            self.sig_abort_workers.emit()
            self.thread.quit()  # this will quit **as soon as thread event loop unblocks**
            self.thread.wait()  # <- so you need to wait for it to *actually* quit


        self.acquiring = False
        settings.cameraON = False
        self.start_exp_action.setEnabled(True)   # already stop can start
        self.LoadfolderAction.setEnabled(True)  #can load
        self.SaveImgAction.setEnabled(True)  #can save
        self.SaveMainImgAction.setEnabled(True)
        self.LoadImgAction.setEnabled(True)
        self.clear_img_stack_action.setEnabled(True)  #can clear stack
        self.clear_main_win_action.setEnabled(True)   #can clear all
        self.detect_action.setEnabled(True)
        self.camsettingAction.setEnabled(True)
        self.camera_setting.cb.setEnabled(True)
        self.camera_setting.further_setting.setEnabled(True)

        if settings.widget_params["Image Display Setting"]["mode"] == 2:
            self.camera_setting.video_mode.setEnabled(True)
            self.camera_setting.hardware_mode.setEnabled(True)
        elif settings.widget_params["Image Display Setting"]["mode"] == 0:
            self.camera_setting.AbsTriger.setEnabled(True)
            self.camera_setting.hardware_mode.setEnabled(True)
        elif settings.widget_params["Image Display Setting"]["mode"] == 1:
            self.camera_setting.AbsTriger.setEnabled(True)
            self.camera_setting.video_mode.setEnabled(True)
        # self.camera_setting.AbsTriger.setEnabled(True)
        # self.camera_setting.video_mode.setEnabled(True)
        # self.camera_setting.hardware_mode.setEnabled(True)

        settings.widget_params["Image Display Setting"]["imgSource"] = "disk"

        # self.camera_setting.video_mode.setChecked(False)
        # self.camera_setting.hardware_mode.setChecked(False)
        # self.camera_setting.video_mode.setEnabled(False)
        # self.camera_setting.hardware_mode.setEnabled(False)

    def detect(self):
        self.camera_setting.detect()

    def queuelkut(self, colormapname):
        from matplotlib import cm
        colormapname = cm.get_cmap(settings.colorlist[colormapname])  # cm.get_cmap("CMRmap")
        colormapname._init()
        lut = (colormapname._lut * 255).view(np.ndarray)
        for i in range(settings.widget_params["Image Display Setting"]["img_stack_num"]):
            plot_win = self.img_queue.plot_wins.get()  #得到 一个作图小窗口
            plot_win.video.setLookupTable(lut)
            self.img_queue.plot_wins.put(plot_win)

    def connect_slot2signal(self):

        # image display widget
        # all parameters' signal are connected to global parameters.
        # self.plot_main_window.tempnum.connect(self.result_dock.change_temp)
        self.plot_main_window.tempnum.connect(self.plot_main_window.change_TEM)
        # self.plot_main_window.temcula.connect(self.addVal)

        self.camera_setting.video_mode.stateChanged.connect(lambda: self.change_camera_mode(self.camera_setting.video_mode))
        self.camera_setting.hardware_mode.stateChanged.connect(lambda: self.change_camera_mode(self.camera_setting.hardware_mode))
        self.camera_setting.AbsTriger.stateChanged.connect(lambda: self.change_camera_mode(self.camera_setting.AbsTriger))

        # image stack widget
        for i in range(settings.widget_params["Image Display Setting"]["img_stack_num"]):
            plot_win = self.img_queue.plot_wins.get()  #得到 一个作图小窗口
            plot_win.img_dict.connect(self.plot_main_window.img_plot) #把发射信号源和接受处连接，点击堆栈小窗口后执行主窗口的显示图像
            plot_win.img_dict.connect(lambda: self.plot_main_window.judrcf(self.img_display_setting.roi, self.img_display_setting.rawdata, self.Fitting_dock.Fitting))
            self.img_queue.plot_wins.put(plot_win)
        # plot main window widget
        self.plot_main_window.atom_number.connect(self.result_dock.change_atom_num)
        self.plot_main_window.Pxatom_num.connect(self.result_dock.change_Pxatom_num)
        self.plot_main_window.atom_numberF.connect(self.FlurenceD.change_atom_num)
        self.plot_main_window.Pxatom_numF.connect(self.FlurenceD.change_Pxatom_num)
        # self.plot_main_window.TotalPhotons_num.connect(self.result_dock.change_TotalPhotons_num)
        self.plot_main_window.roipos.connect(self.img_display_setting.change_roiposition)
        self.plot_main_window.fittingdata.connect(self.Fitting_dock.change_label)
        # self.plot_main_window.HOD.connect(self.result_dock.change_HpeakOD)
        # self.plot_main_window.VOD.connect(self.result_dock.change_VpeakOD)
        self.plot_main_window.HOD.connect(self.plot_main_window.change_HpeakOD)
        self.plot_main_window.VOD.connect(self.plot_main_window.change_VpeakOD)
        self.Fitting_dock.fitting_jud.connect(lambda: self.plot_main_window.add_fitting(self.Fitting_dock.Fitting))
        self.Fitting_dock.fitting_jud.connect(self.plot_main_window.update_ch_fitting_cs)

        # analyse data widget
        self.img_display_setting.roi.stateChanged.connect(lambda: self.plot_main_window.add_roi(self.img_display_setting.roi, self.img_display_setting.rawdata))
        self.img_display_setting.roi.stateChanged.connect(lambda: self.changefit(self.img_display_setting.roi))
        self.img_display_setting.rawdata.stateChanged.connect(lambda: self.plot_main_window.add_rawdata(self.img_display_setting.rawdata))
        self.img_display_setting.colormap.activated.connect(lambda: self.plot_main_window.colormapsel(self.img_display_setting.colormap.currentIndex()))
        self.img_display_setting.colormap.activated.connect(lambda: self.queuelkut(self.img_display_setting.colormap.currentIndex()))
        # self.camera_setting.cross_axes.stateChanged.connect(
        #     lambda: self.plot_main_window.add_cross_axes(self.camera_setting.cross_axes))
        # self.camera_setting.AbsTriger.stateChanged.connect(
        #     lambda: self.absclick(self.camera_setting.AbsTriger))
        self.camera_setting.auto_save.stateChanged.connect(
            lambda: self.autosave(self.camera_setting.auto_save))

        # camera setting widget
        self.camera_setting.apply_button.clicked.connect(self.camera_setting.camera_further_setting.change_exposure)
        self.camera_setting.apply_button.clicked.connect(self.camera_setting.camera_further_setting.change_gain)
        self.camera_setting.apply_button.clicked.connect(self.camera_setting.camera_further_setting.change_shutter)
        self.camera_setting.apply_button.clicked.connect(self.change_camera_params)

    def changefit(self, mode):
        if mode.isChecked():
            pass
        else:
            self.Fitting_dock.Fitting.setCheckState(0)

    def autosave(self,auto_save_state):
        if auto_save_state.isChecked():
            settings.widget_params["Analyse Data Setting"]["autoStatus"] = True

        else:
            settings.widget_params["Analyse Data Setting"]["autoStatus"] = False
        # print(settings.widget_params["Analyse Data Setting"]["autoStatus"])


    def clear_img_stack(self):
        """
        clear image stack
        :return:
        """
        if self.acquiring and settings.widget_params["Image Display Setting"]["mode"] == 0:
            print("video mode can't clear image stack")
            return
        # make sure that queue isn't changing when using qsize()
        for i in range(settings.widget_params["Image Display Setting"]["img_stack_num"]):
            plot_win = self.img_queue.plot_wins.get()
            plot_win.clear_win()
            self.img_queue.plot_wins.put(plot_win)

        settings.absimgData[0] = []  # Erase the last data
        settings.absimgData[1] = []
        settings.absimgData[2] = []
        settings.absimgData[3] = []

    def clear_main_win(self):
        """
              clear main windows
              :return:
              """
        if self.acquiring and settings.widget_params["Image Display Setting"]["mode"] == 0:
            print("video mode can't clear main window")
            return
        self.img_display_setting.roi.setCheckState(0)
        self.plot_main_window.clear_win()
        settings.imgData["Img_photon_range"] = []
        settings.imgData["Img_data"] = []

    def capture_image(self):
        """
        Capture the current image
        :return:
        """
        if self.plot_main_window.img.image is None:
            print("have no image in Mainwindow")
            return
        img_data = np.array(self.plot_main_window.img.image)
        # load image name by path
        img_name2 = (self.plot_main_window.img_label.text())[0:25].replace(' ', '_').replace(':', '').replace('-', '').replace('.', '_')
        #print('def capture_image(self):img_name2',img_name2)
        img_name = str(img_name2)
        img_cap = {'img_data':img_data, 'img_name':img_name}
        self.update_image_queue(img_cap)

    ### LOAD CUSTOM SETTING FOR INSTRUMENT CONNECT AND PARAMETERS ###

    def file_save_imgs(self):
        """
        save image stack's images to disk
        :return:
        """
        # try:
        fpath = IOHelper.get_config_setting('DATA_PATH')
        fpath = Path(fpath)
        dir_path = fpath
        # dir_path = fpath.joinpath(str(datetime.datetime.now()).split('.')[0].replace(' ', '-').replace(':', '_'))
        # print("save images to {}".format(dir_path))
        if settings.m_path != []:
            dir_path = settings.m_path
        if not dir_path.exists():
            dir_path.mkdir()
        for i in range(settings.widget_params["Image Display Setting"]["img_stack_num"]):
            plot_win = self.img_queue.plot_wins.get()
            if plot_win.video.image is not None:
                img_data = np.array(plot_win.video.image)
                # load image name by path
                img_name2 = (plot_win.img_label)[0:25].replace(' ', '_').replace(':', '').replace('-', '').replace('.', '_')
                img_name = str(img_name2)
                # img_data = img_data[::-1]
                # img_data = Image.fromarray(img_data)
                # img_data.save(r"{}\{}.png".format(dir_path, img_name))
                import numpy
                numpy.save(r"{}\{}".format(dir_path, img_name), img_data)
            self.img_queue.plot_wins.put(plot_win)
        print("save images to {}".format(dir_path))
        # print("images have saved.")
        # except OSError:
        #     print("Only new version files can be saved.")

    def Openf(self):
        import os
        fpath = IOHelper.get_config_setting('DATA_PATH')
        if settings.m_path != []:
            fpath = settings.m_path
        # fpath = Path(fpath)
        os.system("start explorer %s" % fpath)

    def Setpath(self):
        mpath = QFileDialog.getExistingDirectory(self, "Set path")
        if len(mpath) > 1:
            settings.m_path = Path(mpath)
            self.path.setTitle('##  Save file to: ' + str(settings.m_path) + '  ##')
            # print(settings.m_path)

    def caminformation2(self):
        
        print("caminformation2")

    def caminformation(self):
        self.caminfdia = QDialog()  # create a dialog
        existingdata = get_camconfig_setting()
        label1 = QLabel('Camera')
        label2 = QLabel('nick name')
        label3 = QLabel('shutter time')
        layt = QGridLayout()
        layt.addWidget(label1,0,0,1,1)
        layt.addWidget(label2,0,1,1,1)
        layt.addWidget(label3,0,2,1,1)
        self.camera_infor = Queue()
        xx = 1
        for n in range(int(len(existingdata)/3)):
            name = QLabel(str(existingdata[3*n]))
            nickname = QLineEdit(str(existingdata[3*n+1]))
            data = QLineEdit(str(existingdata[3*n+2]))
            data.setValidator(QIntValidator(data))
            data.setMaxLength(3)
            # nickname.editingFinished.connect(self.reset_caminf)
            # data.editingFinished.connect(self.reset_caminf)
            layt.addWidget(name,n+1,0,1,1)
            layt.addWidget(nickname,n+1,1,1,1)
            layt.addWidget(data,n+1,2,1,1)
            self.camera_infor.put(name)
            self.camera_infor.put(nickname)
            self.camera_infor.put(data)
            xx = xx+1
        labtip = QLabel('Click apply to change the Settings')
        labtip.setAlignment(Qt.AlignRight)
        button1 = QPushButton('apply')
        layt.addWidget(labtip,xx,0,1,2)
        layt.addWidget(button1,xx,2,1,1)
        button1.clicked.connect(self.reset_caminf)
        self.caminfdia.setFont(QFont("Roman times", 10))

        self.caminfdia.setLayout(layt)
        self.caminfdia.setWindowTitle('Camera information')
        # self.caminfdia.setWindowModality(Qt.ApplicationModal)
        self.caminfdia.exec_()


    def showbkg(self):
        self.showbkgdia = QDialog()  # create a dialog
        self.layt = QVBoxLayout()
        import pyqtgraph as pg
        from pyqtgraph import GraphicsLayoutWidget
        pg.setConfigOptions(imageAxisOrder='row-major')
        self.viewport = GraphicsLayoutWidget()
        self.video_view = self.viewport.addViewBox()
        self.video_view.clear()
        self.video = pg.ImageItem()
        self.video_view.addItem(self.video)
        self.video_view.setMouseEnabled(x=False, y=False)  # make it can not move
        self.layt.addWidget(self.viewport)
        if settings.imgData["BkgImg"] != []:
            self.video.setImage(settings.imgData["BkgImg"])
        self.showbkgdia.setLayout(self.layt)
        self.showbkgdia.setWindowTitle('background image')
        # self.caminfdia.setWindowModality(Qt.ApplicationModal)
        self.showbkgdia.exec_()

    def reset_caminf(self):
        config_path = IOHelper.CAMCONFIGT_FILE_PATH
        caminfo = open(config_path, 'w')

        for m in range(self.camera_infor.qsize()):
            tmpvar = self.camera_infor.get()
            # print(tmpvar.text())
            if m % 3 == 2:
                caminfo.write("{}\n".format(tmpvar.text()))
            else:
                caminfo.write("{},".format(tmpvar.text()))
            self.camera_infor.put(tmpvar)

        caminfo.close()
        self.camera_setting.detect()

    def setcolour(self):
        self.setcol = QDialog()  # create a dialog
        lay = QVBoxLayout()
        self.black = QPushButton("dark", self)
        self.black.clicked.connect(setcoldark)
        self.white = QPushButton("bright", self)
        self.white.clicked.connect(setcollight)
        lay.addWidget(self.black)
        lay.addWidget(self.white)
        self.setcol.setLayout(lay)
        self.setcol.setWindowTitle('colour setting')
        self.setcol.setWindowModality(Qt.NonModal)
        self.setcol.setWindowFlags(Qt.WindowStaysOnTopHint)
        self.setcol.show()


    def Mainwindowfile_save_imgs(self):
        # try:
        if self.plot_main_window.img.image is None:
            print("have no image in Mainwindow")
            return
        fpath = IOHelper.get_config_setting('DATA_PATH')
        fpath = Path(fpath)
        dir_path = fpath
        # dir_path = fpath.joinpath(str(datetime.datetime.now())[2:].split('.')[0].replace(' ', '-').replace(':', '_'))
        # print("save images to {}".format(dir_path))
        if settings.m_path != []:
            dir_path = settings.m_path
        if not dir_path.exists():
            dir_path.mkdir()
        img_data = np.array(self.plot_main_window.img.image)
        # load image name by path
        img_name2 = (self.plot_main_window.img_label.text())[0:25].replace(' ', '_').replace(':', '').replace('-', '').replace('.', '_')
        img_name = str(img_name2)
        # img_data = img_data[::-1]
        # img_data = Image.fromarray(img_data)
        # img_data.save(r"{}\{}.png".format(dir_path, img_name))
        save(r"{}\{}".format(dir_path, img_name), img_data)
        print("save images to {}".format(dir_path))
        # print("images have saved.")
        # except OSError:
        #     print('Only new version files can be saved.')


    def file_load_imgs(self):
        """
        Load previous image to stack.
        :return:
        """
        print("file_load_imgs")
        settings.widget_params["Image Display Setting"]["imgSource"] = "disk"
        fpath = IOHelper.get_config_setting('DATA_PATH')
        img_fpath = QFileDialog.getOpenFileNames(self, "Open File", fpath, "Image files(*.npy *.data *.tif)")  # name path

        if img_fpath[0] != '':
            for i in range(len(img_fpath[0])):
                plot_win = self.img_queue.plot_wins.get()
                try:
                    img_file = str(img_fpath[0][i])
                    img_paths = Path(img_file)
                    plot_win.img_plot(self.load_img_dict(img_paths))
                    self.img_queue.plot_wins.put(plot_win)
                except TypeError:
                    return
                except PermissionError:
                    return
                

    def load_img2stack(self):#load single picture
        """
        load images to image queue, with image name and data
        """
        settings.widget_params["Image Display Setting"]["imgSource"] = "disk"
        fpath = IOHelper.get_config_setting('DATA_PATH')

        img_fpath = QFileDialog.getExistingDirectory(self, "Open File", fpath)  # name path
        img_file = Path(img_fpath)
        # img_path1 = list(img_file.glob('*.png'))
        img_path2 = list(img_file.glob('*.npy'))
        img_path3 = list(img_file.glob('*.data'))
        img_path4 = list(img_file.glob('*.tif'))
        img_paths = img_path2 + img_path3 + img_path4

        for win_index in range(settings.widget_params["Image Display Setting"]["img_stack_num"]):
            if win_index == len(img_paths):
                break
            if img_paths != []:
                plot_win = self.img_queue.plot_wins.get()
                plot_win.img_plot(self.load_img_dict(img_paths[win_index]))
                self.img_queue.plot_wins.put(plot_win)

    ### MISCELLANY ###

    def load_img_dict(self, img_path):
        pathjud = str(img_path)
        pathjud = pathjud[len(pathjud) - 3:]   #Get the version of the file

        if pathjud == 'ata':
            file = open(img_path)
            linescontent = file.readlines()                     #Read the file as a behavior unit
            rows = len(linescontent)                            #get the numbers fo line
            lines = len(linescontent[0].strip().split(' '))
            # print(rows)
            # print(lines)
            img_data = np.zeros((rows, lines))                  # Initialization matrix
            row = 0
            for line in linescontent:
                line = line.strip().split(' ')
                img_data[row, :] = line[:]
                row += 1
            file.close()
        else:
            img_data = np.load(img_path,encoding='bytes')

        # img_data = img_data[::-1]
        img_name = img_path.stem
        img = {
            'img_name': img_name,
            'img_data': img_data
        }
        return img


    def update_console(self, stri):
        MAX_LINES = 50
        # localtime = time.struct_time(time.localtime(time.time()))
        # PRtime = str(localtime.tm_hour) + ':' + str(localtime.tm_min) + ':' + str(localtime.tm_sec) + ' '
        stri = str(stri)
        # new_text = self.result_dock.console_text() + '\n' + stri
        # line_list = new_text.splitlines()
        # N_lines = min(MAX_LINES, len(line_list))
        # # limit output lines
        # new_text = '\n'.join(line_list[-N_lines:])
        # self.result_dock.console_text(new_text)
        # self.result_dock.automatic_scroll()
        # self.prompt.setTitle(stri)

        new_text = self.promptwin_dock.console_text() + '\n' + stri
        line_list = new_text.splitlines()
        N_lines = min(MAX_LINES, len(line_list))
        # limit output lines
        new_text = '\n'.join(line_list[-N_lines:])
        self.promptwin_dock.console_text(new_text)
        self.promptwin_dock.automatic_scroll()


    def update_main_plot_win(self, img_dict): #video_mode do this
        """
        Updates the main plot window at regular intervals. It designs for video mode
        """
        # take the newest image in the queue
        QApplication.processEvents()
        settings.imgData["Img_data"] = []
        if img_dict is None:
            return
        self.plot_main_window.img_plot(img_dict)
        print("update_main_plot_win")
        settings.imgData["Img_data"] = img_dict['img_data']

    def update_image_queue(self, img_dict):   #hardware_mode do this
        # QApplication.processEvents()
        plot_win = self.img_queue.plot_wins.get()
        plot_win.img_plot(img_dict)
        print(plot_win.img_label) #modi
        
        img_name2 = (plot_win.img_label)[0:].replace(' ', '_').replace(':', '').replace('-', '').replace('.', '_')
        
        self.img_queue.plot_wins.put(plot_win)
        print("update image queue")
        if settings.widget_params["Analyse Data Setting"]["autoStatus"] == True:
            QApplication.processEvents()
            fpath = IOHelper.get_config_setting('DATA_PATH')
            fpath = Path(fpath)
            dir_path = fpath
            # dir_path = fpath.joinpath(str(datetime.datetime.now())[2:].split('.')[0].replace(' ', '').replace(':', '_'))
            if settings.m_path != []:
                dir_path = settings.m_path
            if not dir_path.exists():
                dir_path.mkdir()
            img_data = np.array(img_dict['img_data'])
            # load image name by path
            # img_name2 = (self.img_label)[0:20].replace(' ', '~').replace(':', '').replace('-', '')
            
            #modi
            settings.StorageNum = settings.StorageNum + 1
            img_name = str(settings.StorageNum)+str(img_name2)
            # img_name = str(img_name2)
            # img_data = img_data[::-1]
            # from numpy import savetxt
            # save(r"{}\{}".format(dir_path, img_name), img_data)

            save(r"{}\{}".format(dir_path, img_name), img_data)
            print("save images 2 {}".format(dir_path))
        # print("update image queue")

    # def setcol1(self):
    #     # app = QApplication(sys.argv)
    #     palette = QPalette() #调色板
    #     palette.setColor(QPalette.Window, QColor(53, 53, 53))
    #     palette.setColor(QPalette.WindowText, Qt.white)
    #     palette.setColor(QPalette.Base, QColor(25, 25, 25))
    #     palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    #     palette.setColor(QPalette.ToolTipBase, Qt.white)
    #     palette.setColor(QPalette.ToolTipText, Qt.white)
    #     palette.setColor(QPalette.Text, Qt.white)
    #     palette.setColor(QPalette.Button, QColor(53, 53, 53))
    #     palette.setColor(QPalette.ButtonText, Qt.white)
    #     palette.setColor(QPalette.BrightText, Qt.red)
    #     palette.setColor(QPalette.Link, QColor(42, 130, 218))
    #     palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    #     palette.setColor(QPalette.HighlightedText, Qt.black)
    #     app.setPalette(palette)
    #
    # def setcol2(self):
    #     # app = QApplication(sys.argv)
    #     palette = QPalette() #调色板
    #     palette.setColor(QPalette.Window, QColor(232, 232, 232))
    #     app.setPalette(palette)

class Worker(QObject):
    """
    Must derive from QObject in order to emit signals, connect slots to other signals, and operate in a QThread.
    """

    sig_video_mode_img = pyqtSignal(dict)
    sig_hardware_mode_img = pyqtSignal(dict)

    

    #start modi
    def __init__(self):
        super().__init__()
        ## try:
        ### self.camera = Chameleon()
        ###self.camera.initializeCamera(int(settings.instrument_params["Camera"]["index"]))
        print("entered worker")
        
        print("setting params")
        if settings.widget_params["Image Display Setting"]["mode"] == 0:
            mode = 0
        else:
            mode = 2
        ###self.camera.setAcquisitionMode(mode)
        # self.camera.setGRAB_MODE()
        # self.camera.setExposure(settings.instrument_params["Camera"]["exposure time"])
        ###self.camera.setShutter(float(settings.instrument_params["Camera"]["shutter time"]))
        # self.camera.setGain(settings.instrument_params["Camera"]["gain value"])
        # set a low grab timeout to avoid crash when retrieve image.
        ###self.camera.set_grab_timeout(grab_timeout=10)
        ###

        # TUCAM_Dev_Open = TUSDKdll.TUCAM_Dev_Open
        # TUCAMOPEN = TUCAM_OPEN(0, 0)
        # TUCAM_Dev_Open(pointer(TUCAMOPEN));###
        

        # print("TUCAM_Reg_Read")
        
        # # SN
        # TUCAM_Reg_Read = TUSDKdll.TUCAM_Reg_Read
        # cSN = (c_char * 64)() 
        # pSN = cast(cSN, c_char_p)
        # TUCAMREGRW = TUCAM_REG_RW(1, pSN, 64)
        # TUCAM_Reg_Read(c_int64(TUCAMOPEN.hIdxTUCam), TUCAMREGRW)
        # #print(bytes(bytearray(cSN)))
        # print(string_at(pSN))                  ###


        # m_frame = TUCAM_FRAME()
        # m_fs    = TUCAM_FILE_SAVE() 
        # m_format = TUIMG_FORMATS
        # m_frformat= TUFRM_FORMATS
        # m_capmode = TUCAM_CAPTURE_MODES
        # TUCAM_Buf_Alloc = TUSDKdll.TUCAM_Buf_Alloc
        # TUCAM_Cap_Start = TUSDKdll.TUCAM_Cap_Start
##################

        # self.camera.setGRAB_MODE()
        # self.camera.setExposure(settings.instrument_params["Camera"]["exposure time"])
        ###self.camera.setShutter(float(settings.instrument_params["Camera"]["shutter time"]))
        # self.camera.setGain(settings.instrument_params["Camera"]["gain value"])
        # set a low grab timeout to avoid crash when retrieve image.
        ###self.camera.set_grab_timeout(grab_timeout=10)
        ###


        # TUCAM_Buf_WaitForFrame = TUSDKdll.TUCAM_Buf_WaitForFrame

        # TUCAM_Buf_AbortWait = TUSDKdll.TUCAM_Buf_AbortWait
        # TUCAM_Cap_Stop = TUSDKdll.TUCAM_Cap_Stop
        # TUCAM_Buf_Release = TUSDKdll.TUCAM_Buf_Release
        # TUCAM_File_SaveImage = TUSDKdll.TUCAM_File_SaveImage
        # m_fs.nSaveFmt = m_format.TUFMT_TIF.value 

        # m_frame.pBuffer     = 0;
        # m_frame.ucFormatGet = m_frformat.TUFRM_FMT_RAW.value;
        # m_frame.uiRsdSize   = 1;         #需要获取的帧数，1                 
        #print(m_frame.pBuffer)
        #print(m_frame.ucFormatGet) ###




        self.__abort = False
        # a
        # except:
        #     print('Wrong camera selection')
        #     self.stop_exp()

    # @pyqtSlot()
    def work(self):
        # try:
            print("camera start workt")

            # SN
            # Path = './'
            # TUCAM_Api_Init = TUSDKdll.TUCAM_Api_Init
            # TUCAMINIT = TUCAM_INIT(0, Path.encode('utf-8'))
            # TUCAM_Api_Init(pointer(TUCAMINIT));
            TUCAM_Dev_Open = TUSDKdll.TUCAM_Dev_Open
            TUCAMOPEN = TUCAM_OPEN(0, 0)
            TUCAM_Dev_Open(pointer(TUCAMOPEN));

            TUCAM_Reg_Read = TUSDKdll.TUCAM_Reg_Read
            cSN = (c_char * 64)() 
            pSN = cast(cSN, c_char_p)
            TUCAMREGRW = TUCAM_REG_RW(1, pSN, 64)
            TUCAM_Reg_Read(c_int64(TUCAMOPEN.hIdxTUCam), TUCAMREGRW)
            #print(bytes(bytearray(cSN)))

            # Save Image
            m_frame = TUCAM_FRAME()
            m_fs    = TUCAM_FILE_SAVE() 
            m_format = TUIMG_FORMATS
            m_frformat= TUFRM_FORMATS
            m_capmode = TUCAM_CAPTURE_MODES

            TUCAM_Buf_Alloc = TUSDKdll.TUCAM_Buf_Alloc
            TUCAM_Cap_Start = TUSDKdll.TUCAM_Cap_Start
            TUCAM_Buf_WaitForFrame = TUSDKdll.TUCAM_Buf_WaitForFrame
            TUCAM_Buf_AbortWait = TUSDKdll.TUCAM_Buf_AbortWait
            TUCAM_Cap_Stop = TUSDKdll.TUCAM_Cap_Stop
            TUCAM_Buf_Release = TUSDKdll.TUCAM_Buf_Release
            TUCAM_File_SaveImage = TUSDKdll.TUCAM_File_SaveImage
            m_fs.nSaveFmt = m_format.TUFMT_TIF.value 

            m_frame.pBuffer     = 0;
            m_frame.ucFormatGet = m_frformat.TUFRM_FMT_RAW.value;
            m_frame.uiRsdSize   = 1;                          
            print(m_frame.pBuffer)
            print(m_frame.ucFormatGet)
            ###self.camera.startAcquisition()
            ###
            TUCAM_Buf_Alloc(c_int64(TUCAMOPEN.hIdxTUCam), pointer(m_frame))
            #modi
            TUCAM_Cap_Start(c_int64(TUCAMOPEN.hIdxTUCam), m_capmode.TUCCM_TRIGGER_STANDARD.value)
            ###


            if settings.widget_params["Image Display Setting"]["mode"] == 2:
                while True:
                    if self.__abort:
                        break
                    # TestMainWindow().stop_exp_action.setEnabled(False)  # already stop, so connot stop
                    for i in range(3):
                        settings.absimgData[i] = []  #Erase the last data

                    for i in range(3):
                        while True:
                            QApplication.processEvents()  # this could cause change to self.__abort
                            if self.__abort:
                                break
                            #####################
                            #modifrom: img_data = self.camera.retrieveOneImg()  # retrieve image from camera buffer
                            TUCAM_Buf_WaitForFrame(c_int64(TUCAMOPEN.hIdxTUCam), pointer(m_frame))
                            img_name2 = "Image"
                            ImgName = './tif/' + str(img_name2)
                            # ImgName = './' + str(img_name2)
                            print("ImgName is",ImgName)
                            m_fs.pFrame = pointer(m_frame);
                            m_fs.pstrSavePath = ImgName.encode('utf-8');
                            TUCAM_File_SaveImage(c_int64(TUCAMOPEN.hIdxTUCam), m_fs)
                            # 读取tif文件
                            tif_file = "tif\\" + str(img_name2) + ".tif"
                            tif_image = Image.open(tif_file)
                            # 将tif文件转换为numpy数组
                            # tif_array.append(np.array(tif_image))
                            img_data = np.array(tif_image)
                            ############################


                            if img_data is None:
                                continue
                            else:
                                #modi
                                #img_name1 = settings.widget_params["Analyse Data Setting"]["Prefix"]
                                                               
                                if i == 0:
                                    print('i==0')
                                    img_name1 = 'W'
                                elif i == 1:
                                    img_name1 = 'WO'
                                else:
                                    img_name1 = 'BG'
                                timestamp = datetime.datetime.now()
                                #modi
                                # print('timestamp:',str(timestamp)[2:19])
                                #modi
                                # self.sig_hardware_mode_img.emit({'img_name': str(timestamp)[2:19]+str(img_name1), 'img_data': Helper.split_list(img_data)})
                                settings.StackNum = settings.StackNum + 1
                                # self.sig_hardware_mode_img.emit({'img_name': str(settings.StackNum)+str(img_name1), 
                                #                                  'img_data': img_data})
                                ###################
                                # btn_state()

                                ###################
                                img_dict = {'img_data': img_data, 
                                            'img_name': str(settings.StackNum)+str(img_name1)}
                                img_name = str(settings.StackNum)+str(img_name1)
                                # img_name = img_path.stem
                                img = {
                                    'img_name': img_name,
                                    'img_data': img_data}
                                # settings.absimgDatas[self.myserial] = img_data
                                plotMainWindow = PlotMainWindow()
                                plotMainWindow.img_plot(img_dict)
                                # plotMainWindow.btn_state()


                                # import pyqtgraph as pg
                                # self.video = pg.ImageItem()
                                # img_dict = {'img_data': img_data, 
                                #             'img_name': str(settings.StackNum)+str(img_name1)}
                                # self.img_label = img_dict['img_name']
                                
                                # settings.imgData["Img_data"] = img_dict['img_data']
                                # img_dict = pyqtSignal(object)
                                # plotWindow = PlotWindow()
                                # plotWindow.img_dict.emit(img_dict)
                                ###################
                                # self.sig_hardware_mode_img.emit({'img_name': str(settings.StackNum)+str(img_name1), 'img_data': Helper.split_list(img_data)})
                                settings.absimgData[i] = img_data
                                # settings.absimgData[i] = Helper.split_list(img_data)
                                #modi
                                #print('yisiabs',img_name1)
                                break

                    if settings.absimgData[0] != [] and settings.absimgData[1] !=[] and settings.absimgData[2] != []:

                        print('In the calculation')
                        settings.absimgData[3] = Calc_absImg(settings.absimgData[0], settings.absimgData[1], settings.absimgData[2], 0)

                        # modi OD=W
                        # settings.absimgData[3] = settings.absimgData[0]

                        # modi 从窗口中读取命名
                        # img_name1 = settings.widget_params["Analyse Data Setting"]["Prefix"]

                        # modi 常用命名方式
                        # img_name1 = 'OD'
                        # timestamp = datetime.datetime.now()
                        # self.sig_hardware_mode_img.emit({'img_name': str(timestamp)[2:19]+str(img_name1), 'img_data': settings.absimgData[3]})

                        # modi 测试用
                        img_name1 = 'OD'
                        settings.StackNum = settings.StackNum + 1
                        self.sig_hardware_mode_img.emit({'img_name': str(settings.StackNum)+str(img_name1), 
                                                         'img_data': settings.absimgData[3]})
                        ###################
                        # self.video = pg.ImageItem()
                        # self.img_label = img_dict['img_name']
                        # img_dict = {'img_data': np.array(self.video.image), 'img_name': self.img_label}
                        # settings.imgData["Img_data"] = img_dict['img_data']
                        # self.img_dict.emit(img_dict)
                        ###################



            else:
                while True:
                    # check if we need to abort the loop; need to process events to receive signals;
                    # 检查我们是否需要中止回路;需要处理事件以接收信号;
                    QApplication.processEvents()  # this could cause change to self.__abort
                    if self.__abort:    #中断             
                        break

                    img_data = self.camera.retrieveOneImg()  # retrieve image from camera buffer从相机缓冲区检索图像
                    if img_data is None:
                        continue
                    else:
                        img_name1 = settings.widget_params["Analyse Data Setting"]["Prefix"]
                        timestamp = datetime.datetime.now()
                        if settings.widget_params["Image Display Setting"]["mode"] == 1:
                            self.sig_hardware_mode_img.emit({'img_name': str(img_name1)+str(timestamp)[2:], 'img_data': Helper.split_list(img_data)})
                        else:
                            self.sig_video_mode_img.emit({'img_name': str(img_name1)+str(timestamp)[2:], 'img_data': Helper.split_list(img_data)})
                            # set a appropriate refresh value，设置适当的刷新值
                            time.sleep(0.06)#延迟执行给定的秒数。参数可以是一个浮点数，以保证次秒精度。
                            #mark
            ###########self.camera.stopCamera()
            TUCAM_Buf_AbortWait(c_int64(TUCAMOPEN.hIdxTUCam));
            TUCAM_Cap_Stop(c_int64(TUCAMOPEN.hIdxTUCam));
            TUCAM_Buf_Release(c_int64(TUCAMOPEN.hIdxTUCam));

            # CloseCamera
            TUCAM_Dev_Close = TUSDKdll.TUCAM_Dev_Close
            TUCAM_Dev_Close(c_int64(TUCAMOPEN.hIdxTUCam))
            TUCAM_Api_Uninit = TUSDKdll.TUCAM_Api_Uninit
            TUCAM_Api_Uninit
            ####################

        # except:
        #     print(1)

    def abort(self):
        self.__abort = True


def start_main_win():
    global app
    app = QApplication(sys.argv)

    # Force the style to be the same on all OSs:
    app.setStyle("Fusion")

    # palette = QPalette()  # 调色板
    # palette.setColor(QPalette.Window, QColor(53, 53, 53))
    # palette.setColor(QPalette.WindowText, Qt.white)
    # palette.setColor(QPalette.Base, QColor(25, 25, 25))
    # palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    # palette.setColor(QPalette.ToolTipBase, Qt.white)
    # palette.setColor(QPalette.ToolTipText, Qt.white)
    # palette.setColor(QPalette.Text, Qt.white)
    # palette.setColor(QPalette.Button, QColor(53, 53, 53))
    # palette.setColor(QPalette.ButtonText, Qt.white)
    # palette.setColor(QPalette.BrightText, Qt.red)
    # palette.setColor(QPalette.Link, QColor(42, 130, 218))
    # palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    # palette.setColor(QPalette.HighlightedText, Qt.black)
    # app.setPalette(palette)

    light_palette = QPalette()
    # base
    light_palette.setColor(QPalette.WindowText, QColor(0, 0, 0))
    light_palette.setColor(QPalette.Button, QColor(240, 240, 240))
    # light_palette.setColor(QPalette.Light, QColor(180, 180, 180))
    # light_palette.setColor(QPalette.Midlight, QColor(200, 200, 200))
    # light_palette.setColor(QPalette.Dark, QColor(225, 225, 225))
    light_palette.setColor(QPalette.Text, QColor(0, 0, 0))
    light_palette.setColor(QPalette.BrightText, QColor(0, 0, 0))
    light_palette.setColor(QPalette.ButtonText, QColor(0, 0, 0))
    light_palette.setColor(QPalette.Base, QColor(237, 237, 237))
    light_palette.setColor(QPalette.Window, QColor(240, 240, 240))
    # light_palette.setColor(QPalette.Shadow, QColor(20, 20, 20))
    light_palette.setColor(QPalette.Highlight, QColor(76, 163, 224))
    light_palette.setColor(QPalette.HighlightedText, QColor(0, 0, 0))
    light_palette.setColor(QPalette.Link, QColor(0, 162, 232))
    light_palette.setColor(QPalette.AlternateBase, QColor(225, 225, 225))
    light_palette.setColor(QPalette.ToolTipBase, QColor(240, 240, 240))
    light_palette.setColor(QPalette.ToolTipText, QColor(0, 0, 0))

    app.setPalette(light_palette)


    app.setApplicationName("AtomimgSCNU")
    window = MainWindow()
    window.show()
    # print('建议关闭相机后再进行数据处理，以保证操作流畅（正在想办法解决这个问题）')
    # print('It is recommended to turn off the camera before data processing to ensure smoothness')
    print('关闭相机后再进行数据处理')
    print('It is recommended to turn off the camera before data processing')
    
    
    from PyQt5.QtCore import QFileSystemWatcher

                
    class AutoLoadImg:
        def __init__(self):
            self.watcher = QFileSystemWatcher()
            self.watcher.directoryChanged.connect(self.load_latest_img)

        def start_watching(self, dir_path):
            self.dir_path = dir_path
            self.watcher.addPath(dir_path)

        def load_latest_img(self):
            MainWindow.file_load_imgs()
    #modi413
    img_loader = AutoLoadImg()
    dir_path1 = "D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\test\\testimage"  # 监听的图片文件夹路径
    img_loader.start_watching(dir_path1)
    sys.exit(app.exec_())



def setcollight(self):
    light_palette = QPalette()
    # base
    light_palette.setColor(QPalette.WindowText, QColor(0, 0, 0))
    light_palette.setColor(QPalette.Button, QColor(240, 240, 240))
    # light_palette.setColor(QPalette.Light, QColor(180, 180, 180))
    # light_palette.setColor(QPalette.Midlight, QColor(200, 200, 200))
    # light_palette.setColor(QPalette.Dark, QColor(225, 225, 225))
    light_palette.setColor(QPalette.Text, QColor(0, 0, 0))
    light_palette.setColor(QPalette.BrightText, QColor(0, 0, 0))
    light_palette.setColor(QPalette.ButtonText, QColor(0, 0, 0))
    light_palette.setColor(QPalette.Base, QColor(237, 237, 237))
    light_palette.setColor(QPalette.Window, QColor(240, 240, 240))
    # light_palette.setColor(QPalette.Shadow, QColor(20, 20, 20))
    light_palette.setColor(QPalette.Highlight, QColor(76, 163, 224))
    light_palette.setColor(QPalette.HighlightedText, QColor(0, 0, 0))
    light_palette.setColor(QPalette.Link, QColor(0, 162, 232))
    light_palette.setColor(QPalette.AlternateBase, QColor(225, 225, 225))
    light_palette.setColor(QPalette.ToolTipBase, QColor(240, 240, 240))
    light_palette.setColor(QPalette.ToolTipText, QColor(0, 0, 0))

    app.setPalette(light_palette)

def setcoldark(self):#自定义设置
    dark_palette = QPalette()  # 调色板
    dark_palette.setColor(QPalette.Window, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.WindowText, Qt.white)
    dark_palette.setColor(QPalette.Base, QColor(25, 25, 25))
    dark_palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.ToolTipBase, Qt.white)
    dark_palette.setColor(QPalette.ToolTipText, Qt.white)
    dark_palette.setColor(QPalette.Text, Qt.white)
    dark_palette.setColor(QPalette.Button, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.ButtonText, Qt.white)
    dark_palette.setColor(QPalette.BrightText, Qt.red)
    dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.HighlightedText, Qt.black)
    app.setPalette(dark_palette)

    # pass
    # stsetcol2()



if __name__ == "__main__":
    app = QApplication(sys.argv)
    # print(sys.argv[0][-13:-1])
    # print(sys.argv)
    # Force the style to be the same on all OSs:
    app.setStyle("Fusion")

    # Now use a palette to switch to dark colors:
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(53, 53, 53))
    palette.setColor(QPalette.WindowText, Qt.white)
    palette.setColor(QPalette.Base, QColor(25, 25, 25))
    palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    palette.setColor(QPalette.ToolTipBase, Qt.white)
    palette.setColor(QPalette.ToolTipText, Qt.white)
    palette.setColor(QPalette.Text, Qt.white)
    palette.setColor(QPalette.Button, QColor(53, 53, 53))
    palette.setColor(QPalette.ButtonText, Qt.white)
    palette.setColor(QPalette.BrightText, Qt.red)
    palette.setColor(QPalette.Link, QColor(42, 130, 218))
    palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    palette.setColor(QPalette.HighlightedText, Qt.black)
    app.setPalette(palette)


    app.setApplicationName("AtomimgSCNU")
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

