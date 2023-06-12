import pyqtgraph as pg
# from pyqtgraph import ColorMap
from pyqtgraph.Qt import QtGui,QtCore
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from scipy import optimize
from PyQt5.QtGui import *
from Model.DataAnalysis.CaculateAtoms import *
# from decimal import *
from PyQt5.QtGui import QFont
from copy import deepcopy
from matplotlib import cm
# import matplotlib.pyplot as plt
import numpy as np
# getcontext().prec = 4#Set significant number



class PlotMainWindow(QWidget):

    atom_number = pyqtSignal(object)  #from PyQt5.QtCore
    Pxatom_num = pyqtSignal(object)
    atom_numberF = pyqtSignal(object)
    Pxatom_numF = pyqtSignal(object)
    # TotalPhotons_num = pyqtSignal(object)
    fittingdata = pyqtSignal(dict)
    roipos = pyqtSignal(list)
    tempnum = pyqtSignal(object)
    HOD = pyqtSignal(object)
    VOD = pyqtSignal(object)
    # temcula = pyqtSignal(list)

    def __init__(self):
        super(PlotMainWindow, self).__init__()
        self.layout = QGridLayout(self)

        # win = pg.GraphicsView()
        l = pg.GraphicsLayout(border=(100, 100, 100))   #图形布局
        win = pg.GraphicsLayoutWidget()
        win.setCentralItem(l)
        pg.setConfigOptions(imageAxisOrder='row-major')
        # pg.setConfigOptions(leftButtonPan=False)
        self.viewBox = l.addPlot()
        self.viewBox.hideAxis('left')#hide the left and right
        self.viewBox.hideAxis('bottom')
        self.img = pg.ImageItem()

        colormap = cm.get_cmap("jet")  # cm.get_cmap("CMRmap")
        colormap._init()
        lut = (colormap._lut * 255).view(np.ndarray)  # Convert matplotlib colormap from 0-1 to 0 -255 for Qt
        self.img.setLookupTable(lut)

        self.viewBox.setMouseEnabled(x=False, y=False)#make image can not move
        # pg.setConfigOptions(leftButtonPan=False)
        self.viewBox.addItem(self.img)
        self.layout.addWidget(win,0,0,1,8)
        self.img_label = QLabel()
        self.img_label.setFont(QFont("Roman times", 10))
        self.img_label1 = QLabel('| Temperature/K')
        self.img_label1.setFont(QFont("Roman times", 10))
        self.img_label2 = QLabel('| HpeakOD')
        self.img_label2.setFont(QFont("Roman times", 10))
        self.img_label3 = QLabel()
        self.img_label3.setFont(QFont("Roman times", 10))
        self.img_label4 = QLabel('| VpeakOD')
        self.img_label4.setFont(QFont("Roman times", 10))
        self.img_label5 = QLabel()
        self.img_label5.setFont(QFont("Roman times", 10))
        self.layout.addWidget(self.img_label,1,0,1,2)
        self.layout.addWidget(self.img_label1,1,2,1,2)
        self.layout.addWidget(self.img_label2,1,4,1,1)
        self.layout.addWidget(self.img_label3,1,5)
        self.layout.addWidget(self.img_label4,1,6)
        self.layout.addWidget(self.img_label5,1,7)

        self.setLayout(self.layout)
        self.h_axes = None
        self.v_axes = None
        self.data = None
        self.data_shape = None
        screen = QtGui.QDesktopWidget().screenGeometry()
        self.setFixedSize(screen.width()*44/100,screen.width()*(9/16)*63/100)
        # print(self.width(), self.height())

    def change_HpeakOD(self, HpeakOD):
        self.img_label3.setText(str('%.3e' % HpeakOD))

    def change_VpeakOD(self, VpeakOD):
        self.img_label5.setText(str('%.3e' % VpeakOD))

    def change_TEM(self, TEM):
        self.img_label1.setText(str('| '+'%.3e' % TEM + '/K  '))

    def colormapsel(self, colormapname):
        colormapname = cm.get_cmap(settings.colorlist[colormapname])  # cm.get_cmap("CMRmap")
        colormapname._init()
        lut = (colormapname._lut * 255).view(np.ndarray)     # Convert matplotlib colormap from 0-1 to 0 -255 for Qt
        self.img.setLookupTable(lut)

    def add_roi(self, roi_cbk_state, line_cbk_state):
        if roi_cbk_state.isChecked():
            # video mode doesn't have roi statistics
            if settings.widget_params["Image Display Setting"]["imgSource"] == "camera":
                if settings.widget_params["Image Display Setting"]["mode"] == 0:
                    print("video mode doesn't have roi statistics, please choose another mode.")
                    # 0 doesn't check, 2 means check
                    roi_cbk_state.setCheckState(0)
                    settings.widget_params["Analyse Data Setting"]["roiStatus"] = False
                    return
            if self.data is None:
                print("Main plot window doesn't handle image, please load image first")
                roi_cbk_state.setCheckState(0)
                settings.widget_params["Analyse Data Setting"]["roiStatus"] = False
                return
            roisize = int(settings.widget_params["Analyse Data Setting"]["roisize"])
            self.roi = pg.ROI([300, 300], [roisize, roisize], maxBounds=QtCore.QRect(0, 0, self.data_shape[1], self.data_shape[0]))
            # self.roi.setPen(color=QColor(42, 130, 218), width=3)  # set roi width and color
            self.roi.setPen(color=QColor(225, 225, 225), width=3)  # set roi width and color

            self.viewBox.addItem(self.roi)
            # make sure ROI is drawn above image
            self.roi.setZValue(10)
            self.goto_pos()
            self.calculate_roi()
            # if settings.widget_params["Analyse Data Setting"]["add_ten"]:
            # self.vLine = pg.InfiniteLine(angle=90, movable=False)
            # self.hLine = pg.InfiniteLine(angle=0, movable=False)
            # self.vLine.setPen(color='r', width=3)
            # self.hLine.setPen(color='r', width=3)
            # self.vLine.setPos(self.roi.pos()[0]+self.roi.size()[0]/2)
            # self.hLine.setPos(self.roi.pos()[1]+self.roi.size()[1]/2)
            # self.viewBox.addItem(self.vLine, ignoreBounds=True)
            # self.viewBox.addItem(self.hLine, ignoreBounds=True)
            if settings.widget_params["Analyse Data Setting"]["realtime"] == True:
                self.roi.sigRegionChanged.connect(self.update_ch_fitting_cs)
                self.roi.sigRegionChanged.connect(self.calculate_roi)
            else:
                self.roi.sigRegionChangeFinished.connect(self.update_ch_fitting_cs)
                self.roi.sigRegionChangeFinished.connect(self.calculate_roi)

            self.roi.sigRegionChanged.connect(self.changetenpos)


            settings.widget_params["Analyse Data Setting"]["roiStatus"] = True
        else:
            roi_cbk_state.setCheckState(0)
            # axes_cbk_state.setCheckState(0)
            line_cbk_state.setCheckState(0)
            settings.widget_params["Analyse Data Setting"]["roiStatus"] = False
            settings.widget_params["Analyse Data Setting"]["add_rawdata"] = False
            # remove viewBox's items
            # self.viewBox.clear()
            # add image item
            try:
                self.viewBox.removeItem(self.roi)
            except AttributeError:
                pass
    def setroichange(self):
        try:
            if settings.widget_params["Analyse Data Setting"]["realtime"] == True:
                self.roi.sigRegionChangeFinished.disconnect(self.update_ch_fitting_cs)
                self.roi.sigRegionChangeFinished.disconnect(self.calculate_roi)
                self.roi.sigRegionChanged.connect(self.update_ch_fitting_cs)
                self.roi.sigRegionChanged.connect(self.calculate_roi)
            else:
                self.roi.sigRegionChanged.disconnect(self.update_ch_fitting_cs)
                self.roi.sigRegionChanged.disconnect(self.calculate_roi)
                self.roi.sigRegionChangeFinished.connect(self.update_ch_fitting_cs)
                self.roi.sigRegionChangeFinished.connect(self.calculate_roi)
        except AttributeError:
            pass


    def change_roisize(self):
        try:
            self.roi.setSize(int(settings.widget_params["Analyse Data Setting"]["roisize"]))
        except AttributeError:
            pass

    def goto_pos(self):
        try:
            xp = deepcopy(settings.widget_params["Analyse Data Setting"]["xpos"])
            yp = deepcopy(settings.widget_params["Analyse Data Setting"]["ypos"])
            if xp > int(self.data_shape[1] - self.roi.size()[0]/2):
                xp = int(self.data_shape[1] - self.roi.size()[0]/2)-1
            elif xp < int(self.roi.size()[0]/2):
                xp = int(self.roi.size()[0])/2
            if yp > int(self.data_shape[0] - self.roi.size()[1]/2):
                yp = int(self.data_shape[0] - self.roi.size()[1]/2)-1
            elif yp < int(self.roi.size()[0]/2):
                yp = int(self.roi.size()[1])/2
            self.roi.setPos(float(xp - int(self.roi.size()[0]/2)), float(yp - int(self.roi.size()[1]/2)))
        except TypeError:
            print('add ROI')
        except AttributeError:
            print('add ROI')

    def add_rawdata(self, cbk_state):
        if cbk_state.isChecked():
            if settings.widget_params["Analyse Data Setting"]["roiStatus"]:
                settings.widget_params["Analyse Data Setting"]["add_rawdata"] = True
                # add horizontal axes and vertical axes
                self.h_axes = self.viewBox.plot()
                self.h_axes.setPen(color='y', width=2)#x
                # TODO: vertical axes hasn't finishe
                self.v_axes = self.viewBox.plot()
                self.v_axes.setPen(color='y', width=2)
                settings.widget_params["Analyse Data Setting"]["add_cross_axes"] = True
                self.vLine = pg.InfiniteLine(angle=90, movable=False)
                self.hLine = pg.InfiniteLine(angle=0, movable=False)
                self.vLine.setPen(color='r', width=2, style=QtCore.Qt.DashLine) #设置画线的笔
                self.hLine.setPen(color='r', width=2, style=QtCore.Qt.DashLine)
                self.vLine.setPos(self.roi.pos()[0] + self.roi.size()[0] / 2)
                self.hLine.setPos(self.roi.pos()[1] + self.roi.size()[1] / 2)
                self.viewBox.addItem(self.vLine, ignoreBounds=True)
                self.viewBox.addItem(self.hLine, ignoreBounds=True)
                self.update_ch_fitting_cs()

            else:
                print("please add roi first.")
                # 0 doesn't check, 2 means check
                cbk_state.setCheckState(0)
                settings.widget_params["Analyse Data Setting"]["add_rawdata"] = False
                settings.widget_params["Analyse Data Setting"]["add_cross_axes"] = False
                return
        else:
            cbk_state.setCheckState(0)
            self.viewBox.removeItem(self.h_axes)
            self.viewBox.removeItem(self.v_axes)
            try:
                self.viewBox.removeItem(self.vLine)
                self.viewBox.removeItem(self.hLine)
            except AttributeError:
                pass

            if settings.widget_params["Analyse Data Setting"]["roiStatus"]:
                self.viewBox.addItem(self.roi)
            settings.widget_params["Analyse Data Setting"]["add_rawdata"] = False
            # remove plotItem if cross axes has added
            if self.h_axes is not None and self.v_axes is not None:
                self.viewBox.removeItem(self.h_axes)
                self.viewBox.removeItem(self.v_axes)

    def add_fitting(self, mode):
        # if mode.isChecked():
        if settings.widget_params["Fitting Setting"]["mode"] == 1:
            # pass
            try:
                self.viewBox.removeItem(self.h_axes2)
                self.viewBox.removeItem(self.v_axes2)
            except AttributeError:
                pass
            self.h_axes2 = self.viewBox.plot()
            self.h_axes2.setPen(color='r', width=2)  # x
            self.v_axes2 = self.viewBox.plot()
            self.v_axes2.setPen(color='r', width=2)
            self.update_ch_fitting_cs()
            # print(111)
        else:
            try:
                self.viewBox.removeItem(self.h_axes2)
                self.viewBox.removeItem(self.v_axes2)
            except AttributeError:
                pass
        # else:
        #     try:
        #         self.viewBox.removeItem(self.h_axes2)
        #         self.viewBox.removeItem(self.v_axes2)
        #     except AttributeError:
        #         pass

    # def add_cross_axes(self, cbk_state):
    #     if cbk_state.isChecked():
    #         if settings.widget_params["Analyse Data Setting"]["roiStatus"]:
    #             settings.widget_params["Analyse Data Setting"]["add_cross_axes"] = True
    #             self.vLine = pg.InfiniteLine(angle=90, movable=False)
    #             self.hLine = pg.InfiniteLine(angle=0, movable=False)
    #             self.vLine.setPen(color='r', width=3)
    #             self.hLine.setPen(color='r', width=3)
    #             self.vLine.setPos(self.roi.pos()[0] + self.roi.size()[0] / 2)
    #             self.hLine.setPos(self.roi.pos()[1] + self.roi.size()[1] / 2)
    #             self.viewBox.addItem(self.vLine, ignoreBounds=True)
    #             self.viewBox.addItem(self.hLine, ignoreBounds=True)
    #         else:
    #             settings.widget_params["Analyse Data Setting"]["add_cross_axes"] = False
    #             print("please add roi first.")
    #             cbk_state.setCheckState(0)
    #             return
    #     else:
    #         cbk_state.setCheckState(0)
    #         try:
    #             self.viewBox.removeItem(self.vLine)
    #             self.viewBox.removeItem(self.hLine)
    #         except AttributeError:
    #             pass
    def changetenpos(self):
        try:
            if settings.widget_params["Analyse Data Setting"]["add_cross_axes"]:
                self.vLine.setPos(self.roi.pos()[0] + self.roi.size()[0] / 2)
                self.hLine.setPos(self.roi.pos()[1] + self.roi.size()[1] / 2)
        except AttributeError:
            pass

    def update_ch_fitting_cs(self):
        if settings.cameraON == False:
            try:
                # if settings.widget_params["Analyse Data Setting"]["add_cross_axes"]:
                #     self.vLine.setPos(self.roi.pos()[0]+self.roi.size()[0]/2)
                #     self.hLine.setPos(self.roi.pos()[1]+self.roi.size()[1]/2)
                if settings.widget_params["Analyse Data Setting"]["add_rawdata"] or settings.widget_params["Fitting Setting"]["mode"] == 1:

                    self.datals = self.data
                    h_dataraw = self.datals[int(self.roi.pos()[1] + self.roi.size()[1] / 2), :]
                    h_data = copy.deepcopy(h_dataraw)
                    num_h = range(len(h_data))
                    num_h_data = list(num_h)#水平横坐标范围

                    v_dataraw = self.datals[:, int(self.roi.pos()[0] + self.roi.size()[0] / 2)]
                    v_data = v_dataraw
                    num_v = range(len(v_data))
                    num_v_data = list(num_v)#竖直横坐标范围

                    # vlen = np.ones(len(v_data))# make it at right
                    # vlenlist = list(len(h_data) * vlen)
                    # v_data = list(map(lambda x: x[0] - x[1], zip(vlenlist, v_data)))
                    # print(self.roi.pos())
                    # num_v_data = list([x*len(h_data) for x in num_v_data])
                    if settings.widget_params["Fitting Setting"]["mode"] == 1:
                        # num_h_data2 = num_h_data[int(self.roi.pos()[0]): int(self.roi.pos()[0] + self.roi.size()[0])]
                        num_h_data2 = np.array(num_h_data)
                        import warnings
                        warnings.filterwarnings("ignore")#撤销runtimewarning提醒
                        avalue = (h_data[int(self.roi.pos()[0])] + h_data[int(self.roi.pos()[0])-1] + h_data[int(self.roi.pos()[0])+1])/3
                        p0 = [avalue, int(self.roi.pos()[0] + self.roi.size()[0] / 2), 100, 0]
                        plesq = optimize.leastsq(residuals, p0, args=(h_data, num_h_data2))
                        list1 = [plesq[0][0],plesq[0][1],plesq[0][2],plesq[0][3]]
                        settings.data3 = round(plesq[0][2],2)
                        # list1[3] = 0  #tuple cannot change, so list it
                        # if list1[0] < 0 or list1[0] > 1500 or list1[1] > 2000 or abs(list1[2]) > 100:
                        #     list1 = [0,0,0,0]       #When the fitting data deviates too much, do not draw, Same thing down here.
                        list1 = np.array(tuple(list1))
                        h_data2 = peval(num_h_data, list1)

                        num_v_data2 = np.array(num_v_data)
                        # v_data2 = self.datals[: , int(self.roi.pos()[0] + self.roi.size()[0] / 2)]
                        # v_data2 = self.data[int(self.roi.pos()[1]): int(self.roi.pos()[1] + self.roi.size()[1]), int(self.roi.pos()[0] + self.roi.size()[0] / 2)]
                        # v_data2 = np.array(v_data2)
                        avalue2 = (v_data[int(self.roi.pos()[1])] + v_data[int(self.roi.pos()[1]) - 1] + v_data[int(self.roi.pos()[1]) + 0]) / 3
                        p1 = [avalue2, int(self.roi.pos()[1] + self.roi.size()[1] / 2), 100, 0]
                        plesq2 = optimize.leastsq(residuals, p1, args=(v_data, num_v_data2))
                        list2 = [plesq2[0][0], plesq2[0][1], plesq2[0][2], plesq2[0][3]]
                        settings.data6 = round(plesq2[0][2],2)
                        # if list2[2] > 1500:
                        # list2[3] = 0
                        # if list2[0] < 0 or list2[0] > 1500 or list2[1] > 2000 or abs(list2[2]) > 100:
                        #     list2 = [0, 0, 0, 0]   #

                        list2 = np.array(tuple(list2))
                        v_data2 = peval(num_v_data, list2)
                        data7 = abs(2*np.pi*np.sqrt(plesq[0][0]*plesq2[0][0])*plesq[0][2]*plesq2[0][2])
                        CCDPlSize = [3.75, 3.75]
                        Magnification = float(settings.widget_params["Analyse Data Setting"]["magValue"])
                        data8 = plesq[0][2]*CCDPlSize[0]*Magnification*1E-3
                        data9 = plesq2[0][2]*CCDPlSize[1]*Magnification*1E-3
                        self.fittingdata.emit({'data1': plesq[0][0] , 'data2': plesq[0][1], 'data3': plesq[0][2], 'data4': plesq2[0][0],'data5': plesq2[0][1], 'data6': plesq2[0][2],'data7': data7, 'data8': data8,'data9': data9})
                        # h_data2[h_data2 < 0] = 0
                        # v_data2[v_data2 < 0] = 0
                        self.h_axes2.setData(num_h_data2, self.Modify_data_size(h_data2, self.datals.shape[0]))
                        self.v_axes2.setData(self.Modify_data_size(v_data2, self.datals.shape[1]), num_v_data2)

                    if settings.widget_params["Analyse Data Setting"]["add_rawdata"]:
                        # h_data[h_data < 0] = 0
                        # v_data[v_data < 0] = 0
                        self.h_axes.setData(num_h_data, self.Modify_data_size(h_data, self.datals.shape[0]))
                        self.v_axes.setData(self.Modify_data_size(v_data, self.datals.shape[1]), num_v_data)
            except AttributeError:
                pass
        else:
            print('It is recommended to turn off the camera first')

    def Modify_data_size(self, data, PictureMax):
        #修改数据的大小，确保数据不超出图像范围
        Modify = max(data) / (PictureMax*0.2)

        return data/Modify

    def calculate_roi(self):
        # calculate atom number
        try:
            if self.roi.pos()[0] < 0 or self.roi.pos()[1] < 0 or self.roi.size()[1] > self.data_shape[1] or self.roi.size()[0] > self.data_shape[0]:
                return
            
            #print('before TotalPhotons is ',(sum(self.data[int(self.roi.pos()[1]):int(self.roi.pos()[1] + self.roi.size()[0]), 
            #                                 int(self.roi.pos()[0]):int(self.roi.pos()[0] + self.roi.size()[1])])))

            TotalPhotons = sum(sum(self.data[int(self.roi.pos()[1]):int(self.roi.pos()[1] + self.roi.size()[0]), 
                                             int(self.roi.pos()[0]):int(self.roi.pos()[0] + self.roi.size()[1])]))
            #modi
            print('TotalPhotons is ',TotalPhotons)
            ROIsize = self.roi.size()[0]*self.roi.size()[1]

            #以吸收成像的参数计算原子数
            calculatedata = calculateAtom(TotalPhotons, ROIsize, 1)
            self.atom_number.emit(calculatedata[0])
            self.Pxatom_num.emit(calculatedata[1])

            #以荧光成像的参数计算原子数
            calculatedataF = calculateAtom(TotalPhotons, ROIsize, 0)
            self.atom_numberF.emit(calculatedataF[0])
            self.Pxatom_numF.emit(calculatedataF[1])

            xpos = self.roi.pos()[1] + self.roi.size()[1]/2
            ypos = self.roi.pos()[0] + self.roi.size()[0]/2
            position = [ypos, xpos]
            self.roipos.emit(position)
        except:
            pass

    def judrcf(self, mode1, mode2, mode3):
        if settings.widget_params["Analyse Data Setting"]["roiStatus"] == True:

            self.add_roi(mode1, mode2)
            if settings.widget_params["Analyse Data Setting"]["add_rawdata"] == True:
                self.add_rawdata(mode2)
            if settings.widget_params["Fitting Setting"]["mode"] == 1:
                # print(1)
                self.add_fitting(mode3)

    def measureT(self):
        h_data = self.data[int(settings.ceny), :]
        num_h = range(len(h_data))
        num_h_data = list(num_h)
        num_h_data2 = np.array(num_h_data)
        h_data2 = self.data[int(settings.ceny), :]
        h_data2 = np.array(h_data2)
        import warnings
        warnings.filterwarnings("ignore")  # 撤销runtimewarning提醒
        avalue01 = (h_data[int(settings.cenx)] + h_data[int(settings.cenx) - 1] + h_data[int(settings.cenx) + 1]) / 3
        # p0 = [avalue, int(self.roi.pos()[0] + self.roi.size()[0] / 2), 100, 0]
        p00 = [avalue01, int(settings.cenx), int(200 / 2), 0]
        plesq3 = optimize.leastsq(residuals, p00, args=(h_data2, num_h_data2))
        h = [0, abs(round(plesq3[0][2], 2)) ** 2]
        list3 = [plesq3[0][0], plesq3[0][1], plesq3[0][2], plesq3[0][3]]
        list3 = np.array(tuple(list3))
        h_peakOD = max(peval(num_h_data2, list3))
        self.HOD.emit(h_peakOD)

        v_data = self.data[:, int(settings.cenx)]
        num_v = range(len(v_data))
        num_v_data = list(num_v)
        num_v_data2 = np.array(num_v_data)
        v_data2 = self.data[:, int(settings.cenx)]
        v_data2 = np.array(v_data2)
        avalue02 = (h_data[int(settings.ceny)] + h_data[int(settings.ceny) - 1] + h_data[int(settings.ceny) + 1]) / 3
        p11 = [avalue02, int(settings.ceny), int(200 / 2), 0]
        plesq4 = optimize.leastsq(residuals, p11, args=(v_data2, num_v_data2))
        v = [0, abs(round(plesq4[0][2], 2)) ** 2]
        list4 = [plesq4[0][0], plesq4[0][1], plesq4[0][2], plesq4[0][3]]
        list4 = np.array(tuple(list4))
        v_peakOD = max(peval(num_v_data2, list4))
        self.VOD.emit(v_peakOD)

        t = [0,(float(settings.widget_params["Analyse Data Setting"]["TOF"]) / 1000) ** 2]

        temp_fit = np.polyfit(t, h, 1)
        temp_fit2 = np.polyfit(t, v, 1)
        mass = 85.4678 * 1.66053886 * 10 ** (-27)
        kbol = 1.38 * 1e-23
        CCDPlSize = [3.75, 3.75]
        pixelArea = CCDPlSize[0] * CCDPlSize[1] * 1e-12  # μm**2
        M = float(settings.widget_params["Analyse Data Setting"]["magValue"])
        # tem = (temp_fit[0] ** 2 * M ** 2 * mass / kbol / pixelArea + temp_fit2[0] ** 2 * M ** 2 * mass / kbol / pixelArea) / 2
        tem = (temp_fit[0] * mass / M**2 / kbol * pixelArea + temp_fit2[0] * mass / M**2 / kbol * pixelArea) / 2

        # self.tempnum.connect(self.result_dock.change_temp)
        self.tempnum.emit(tem)

    def img_plot(self, img_dict):
        """
        design for software mode and hardware mode, choose image from image stack to display in main window
        :param img_dict:
        :return:
        """
        self.viewBox.clear()
        self.viewBox.addItem(self.img)
        self.img.clear()
        self.img.setImage(img_dict['img_data'])
        # print(img_dict['img_data'][794,420:450])
        self.data = img_dict['img_data']
        print('self.data_shape is',self.data_shape)
        print(img_dict) #改了还未运行
        # print('self.data is',self.data)


        import cv2
        # gray_image = cv2.cvtColor(img_dict['img_data'], cv2.COLOR_BGR2GRAY)
        # ret, thresh = cv2.threshold(img_dict['img_data'], 0, 255, 3)
        # image, contours, hierarchy = cv2.findContours(thresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE) 旧版本返回三个参数，新版本返回2个
        # contours, hierarchy = cv2.findContours(img_dict['img_data'], cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

        # cnt = contours[0]
        m = cv2.moments(img_dict['img_data'])
        #modi
        # plt.imshow(m)
        # plt.show()
        # print(m)
        x = int(m['m10'] / m['m00'])
        y = int(m['m01'] / m['m00'])
        mp = cv2.moments(img_dict['img_data'][y-200:y+200, x-200:x+200])
        xp = int(mp['m10'] / mp['m00'])
        yp = int(mp['m01'] / mp['m00'])
        mp2 = cv2.moments(img_dict['img_data'][yp+y-200-200:yp+y-200+200, xp+x-200-200:xp+x-200+200])
        xp2 = int(mp2['m10'] / mp2['m00'])
        yp2 = int(mp2['m01'] / mp2['m00'])
        mp3 = cv2.moments(img_dict['img_data'][yp+y-200+yp2-200-200:yp+y-200+yp2-200+200, xp+x-200+xp2-200-200:xp+x-200+xp2-200+200])
        xp3 = mp3['m10'] / mp3['m00']
        yp3 = mp3['m01'] / mp3['m00']
        # print('x=', xp+x-200, 'y=', yp+y-200)
        settings.widget_params["Analyse Data Setting"]["xpos"] = x+xp-200+xp2-200+xp3-200
        settings.widget_params["Analyse Data Setting"]["ypos"] = y+yp-200+yp2-200+yp3-200
        settings.cenx = x+xp-200+xp2-200+xp3-200
        settings.ceny = y+yp-200+yp2-200+yp3-200

        self.data_shape = self.data.shape
        self.img_label.setText(img_dict['img_name'])
        #Set the initial axis so that the size of the image remains the same when the axis is added
        numh = range(self.data_shape[1])
        numh = list(numh)
        numhd = np.zeros(self.data_shape[1])
        numv = range(self.data_shape[0])
        numv = list(numv)
        numvd = np.zeros(self.data_shape[0])
        self.h_axesm = self.viewBox.plot()
        self.h_axesm.setPen(color='k', width=2)  # x
        self.h_axesm.setData(numh, numhd)
        self.v_axesm = self.viewBox.plot()
        self.v_axesm.setPen(color='k', width=2)
        self.v_axesm.setData(numvd, numv)
        # if settings.widget_params["Analyse Data Setting"]["roiStatus"] == True:
            # roi_cbk_state.setCheckState(1)
        if settings.cameraON == False:
            self.measureT()


    def img_plot2(self):
        if settings.imgData["BkgImg"] !=[] and settings.imgData["Img_data"] !=[]:
            self.viewBox.clear()
            self.viewBox.addItem(self.img)
            self.img.clear()
            # Reloading improves operational efficiency
            subt = deepcopy(settings.imgData["Img_data"])
            # print(settings.imgData["Img_data"][0, 0:10])
            # print(settings.imgData["BkgImg"][0, 0:10])
            settings.imgData["Img_data"] = settings.imgData["Img_data"] - settings.imgData["BkgImg"]
            settings.imgData["Img_data"][settings.imgData["Img_data"] > subt] = 0
            self.img.setImage(settings.imgData["Img_data"])
            self.data = settings.imgData["Img_data"]
            self.data_shape = settings.imgData["Img_data"].shape
            self.update_ch_fitting_cs()
            try:
                self.calculate_roi()
            except AttributeError:
                pass
        elif settings.imgData["Img_data"] ==[]:
            print('No image')
        elif settings.imgData["BkgImg"] ==[]:
            print('No background image')


    def img_plot3(self):
        if settings.imgData["Img_data"] != []:
            self.viewBox.clear()
            self.viewBox.addItem(self.img)
            self.img.clear()
            # Reloading improves operational efficiency
            settings.imgData["Img_photon_range"] = deepcopy(settings.imgData["Img_data"])
            try:
                settings.imgData["Img_photon_range"][settings.imgData["Img_data"] <= float(settings.widget_params["Image Display Setting"]["pfMin"])] = settings.widget_params["Image Display Setting"]["pfMin"]
                print('imgDatamax',max(settings.imgData))
                print('Img_photon_range',settings.imgData["Img_photon_range"])
                settings.imgData["Img_photon_range"][settings.imgData["Img_data"] >= float(settings.widget_params["Image Display Setting"]["pfMax"])] = settings.widget_params["Image Display Setting"]["pfMax"]
            except ValueError:
                print('The edit box cannot be empty')
            self.img.setImage(settings.imgData["Img_photon_range"])
            self.data = settings.imgData["Img_photon_range"]
            self.data_shape = settings.imgData["Img_photon_range"].shape
            self.update_ch_fitting_cs()
            try:
                self.calculate_roi()
            except AttributeError:
                pass
            # print('photon filter ： finish.')
        else:
            print('No image')


    def clear_win(self):
        self.viewBox.clear()
        # add image item
        self.viewBox.addItem(self.img)
        if self.img is None:
            return
        self.img.clear()
        self.img_label.setText('')
        self.img_label1.setText('| Tem/K')
        self.img_label2.setText('| HpeakOD')
        self.img_label3.setText('')
        self.img_label4.setText('| VpeakOD')
        self.img_label5.setText('')
        self.data = None
        self.data_shape = None

    def ClearPlotwin(self):
        self.viewBox.clear()
        # add image item
        self.viewBox.addItem(self.img)
        self.img.clear()
        self.img_label.setText('Waiting for the trigger')
        self.data = None
        self.data_shape = None

def func(xx, aa, bb, cc, dd):
    return aa * np.e ** (-((xx - bb) ** 2) / 2 / cc ** 2) + dd

# def logistic4(x, A, B, C, D):
#     return (A-D)/(1+(x/C)**B)+D

def residuals(p, y, x):
    [A, B, C, D] = p
    return y - func(x, A, B, C, D)

def peval(x, p):
    [A, B, C, D] = p
    return func(x, A, B, C, D)
