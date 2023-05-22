# import pyqtgraph as pg
# from pyqtgraph import GraphicsLayoutWidget
# from Utilities.IO import IOHelper
# from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from Utilities.Helper import settings
# from pathlib import Path
# import numpy as np
# from PIL import Image
# import datetime
# from queue import Queue
# from PyQt5 import QtGui
from PyQt5.QtGui import QFont, QDoubleValidator

class ImgParameters(QWidget):
    # abs_img = pyqtSignal(dict)

    def __init__(self, parent=None):
        super(ImgParameters, self).__init__(parent=parent)
        self.parent = parent

        self.horizontalGroupBox1 = QGroupBox()
        # self.horizontalGroupBox2 = QGroupBox("flurence Parameters")
        layout1 = QGridLayout()
        # layout2 = QHBoxLayout()

        self.setFont(QFont("Roman times", 12))
        self.mag = QLabel('CamX',self)
        self.magValue = QLineEdit()
        self.magValue.setValidator(QDoubleValidator(self.magValue))
        self.magValue.setMaxLength(6)
        self.tof = QLabel('TOF/ms')
        self.TOF = QLineEdit('')
        self.TOF.setValidator(QDoubleValidator(self.TOF))
        self.TOF.setMaxLength(9)
        layout1.addWidget(self.mag,0,0,1,1)
        layout1.addWidget(self.magValue,0,1,1,1)
        layout1.addWidget(self.tof,0,2,1,1)
        layout1.addWidget(self.TOF,0,3,1,1)

        # ToPwrLabel = QLabel('P/mW')
        # self.ToPwr = QLineEdit()
        # self.ToPwr.setValidator(QDoubleValidator(self.ToPwr))
        # self.ToPwr.setMaxLength(6)
        # DetuLabel = QLabel('Det/MHz')
        # self.Detu = QLineEdit()
        # self.Detu.setValidator(QDoubleValidator(self.Detu))
        # self.Detu.setMaxLength(6)
        # DiaLabel = QLabel('d/mm')
        # self.Dia = QLineEdit()
        # self.Dia.setValidator(QDoubleValidator(self.Dia))
        # self.Dia.setMaxLength(6)

        # layout2.addWidget(ToPwrLabel)
        # layout2.addWidget(self.ToPwr)
        # layout2.addWidget(DetuLabel)
        # layout2.addWidget(self.Detu)
        # layout2.addWidget(DiaLabel)
        # layout2.addWidget(self.Dia)

        self.horizontalGroupBox1.setLayout(layout1)
        # self.horizontalGroupBox2.setLayout(layout2)

        self.vertical_layout = QGridLayout()
        # self.vertical_layout.addLayout(layout2,1,0,1,6)
        self.vertical_layout.addWidget(self.horizontalGroupBox1,0,0,1,6)
        self.setLayout(self.vertical_layout)

        self.default_setting()

        # self.magValue.editingFinished.connect(self.change_mag)
        # self.TOF.editingFinished.connect(self.change_TOF)
        # self.Detu.editingFinished.connect(self.change_Detu)
        # self.Dia.editingFinished.connect(self.change_Dia)
        # self.ToPwr.editingFinished.connect(self.change_ToPwr)

        screen = QDesktopWidget().screenGeometry()  # Control window size
        self.setFixedHeight(screen.width()*(9/16)*13/100)
        self.setMinimumWidth(480)

    # def update_image2(self,img_dict):
    #     self.img.setImage(img_dict['img_data'])
    #     self.img_labelt1.setText(img_dict['img_name'])
    #     # self.data = img_dict['img_data']
    #     # self.data_shape = self.data.shape

    # def push2_state(self):
    #     fpath = IOHelper.get_config_setting('DATA_PATH')
    #     fpath = Path(fpath)
    #     dir_path = fpath.joinpath(str(datetime.datetime.now())[2:].split('.')[0].replace(' ', '-').replace(':', '_'))
    #     # print("save images to {}".format(dir_path))
    #     if settings.m_path != []:
    #         dir_path = settings.m_path
    #     if not dir_path.exists():
    #         dir_path.mkdir()
    #     img_data = np.array(self.img.image)
    #     # load image name by path
    #     img_name = (self.img_labelt1.text())[0:20].replace(' ', '~').replace(':', '_').replace('-', '')
    #     img_data = img_data[::-1]
    #     import numpy
    #     numpy.savetxt(r"{}\{}.ndata".format(dir_path, img_name), img_data, fmt='%.2e', delimiter=' ', newline='\n',
    #                   header='', footer='', comments=' ', encoding=None)
    #     print("save images to {}".format(dir_path))


    # def push_state(self):
    #     if settings.absimgDatas[0] != [] and settings.absimgDatas[1] != [] and settings.absimgDatas[2] != []:
    #         withatom = np.zeros((settings.absimgDatas[0].shape[0], settings.absimgDatas[0].shape[1]))
    #         withoutatom = np.zeros((settings.absimgDatas[1].shape[0], settings.absimgDatas[1].shape[1]))
    #         totalmat = np.zeros((settings.absimgDatas[1].shape[0], settings.absimgDatas[1].shape[1]))
    #         settings.absimgDatas[3] = np.zeros((settings.absimgDatas[1].shape[0], settings.absimgDatas[1].shape[1]))
    #
    #         # print('In the calculation')
    #         import warnings
    #         warnings.filterwarnings("ignore")
    #         for ii in range(settings.absimgDatas[1].shape[0]):
    #             for jj in range(settings.absimgDatas[1].shape[1]):
    #                 withatom[ii, jj] = settings.absimgDatas[0][ii, jj] - settings.absimgDatas[2][ii, jj]  ####
    #                 withoutatom[ii, jj] = settings.absimgDatas[1][ii, jj] - settings.absimgDatas[2][ii, jj]  ###
    #
    #                 if withoutatom[ii, jj] != 0:
    #                     totalmat[ii, jj] = withatom[ii, jj] / withoutatom[ii, jj]  ##########
    #                 else:
    #                     totalmat[ii, jj] = 1
    #
    #                 if totalmat[ii, jj] >= 1 or totalmat[ii, jj] <= 0:
    #                     totalmat[ii, jj] = 1
    #
    #                 settings.absimgDatas[3][ii, jj] = -np.log(totalmat[ii, jj])
    #         # print(settings.absimgData[3][0:20,0:20])
    #
    #         timestamp = datetime.datetime.now()
    #         # self.abs_img.emit({'img_name': str(timestamp)[2:], 'img_data': settings.absimgDatas[3]})
    #         self.img_Push2.setEnabled(True)
    #     else:
    #         print('Please add images')


    def default_setting(self):

        self.magValue.setText(str(settings.widget_params["Analyse Data Setting"]["magValue"]))
        self.TOF.setText(str(settings.widget_params["Analyse Data Setting"]["TOF"]))
        # self.Detu.setText(str(settings.widget_params["Analyse Data Setting"]["Detu"]))
        # self.Dia.setText(str(settings.widget_params["Analyse Data Setting"]["Dia"]))
        # self.ToPwr.setText(str(settings.widget_params["Analyse Data Setting"]["ToPwr"]))

    # def change_mag(self):
    #     settings.widget_params["Analyse Data Setting"]["magValue"] = self.magValue.text()
    #
    # def change_TOF(self):
    #     settings.widget_params["Analyse Data Setting"]["TOF"] = self.TOF.text()


# class PlotWindow(QWidget):
#
#     img_dict = pyqtSignal(object)
#
#     myserial = 5
#
#     def __init__(self):
#         super(PlotWindow, self).__init__()
#         self.layout = QHBoxLayout(self)
#
#         pg.setConfigOptions(imageAxisOrder='row-major')
#         self.viewport = GraphicsLayoutWidget()
#         self.video_view = self.viewport.addViewBox()
#         self.video = pg.ImageItem()
#         self.video_view.addItem(self.video)
#         self.video_view.setMouseEnabled(x=False, y=False)#make it can not move
#
#         self.setLayout(self.layout)
#
#         self.layout.addWidget(self.viewport)
#         self.img_label = QLabel()
#
#         # self.horizontalLayout = QVBoxLayout()
#         # self.horizontalLayout.addWidget(self.img_label)
#         # self.layout.addLayout(self.horizontalLayout)
#         screen = QtGui.QDesktopWidget().screenGeometry()
#         self.setFixedSize(screen.width() * 15 / 100, screen.height() * 14.5 / 100)
#
#
#     def mousePressEvent(self, event):
#         if event.button() == Qt.LeftButton:
#             try:
#                 fpath = IOHelper.get_config_setting('DATA_PATH')
#
#                 img_fpath = QFileDialog.getOpenFileName(self, "Open File", fpath)  # name path
#                 strimg_fpath = str(img_fpath)
#                 img_file = strimg_fpath[2:len(strimg_fpath) - 19]
#                 img_path = Path(img_file)
#
#                 file = open(img_path)
#                 linescontent = file.readlines()  # Read the file as a behavior unit
#                 rows = len(linescontent)  # get the numbers fo line
#                 lines = len(linescontent[0].strip().split(' '))
#                 img_data = np.zeros((rows, lines))  # Initialization matrix
#                 row = 0
#                 for line in linescontent:
#                     line = line.strip().split(' ')
#                     img_data[row, :] = line[:]
#                     row += 1
#                 file.close()
#
#                 img_data = img_data[::-1]
#                 img_name = img_path.stem
#                 img = {
#                     'img_name': img_name,
#                     'img_data': img_data}
#                 settings.absimgDatas[self.myserial] = img_data
#
#                 self.img_plot(img)
#             except TypeError:
#                 return
#             except PermissionError:
#                 return
#
#     def update_console(self, stri):
#         MAX_LINES = 50
#         stri = str(stri)
#         new_text = self.prompt_dock.console_text() + '\n' + stri
#         line_list = new_text.splitlines()
#         N_lines = min(MAX_LINES, len(line_list))
#         # limit output lines
#         new_text = '\n'.join(line_list[-N_lines:])
#         self.prompt_dock.console_text(new_text)
#         self.prompt_dock.automatic_scroll()
#
#     def save_image(self):
#         try:
#             if self.video.image is None:
#                 print("have no image in window")
#                 return
#             fpath = IOHelper.get_config_setting('DATA_PATH')
#             fpath = Path(fpath)
#             dir_path = fpath.joinpath(str(datetime.datetime.now()).split('.')[0].replace(' ', '-').replace(':', '_'))
#             # print("save images to {}".format(dir_path))
#             if not dir_path.exists():
#                 dir_path.mkdir()
#                 img_data = np.array(self.video.image)
#                 # load image name by path
#                 img_name1 = settings.widget_params["Analyse Data Setting"]["Prefix"]
#                 img_name2 = (self.img_label.text())[0:20].replace(' ', '~').replace(':', '').replace('-', '')
#                 img_name = str(img_name1) + str(img_name2)
#                 img_data = img_data[::-1]
#                 img_data = Image.fromarray(img_data)
#                 img_data.save(r"{}\{}.png".format(dir_path, img_name))
#             print("save images to {}".format(dir_path))
#             # print("images have saved.")
#         except OSError:
#             print('Only new version files can be saved.')
#
#     def img_plot(self, img_dict):
#         self.video.setImage(img_dict['img_data'])
#         self.img_label.setText(img_dict['img_name'])
#
#     def clear_win(self):
#         self.video.clear()
#         self.img_label.setText('')
