import pyqtgraph as pg
from pyqtgraph import GraphicsLayoutWidget
from Utilities.IO import IOHelper
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from Utilities.Helper import settings
from pathlib import Path
import numpy as np
from numpy import save
from PIL import Image
import datetime
from queue import Queue
from PyQt5 import QtGui
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from pyqtgraph.Qt import *
from Widget.CoreWidget import ParametersWidget

# w = QWidget()
# img_queue = ImgQueueWidget()
# scroll = QScrollArea()
# scroll.setWidget(img_queue)
#
# vbox = QVBoxLayout()
# vbox.addWidget(scroll)
# w.setLayout(vbox)

class ImgQueueWidget(QWidget):

    def __init__(self, parent=None):
        super(ImgQueueWidget, self).__init__(parent)
        # plot image history
        self.verticalLayout = QVBoxLayout()
        self.plot_wins = Queue(settings.widget_params['Image Display Setting']['img_stack_num'])  # 10
        for i in range(settings.widget_params['Image Display Setting']['img_stack_num']):
            # 添加十个小窗口，小窗口类在下面
            plot_win = PlotWindow()
            plot_win.video.image = None
            self.plot_wins.put(plot_win)
            self.verticalLayout.addWidget(plot_win)
            if (i+1) % 4 == 0:
                lab = QLabel('***', self)
                self.verticalLayout.addWidget(lab)

        # self.s1 = QSlider(Qt.Vertical)
        # self.s1.setTickPosition(QSlider.TicksLeft)
        # self.verticalLayout.addWidget(self.s1)

        self.setLayout(self.verticalLayout)
        screen = QtGui.QDesktopWidget().screenGeometry()
        # self.setFixedSize(screen.width()*14/100,screen.width()*(9/16)*(60/100))
        # print(self.width(), self.height())

#小窗口类
class PlotWindow(QWidget):

    img_dict = pyqtSignal(object)

    def __init__(self):
        super(PlotWindow, self).__init__()
        self.layout = QHBoxLayout(self)

        pg.setConfigOptions(imageAxisOrder='row-major')
        self.viewport = GraphicsLayoutWidget()
        self.video_view = self.viewport.addViewBox()
        self.video = pg.ImageItem()
        from matplotlib import cm
        colormapname = cm.get_cmap('jet')  # cm.get_cmap("CMRmap")
        colormapname._init()
        lut = (colormapname._lut * 255).view(np.ndarray)
        self.video.setLookupTable(lut)
        # self.video_view.clicked.connect(self.btn_state)
        self.video_view.addItem(self.video)
        self.video_view.setMouseEnabled(x=False, y=False)#make it can not move
        self.img_label = 'no image'
        self.viewport.setToolTip(str(self.img_label)) #str(self.img_label)
        self.setLayout(self.layout)

        self.layout.addWidget(self.viewport)

        palette = QPalette()  # 调色板
        # palette.setColor(QPalette.Window, QColor(193, 205, 205))
        # palette.setColor(QPalette.WindowText, Qt.black)
        palette.setColor(QPalette.Base, QColor(53, 53, 53))
        # palette.setColor(QPalette.AlternateBase, QColor(25, 25, 25))
        # palette.setColor(QPalette.ToolTipBase, Qt.black)
        # palette.setColor(QPalette.ToolTipText, Qt.black)
        # palette.setColor(QPalette.Text, Qt.black)
        # palette.setColor(QPalette.Button, QColor(25, 25, 25))
        # palette.setColor(QPalette.ButtonText, Qt.black)
        # palette.setColor(QPalette.BrightText, Qt.red)
        # palette.setColor(QPalette.Link, QColor(42, 130, 218))
        # palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        # palette.setColor(QPalette.HighlightedText, Qt.white)
        self.setPalette(palette)

        screen = QDesktopWidget().screenGeometry()
        # print(screen)
        self.setFixedSize(screen.width() * (9/16)*(14 / 100) * (44/33.5), screen.width() * (9/16)*(14 / 100))

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.btn_state()
        self.vLine = pg.InfiniteLine(angle=90, movable=False)
        self.hLine = pg.InfiniteLine(angle=0, movable=False)
        self.vLine.setPen(color='r', width=2, style=QtCore.Qt.DashLine)  # 设置画线的笔
        self.hLine.setPen(color='r', width=2, style=QtCore.Qt.DashLine)
        # self.vLine.setPos(self.roi.pos()[0] + self.roi.size()[0] / 2)
        # self.hLine.setPos(self.roi.pos()[1] + self.roi.size()[1] / 2)

    # def mouseReleaseEvent(self, event):
    #     if event.button() == Qt.LeftButton:
    #         print('r')

    def mouseDoubleClickEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.save_image()

    # def enterEvent(self, event):
    #     pass
    #
    # def leaveEvent(self, event):
    #     pass

    def btn_state(self):
        if self.video.image is None:
            print("have no image in window")
            # from MainWindow import TestMainWindow
            # TestMainWindow.path.setTitle(str('have no image in window'))
            return
        # img_analyse_setting.roi.setChecked(False)
        img_dict = {'img_data': np.array(self.video.image), 'img_name': self.img_label}
        settings.imgData["Img_data"] = img_dict['img_data']
        self.img_dict.emit(img_dict)

    def save_image(self):
        # try:
        if self.video.image is None:
            print("have no image in window")
            return
        fpath = IOHelper.get_config_setting('DATA_PATH')
        fpath = Path(fpath)
        dir_path = fpath
        # dir_path = fpath.joinpath(str(datetime.datetime.now())[2:].split('.')[0].replace(' ', '').replace(':', '_'))
        if settings.m_path != []:
            dir_path = settings.m_path
        # print("save images to {}".format(dir_path))
        if not dir_path.exists():
            dir_path.mkdir()
        img_data = np.array(self.video.image)
        # load image name by path
        img_name2 = (self.img_label)[0:25].replace(' ', '_').replace(':', '').replace('-', '').replace('.', '_')
        img_name = str(img_name2)
        img_data = img_data[::-1]
        # img_data = Image.fromarray(img_data)
        # img_data.save(r"{}\{}.png".format(dir_path, img_name))
        save(r"{}\{}".format(dir_path, img_name+'.DATA'), img_data)
        print("save images to {}".format(dir_path))
            # print("images have saved.")
        # except OSError:
        #     print('Image cannot be saved.')

#############################################################
    def img_plot(self, img_dict):
        self.video.clear()
        self.video.setImage(img_dict['img_data'])
        self.img_label = img_dict['img_name']
        self.viewport.setToolTip(str(self.img_label))


    def clear_win(self):
        self.video.clear()
        self.img_label = 'no image'
        self.viewport.setToolTip(str(self.img_label))


