import sys
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import *
import pyqtgraph as pg

class DemoApp(QMainWindow):                 #创建一个主窗口类
    def __init__(self):
        super(DemoApp, self).__init__()
        self.resize(500, 500)               # 设置窗口的大小
        self.setWindowTitle('DemoApp')      # 设置窗口标题

        #创建一个用于显示图像的窗口部件win
        pglayout = pg.GraphicsLayout(border=(100, 100, 100))
        win = pg.GraphicsLayoutWidget()
        win.setCentralItem(pglayout)
        self.viewBox = pglayout.addPlot()
        self.viewBox.hideAxis('left')        # 隐藏左坐标轴
        self.viewBox.hideAxis('bottom')      # 隐藏下坐标轴
        self.img = pg.ImageItem()
        self.viewBox.addItem(self.img)

        self.setCentralWidget(win)           # 把win设为中心窗口

        DemoWidget_one = DemoWidget()                           # 实例化DemoWidget
        DockWidget = QDockWidget("DemoWidget", self)            # 创建一个停靠窗口
        DockWidget.setWidget(DemoWidget_one)                    # 把DemoWidget放入停靠窗口
        DockWidget.setAllowedAreas(Qt.BottomDockWidgetArea)     # 设置允许停靠位置
        self.addDockWidget(Qt.BottomDockWidgetArea, DockWidget) # 把窗口停靠在中心窗口下方

        self.move_center()         # 把窗口置于屏幕中央

    def move_center(self):
        screen = QDesktopWidget().screenGeometry()
        form = self.geometry()
        x_move_step = (screen.width() - form.width()) / 2
        y_move_step = (screen.height() - form.height()) / 2
        self.move(x_move_step, y_move_step)

class DemoWidget(QWidget):                  #创建一个小窗口类

    def __init__(self):
        super(DemoWidget, self).__init__()

        self.layout = QGridLayout(self)              #创建一个网格布局控件

        self.Label_one = QLabel('This is a label', self)                  #创建一个标签
        self.LineEdit_one = QLineEdit('This is a LineEdit', self)         #创建一个单行文本框
        self.PushBotton_one = QPushButton('This is a PushButton', self)   #创建一个按钮
        self.RadioButton_one = QRadioButton('This is a RadioButton', self)#创建一个单选按钮
        self.CheckBox_one = QCheckBox('This is a CheckBox', self)         #创建一个复选按钮
        self.ComboBox_one = QComboBox(self)                               #创建一个下拉条
        self.Slider_one = QSlider(self)                                   #创建一个滑动条

        list = ['Combo_one', 'Combo_two', 'Combo_three']               #创建一个列表
        self.ComboBox_one.addItems(list)                               #把列表设为下拉条的项目
        self.PushBotton_one.clicked.connect(lambda: self.ActionFunction())  #为按钮连接响应函数

        self.layout.addWidget(self.Label_one,0,0,1,1)           #把各控件置于网格布局的某各位置
        self.layout.addWidget(self.LineEdit_one,0,1)            #前两位数字代表行列位置
        self.layout.addWidget(self.PushBotton_one,1,0)          #后两位数字代表所占位置大小，默认为1*1
        self.layout.addWidget(self.RadioButton_one,1,1)
        self.layout.addWidget(self.CheckBox_one,2,0)
        self.layout.addWidget(self.ComboBox_one,2,1)
        self.layout.addWidget(self.Slider_one,0,2,3,1)          #把滑动条所占位置设为3行1列

    def ActionFunction(self):                              #创建按钮的响应函数
        print('The button worked')

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = DemoApp()
    window.show()
    sys.exit(app.exec_())


################


from Utilities.Helper import Helper

class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.fileMenu = self.menuBar().addMenu("File")   #在主窗口中创建一个菜单

        #创建一个功能操作
        self.LoadImgAction = Helper.create_action(self,"Load Images",
                                                slot=self.file_load_imgs,
                                                shortcut=None,
                                                icon=None,
                                                tip="Load previous images to image stack")

        self.fileMenu.addAction(self.LoadImgAction)    #在该菜单中添加以上操作


from PyQt5.QtGui import QIcon
from Utilities.IO import IOHelper

#基于QAction定义一个创建功能操作的函数
def create_action(parent, text, slot=None, shortcut=None, icon=None, tip=None, checkable=False):
    action = QAction(text, parent)
    if icon is not None:
        icpath = IOHelper.get_configt_setting('DATA_PATH')
        # icpath = Path(icpath)
        icpath = str(icpath).replace('\\','/')
        action.setIcon(QIcon(icpath + "/%s.png" % icon))
    if shortcut is not None:
        action.setShortcut(shortcut)
    if tip is not None:
        action.setToolTip(tip)
        action.setStatusTip(tip)
    if slot is not None:
        action.triggered.connect(slot)
    if checkable:
        action.setCheckable(True)
    return action

########################