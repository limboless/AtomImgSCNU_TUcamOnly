import PyQt5.QtWidgets as QtWidgets
from PyQt5 import QtGui
from PyQt5.QtWidgets import *
from Utilities.Helper import settings
from Utilities.IO import IOHelper
from pathlib import Path
from PyQt5.QtCore import *
import warnings
import numpy as np
import cv2
from scipy import optimize
from PyQt5.QtGui import QFont, QDoubleValidator, QIntValidator



class Absorption(QWidget):

    """
        这是一个窗口部件类，用于显示温度拟合的窗口控件，控件的响应函数计算原子团温度
    """

    MoveTEM = pyqtSignal(object)


    def __init__(self, parent=None):
        super(Absorption, self).__init__(parent)
        # self.setToolTip('Double click to bring up the prompt window')

        # self.typecb = QComboBox()
        # self.typecb.setFont(QFont("Roman times", 13))
        # self.typecb.activated.connect(self.Changet)
        # Process_data = ['Absorbtion Image','Flurence Image']
        # self.typecb.addItems(Process_data)  # get camera num

        self.temp_label = QtWidgets.QLabel('Tempreature')
        self.temp_label.setFont(QFont("Roman times", 18))
        self.temp_num = QtWidgets.QLabel(str('N'))
        self.temp_num.setFont(QFont("Roman times", 24))
        self.hLayout5 = QtWidgets.QHBoxLayout()
        self.hLayout5.addWidget(self.temp_label)
        self.hLayout5.addWidget(self.temp_num)

        self.atom_num_label = QtWidgets.QLabel('Atom')
        self.atom_num_label.setFont(QFont("Roman times", 18))
        self.atom_num = QtWidgets.QLabel(str('N'))
        self.atom_num.setFont(QFont("Roman times", 24))
        self.hLayout1 = QtWidgets.QHBoxLayout()
        self.hLayout1.addWidget(self.atom_num_label)
        self.hLayout1.addWidget(self.atom_num)

        # self.TotalPhotons_num_label = QtWidgets.QLabel('TotalPhotons')
        # self.TotalPhotons_num_label.setFont(QFont("Roman times", 18))
        # self.TotalPhotons_num = QtWidgets.QLabel(str('N'))
        # self.TotalPhotons_num.setFont(QFont("Roman times", 24))
        # self.hLayout3 = QtWidgets.QHBoxLayout()
        # self.hLayout3.addWidget(self.TotalPhotons_num_label)
        # self.hLayout3.addWidget(self.TotalPhotons_num)

        self.atom_numpx_label = QtWidgets.QLabel('Atom/px')
        self.atom_numpx_label.setFont(QFont("Roman times", 18))
        self.atom_numpx = QtWidgets.QLabel(str('N'))
        self.atom_numpx.setFont(QFont("Roman times", 24))
        self.hLayout2 = QtWidgets.QHBoxLayout()
        self.hLayout2.addWidget(self.atom_numpx_label)
        self.hLayout2.addWidget(self.atom_numpx)

        self.layouttemp = QGridLayout()
        self.NUM = QLabel('Number of images')
        self.QTIM = QLabel('Time starting/ms')
        self.JTIM = QLabel('Time interval/ms')
        self.NUME = QLineEdit('5')
        self.QTIME = QLineEdit('1')
        self.JTIME = QLineEdit('1')
        self.NUME.setValidator(QIntValidator(self.NUME))
        self.NUME.setMaxLength(8)
        self.QTIME.setValidator(QDoubleValidator(self.QTIME))
        self.QTIME.setMaxLength(8)
        self.JTIME.setValidator(QDoubleValidator(self.JTIME))
        self.JTIME.setMaxLength(8)
        self.layouttemp.addWidget(self.NUM, 0, 0, 1, 1)
        self.layouttemp.addWidget(self.QTIM, 1, 0, 1, 1)
        self.layouttemp.addWidget(self.JTIM, 2, 0, 1, 1)
        self.layouttemp.addWidget(self.NUME, 0, 1, 1, 1)
        self.layouttemp.addWidget(self.QTIME, 1, 1, 1, 1)
        self.layouttemp.addWidget(self.JTIME, 2, 1, 1, 1)

        self.Tem_pushbutton = QPushButton('Load')
        self.Tem_pushbutton.setFont(QFont("Roman times", 10))
        self.hLayoutP = QtWidgets.QVBoxLayout()
        self.hLayoutP.addWidget(self.Tem_pushbutton)
        self.Tem_pushbutton.clicked.connect(self.TMeasure)


        self.wLayout = QtWidgets.QGridLayout()
        # self.wLayout.addWidget(self.typecb,0,0,1,1)
        self.wLayout.addLayout(self.hLayout5,1,0,1,2)
        # self.wLayout.addLayout(self.hLayout3,2,0,1,2)
        self.wLayout.addLayout(self.hLayout1,3,0,1,2)
        self.wLayout.addLayout(self.hLayout2,4,0,1,2)
        self.wLayout.addLayout(self.hLayoutP,0,0,1,1)
        self.wLayout.addLayout(self.layouttemp,0,1,1,1)
        # self.wLayout.addLayout(self.hLayout6,5,0,1,2)
        # self.wLayout.addLayout(self.hLayout7,6,0,1,2)

        # self.totalLayout = QtWidgets.QHBoxLayout()
        # self.totalLayout.addLayout(self.fLayout)
        # self.totalLayout.addLayout(self.wLayout)
        # self.totalLayout.setStretchFactor(self.fLayout, 1)
        # self.totalLayout.setStretchFactor(self.wLayout, 4)
        self.setLayout(self.wLayout)

        self.prompt = QDialog()  # create a dialog
        self.consoleTextEdit = QtWidgets.QTextEdit()
        self.consoleTextEdit.setReadOnly(True)

        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.addWidget(self.consoleTextEdit)
        self.prompt.setLayout(self.verticalLayout)

        self.MoveTEM.connect(self.change_temp)


        screen = QDesktopWidget().screenGeometry()
        self.setFixedHeight(screen.width()*(9/16)*40/100)
        # self.prompt.setFixedWidth(screen.width()*50/100)

    def TMeasure(self):
        settings.widget_params["Image Display Setting"]["imgSource"] = "disk"
        self.fpath = IOHelper.get_config_setting('DATA_PATH')
        img_fpath = QFileDialog.getOpenFileNames(self, "Open File", self.fpath, "Image files(*.npy *.data)")  # name path
        sigmaTotal_H = []
        sigmaTotal_V = []
        Timeof = []
        try:
            for i in range(int(self.NUME.text())):
                Convertingunits = (float(self.QTIME.text()) + i * float(self.JTIME.text())) / 1000
                Timeof.append((Convertingunits) ** 2)

            if img_fpath[0] != []:
                if len(img_fpath[0]) == int(self.NUME.text()):
                    for i in range(len(img_fpath[0])):
                        try:
                            img_file = str(img_fpath[0][i])
                            img_paths = Path(img_file)
                            sigmasquare = self.Getsigma(img_paths)
                            sigmaTotal_H.append(sigmasquare[0])
                            sigmaTotal_V.append(sigmasquare[1])
                        except TypeError:
                            return
                        except PermissionError:
                            return
                    sigmaTotal_H = sorted(sigmaTotal_H)
                    sigmaTotal_V = sorted(sigmaTotal_V)
                    TEMP = self.GetTem(sigmaTotal_H, sigmaTotal_V, Timeof)
                    self.MoveTEM.emit(TEMP)
                else:
                    print('Abnormal number of images')
        except:
            print('Please make sure that your input is correct')

    def Getsigma(self, img_path):
        self.pathjud = str(img_path)
        self.pathjud = self.pathjud[len(self.pathjud) - 3:]   #Get the version of the file

        if self.pathjud == 'ata':
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

        m = cv2.moments(img_data)

        x = int(m['m10'] / m['m00'])
        y = int(m['m01'] / m['m00'])
        mp = cv2.moments(img_data[y - 200:y + 200, x - 200:x + 200])
        xp = int(mp['m10'] / mp['m00'])
        yp = int(mp['m01'] / mp['m00'])
        mp2 = cv2.moments(img_data[yp + y - 200 - 200:yp + y - 200 + 200, xp + x - 200 - 200:xp + x - 200 + 200])
        xp2 = int(mp2['m10'] / mp2['m00'])
        yp2 = int(mp2['m01'] / mp2['m00'])
        mp3 = cv2.moments(img_data[yp + y - 200 + yp2 - 200 - 200:yp + y - 200 + yp2 - 200 + 200,
                          xp + x - 200 + xp2 - 200 - 200:xp + x - 200 + xp2 - 200 + 200])
        xp3 = mp3['m10'] / mp3['m00']
        yp3 = mp3['m01'] / mp3['m00']

        cenx = x + xp - 200 + xp2 - 200 + xp3 - 200
        ceny = y + yp - 200 + yp2 - 200 + yp3 - 200

        h_data = img_data[int(ceny), :]
        num_h = range(len(h_data))
        num_h_data = list(num_h)
        num_h_data2 = np.array(num_h_data)
        h_data2 = np.array(h_data)

        warnings.filterwarnings("ignore")  # 撤销runtimewarning提醒
        avalue01 = (h_data[int(cenx)] + h_data[int(cenx) - 1] + h_data[int(cenx) + 1]) / 3
        p00 = [avalue01, int(cenx), int(200 / 2), 0]
        plesq3 = optimize.leastsq(residuals, p00, args=(h_data2, num_h_data2))
        h = abs(round(plesq3[0][2], 2)) ** 2

        v_data = img_data[:, int(cenx)]
        num_v = range(len(v_data))
        num_v_data = list(num_v)
        num_v_data2 = np.array(num_v_data)
        v_data2 = np.array(v_data)
        avalue02 = (h_data[int(ceny)] + h_data[int(ceny) - 1] + h_data[int(ceny) + 1]) / 3
        p11 = [avalue02, int(ceny), int(200 / 2), 0]
        plesq4 = optimize.leastsq(residuals, p11, args=(v_data2, num_v_data2))
        v = abs(round(plesq4[0][2], 2)) ** 2

        sigmasquare = [h, v]

        return sigmasquare

    def GetTem(self, sigmasquareH, sigmasquareV, Timeof):
        temp_fit = np.polyfit(Timeof, sigmasquareH, 1)
        mass = 85.4678 * 1.66053886 * 10 ** (-27)
        kbol = 1.38 * 1e-23
        CCDPlSize = [3.75, 3.75]
        pixelArea = CCDPlSize[0] * CCDPlSize[1] * 1e-12  # μm**2
        M = float(settings.widget_params["Analyse Data Setting"]["magValue"])
        TemperatureH = temp_fit[0] * mass / M ** 2 / kbol * pixelArea

        temp_fitV = np.polyfit(Timeof, sigmasquareV, 1)
        mass = 85.4678 * 1.66053886 * 10 ** (-27)
        kbol = 1.38 * 1e-23
        CCDPlSize = [3.75, 3.75]
        pixelArea = CCDPlSize[0] * CCDPlSize[1] * 1e-12  # μm**2
        M = float(settings.widget_params["Analyse Data Setting"]["magValue"])
        TemperatureV = temp_fitV[0] * mass / M ** 2 / kbol * pixelArea

        Temperature = (TemperatureH + TemperatureV) / 2

        return Temperature


    # def Changet(self):
    #     if self.typecb.currentText() == 'Flurence Image':
    #         settings.widget_params["calculate Setting"]["mode"] = 0
    #     elif self.typecb.currentText() == 'Absorbtion Image':
    #         settings.widget_params["calculate Setting"]["mode"] = 1

    # def promptset(self):
    #     self.prompt.setWindowTitle('prompt')
    #     self.prompt.setWindowModality(Qt.NonModal)
    #     self.prompt.setWindowFlags(Qt.WindowStaysOnTopHint)
    #     self.prompt.show()

    # def mouseDoubleClickEvent(self, event):
    #     if event.button() == Qt.LeftButton:
    #         self.promptset()

    # def console_text(self, new_text=None):
    #
    #     """get/set method for the text in the console"""
    #
    #     if new_text == None:
    #
    #         return str((self.consoleTextEdit.toPlainText())).rstrip()
    #
    #     else:
    #
    #         self.consoleTextEdit.setPlainText(new_text)
    #
    # def automatic_scroll(self):
    #     """
    #     performs an automatic scroll up
    #     the latest text shall always be in view
    #     """
    #     sb = self.consoleTextEdit.verticalScrollBar()
    #     sb.setValue(sb.maximum())

    def change_temp(self, temp_num):
        self.temp_num.setText(str('%.3e' % temp_num))

    def change_atom_num(self, atom_num):
        self.atom_num.setText(str('%.3e' % atom_num))


    # def change_TotalPhotons_num(self, TotalPhotons_num):
    #     self.TotalPhotons_num.setText(str('%.3e' % TotalPhotons_num))

    def change_Pxatom_num(self, Pxatom_num):
        self.atom_numpx.setText(str('%.3e' % Pxatom_num))

class Flurence(QWidget):

    def __init__(self, parent=None):
        super(Flurence, self).__init__(parent=parent)
        self.parent = parent
        self.layout2 = QGridLayout()
        # self.setFont(QFont("Roman times", 12))

        self.ToPwrLabel = QLabel('P/mW')
        self.ToPwr = QLineEdit()
        self.ToPwr.setValidator(QDoubleValidator(self.ToPwr))
        self.ToPwr.setMaxLength(8)
        self.DetuLabel = QLabel('Det/MHz')
        self.Detu = QLineEdit()
        self.Detu.setValidator(QDoubleValidator(self.Detu))
        self.Detu.setMaxLength(8)
        self.DiaLabel = QLabel('d/mm')
        self.Dia = QLineEdit()
        self.Dia.setValidator(QDoubleValidator(self.Dia))
        self.Dia.setMaxLength(8)

        self.layout2.addWidget(self.ToPwrLabel,0,0,1,1)
        self.layout2.addWidget(self.ToPwr,     0,1,1,1)
        self.layout2.addWidget(self.DetuLabel, 1,0,1,1)
        self.layout2.addWidget(self.Detu,      1,1,1,1)
        self.layout2.addWidget(self.DiaLabel,  2,0,1,1)
        self.layout2.addWidget(self.Dia,       2,1,1,1)

        self.temp_label = QtWidgets.QLabel('Tempreature')
        self.temp_label.setFont(QFont("Roman times", 18))
        self.temp_num = QtWidgets.QLabel(str('N'))
        self.temp_num.setFont(QFont("Roman times", 24))

        self.atom_num_label = QtWidgets.QLabel('Atom')
        self.atom_num_label.setFont(QFont("Roman times", 18))
        self.atom_num = QtWidgets.QLabel(str('N'))
        self.atom_num.setFont(QFont("Roman times", 24))

        self.atom_numpx_label = QtWidgets.QLabel('Atom/px')
        self.atom_numpx_label.setFont(QFont("Roman times", 18))
        self.atom_numpx = QtWidgets.QLabel(str('N'))
        self.atom_numpx.setFont(QFont("Roman times", 24))

        self.hLayout1 = QtWidgets.QGridLayout()
        self.hLayout1.addWidget(self.temp_label,0,0,1,1)
        self.hLayout1.addWidget(self.temp_num,0,1,1,1)
        self.hLayout1.addWidget(self.atom_num_label,1,0,1,1)
        self.hLayout1.addWidget(self.atom_num,1,1,1,1)
        self.hLayout1.addWidget(self.atom_numpx_label,2,0,1,1)
        self.hLayout1.addWidget(self.atom_numpx,2,1,1,1)

        self.vertical_layout = QGridLayout()
        self.vertical_layout.addLayout(self.layout2, 0,0,1,1)
        self.vertical_layout.addLayout(self.hLayout1,1,0,1,2)
        self.setLayout(self.vertical_layout)

        self.default_setting()

        # self.Detu.editingFinished.connect(self.change_Detu)
        # self.Dia.editingFinished.connect(self.change_Dia)
        # self.ToPwr.editingFinished.connect(self.change_ToPwr)

        screen = QtGui.QDesktopWidget().screenGeometry()  # Control window size
        self.setFixedHeight(screen.width()*(9/16)*40/100)
        self.setMinimumWidth(480)
        self.DetuLabel.setFixedWidth(screen.width() * (9 / 16) * 10 / 100)
        self.DiaLabel.setFixedWidth(screen.width() * (9 / 16) * 10 / 100)
        self.ToPwrLabel.setFixedWidth(screen.width() * (9 / 16) * 10 / 100)
        self.Detu.setFixedWidth(screen.width()*(9/16)*15/100)
        self.Dia.setFixedWidth(screen.width()*(9/16)*15/100)
        self.ToPwr.setFixedWidth(screen.width()*(9/16)*15/100)

    def default_setting(self):
        self.Detu.setText(str(settings.widget_params["Analyse Data Setting"]["Detu"]))
        self.Dia.setText(str(settings.widget_params["Analyse Data Setting"]["Dia"]))
        self.ToPwr.setText(str(settings.widget_params["Analyse Data Setting"]["ToPwr"]))

    # def change_Detu(self):
    #     settings.widget_params["Analyse Data Setting"]["Detu"] = self.Detu.text()
    #
    # def change_Dia(self):
    #     settings.widget_params["Analyse Data Setting"]["Dia"] = self.Dia.text()
    #
    # def change_ToPwr(self):
    #     settings.widget_params["Analyse Data Setting"]["ToPwr"] = self.ToPwr.text()

    def change_temp(self, temp_num):
        self.temp_num.setText(str('%.3e' % temp_num))

    def change_atom_num(self, atom_num):
        self.atom_num.setText(str('%.3e' % atom_num))

    def change_Pxatom_num(self, Pxatom_num):
        self.atom_numpx.setText(str('%.3e' % Pxatom_num))

# class promptwin(QtWidgets.QWidget):
#
#     def __init__(self, parent=None):
#         super(promptwin, self).__init__(parent)
#         self.consoleTextEditp = QtWidgets.QTextEdit()
#         self.consoleTextEditp.setReadOnly(True)
#
#         self.verticalLayoutp = QtWidgets.QVBoxLayout()
#         self.verticalLayoutp.addWidget(self.consoleTextEditp)
#         self.setLayout(self.verticalLayoutp)
#
#     def console_text(self, new_text=None):
#
#         """get/set method for the text in the console"""
#
#         if new_text == None:
#
#             return str((self.consoleTextEditp.toPlainText())).rstrip()
#
#         else:
#             self.consoleTextEditp.setPlainText(new_text)
#
#     def automatic_scroll(self):
#         """
#         performs an automatic scroll up
#         the latest text shall always be in view
#         """
#         sbp = self.consoleTextEditp.verticalScrollBar()
#         sbp.setValue(sbp.maximum())

def func(xx, aa, bb, cc, dd):
    return aa * np.e ** (-((xx - bb) ** 2) / 2 / cc ** 2) + dd

def residuals(p, y, x):
    [A, B, C, D] = p
    return y - func(x, A, B, C, D)

def peval(x, p):
    [A, B, C, D] = p
    return func(x, A, B, C, D)