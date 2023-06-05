from PyQt5.QtCore import QFileSystemWatcher
from Utilities.Helper import settings, Helper
from Utilities.IO import IOHelper

from Model.DataAnalysis.CaculateAtoms import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtGui import QIntValidator
def file_load_imgs(self):
        """
        Load previous image to stack.
        :return:
        """
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
                
class AutoLoadImg:
    def __init__(self):
        self.watcher = QFileSystemWatcher()
        self.watcher.directoryChanged.connect(self.load_latest_img)

    def start_watching(self, dir_path):
        self.dir_path = dir_path
        self.watcher.addPath(dir_path)

    def load_latest_img(self):
        self.file_load_imgs()


