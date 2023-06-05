import sys
sys.path.append(".")
# import PyCapture2
from multiprocessing import Process
from MainWindow import start_main_win
from Utilities.IO.IOHelper import create_config_file, create_configt_file, create_camconfig_file
#modi
# global StackNum
# StackNum = 0
# global StorageNum
# StorageNum = 0

if __name__ == '__main__':
    create_config_file()
    create_configt_file()
    create_camconfig_file()

    p = Process(target=start_main_win)
    p.start()
    p.join()
    #start_main_win()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def print_camera_info(cam):
    cam_info = cam.getCameraInfo()
    print('\n*** CAMERA INFORMATION ***\n')
    print('Serial number - %d' % cam_info.serialNumber)
    print('Camera model - %s' % cam_info.modelName)
    print('Camera vendor - %s' % cam_info.vendorName)
    print('Sensor - %s' % cam_info.sensorInfo)
    print('Resolution - %s' % cam_info.sensorResolution)
    print('Firmware version - %s' % cam_info.firmwareVersion)
    print('Firmware build time - %s' % cam_info.firmwareBuildTime)
    print()




####################################
# 增加video模式
# 堆栈数量改成4√ 每次在大窗口更新最新图片（OD）
# 修改导入图片路径——看一下为什么手动导入图片无法正常导入
# 原子数 字号调大 （先不搞）
# 十字表 像素量坐标轴显示 进阶：显示到图像外围
# filter 输入 整数改小数
# 修改命名方式
# 修改检测相机函数，把开始循环遍历用的回调参数修改，用成不用开机使用的参数，若不行则测试是否能用不开机使用参数
# 上面解决无法在无相机时正常启动
