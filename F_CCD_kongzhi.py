import PyCapture2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def print_build_info():
    lib_ver = PyCapture2.getLibraryVersion()
    print('PyCapture2 library version: %d %d %d %d' % (lib_ver[0], lib_ver[1], lib_ver[2], lib_ver[3]))
    print()

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

'设置是否可以调整时间戳，曝光时间和增益'
def enable_embedded_timestamp(cam, enable_timestamp, enable_exposure, enable_gain):
    embedded_info = cam.getEmbeddedImageInfo()
    if embedded_info.available.timestamp:
        cam.setEmbeddedImageInfo(timestamp = enable_timestamp, exposure = enable_exposure, gain = enable_gain)
        if enable_timestamp:
            print('\nTimeStamp、exposure、gain is enabled.\n')
        else:
            print('\nTimeStamp、exposure、gain is disabled.\n')

def grab_images(cam, num_images_to_grab):
    prev_ts = None
    '修改曝光时间和增益'
    PyCapture2.PROPERTY_TYPE(AUTO_EXPOSURE = 1, GAIN = 13)
    for i in range(num_images_to_grab):
        try:
            image = cam.retrieveBuffer()

            alist = []
            img_rows = image.getRows()
            img_cols = image.getCols()
            img_data = image.getData()[::-1]
            for i in range(img_rows):
                alist.append(img_data[i * img_cols:(i + 1) * img_cols])

        except PyCapture2.Fc2error as fc2Err:
            print('Error retrieving buffer : %s' % fc2Err)
            continue

        ts = image.getTimeStamp()
        if prev_ts:
            diff = (ts.cycleSeconds - prev_ts.cycleSeconds) * 8000 + (ts.cycleCount - prev_ts.cycleCount)
            print('Timestamp [ %d %d ] - %d' % (ts.cycleSeconds, ts.cycleCount, diff))
        prev_ts = ts
    return alist

'计算质心'
def zhixin(data):
    b = data.sum()
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            a = data(i,j)*np.array([i,j])/b
    return a



# Example Main
#

# Print PyCapture2 Library Information
print_build_info()

# Ensure sufficient cameras are found
bus = PyCapture2.BusManager()
num_cams = bus.getNumOfCameras()
print('Number of cameras detected: ', num_cams)
if not num_cams:
    print('Insufficient number of cameras. Exiting...')
    exit()

# Select camera on 0th index
c = PyCapture2.Camera()
uid = bus.getCameraFromIndex(0)
c.connect(uid)
print_camera_info(c)

# Enable camera embedded timestamp
enable_embedded_timestamp(c, True, True, True)

print('Starting image capture...')
c.startCapture()
data = grab_images(c, 1)
print(data)

pd_data = pd.DataFrame(data)
pd_data.to_csv('ccd数据.csv', header=None)
c.stopCapture()

# Disable camera embedded timestamp
enable_embedded_timestamp(c, False, False, False)
c.disconnect()

print('Done!\n')


data = np.loadtxt(open("ccd数据.csv"), delimiter=",", skiprows=0)
print("已读入ccd数据\n")

plt.imshow(data, cmap='gray', interpolation='nearest', vmin=0, vmax=255)
plt.show()

zhi = zhixin(data)
print('质心坐标 : %s' % zhi)
