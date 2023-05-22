import numpy as np
import copy
from PIL import Image

#########################
# 一次处理存放在相同文件夹下的
# 将tif文件与此python文件放在同一文件夹下
# 读取tif文件
for i in range(100) :
    img_name1 = "Image"  #将Image修改为图片命名。例如图片名为843_00001,则修改Image为843
    if i < 10 :
        ImgName = './' + str(img_name1) + '_' +  '0000' + str(i) + '.tif'
    elif i == 100:
        ImgName = './' + str(img_name1) + '_' +  '00' + str(i) + '.tif'
    else:
        ImgName = './' + str(img_name1) + '_' + '000' + str(i) + '.tif'
    # tmpWI1 = "D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\test\\843_00002.tif"
    tmpWI = Image.open(ImgName)
    
    # 将tif文件转换为numpy数组
    tif_array = np.array(tmpWI)
    # 保存为npy文件
    if i < 10 :
        ImgNamen = './' + str(img_name1) + '_' +  '0000' + str(i) + '.npy'
    elif i == 100:
        ImgName = './' + str(img_name1) + '_' +  '00' + str(i) + '.npy'
    else:
        ImgNamen = './' + str(img_name1) + '_' + '000' + str(i) + '.npy'
    # ImgNamen = "D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\test\\843_00002.npy"
    np.save(ImgNamen, tif_array)




# # 读取npy格式深度图文件
# # depth_map = np.load('depth_map.npy'
# # depth_tif2npy = np.load('AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\pic\\ODimage\\OD230311_130506_923280.npy')
# data_tif2npy = np.load('D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\pic\\npyODimage\\OD230311_130506_923280.npy')
# depth_npy = np.load('D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\pic\\yuan_npy\\Data221027_174319_111.npy')
# # 读取浮点类型的npy文件
# data_float = data_tif2npy
# # 将浮点类型转换为整型
# data_int = data_float.astype(np.uint16)
# # 保存整型的npy文件
# np.save('D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\pic\\float2int\\uint16ODData221027_111.npy', data_int)

# modi 修改图像
#settings.absimgData[3] = Calc_absImg(settings.absimgData[0], settings.absimgData[1], settings.absimgData[2], 0)

# def Calc_absImg(tmpWI, tmpWO, tmpBG):
#    withoutAtom = tmpWO - tmpBG
#    withoutAtom[withoutAtom == 0] = 1 #分母为零将导致计算无意义，先改为 1，后将这些点的数值改为 1，取对数后为 0
#    img = (tmpWI - tmpBG) / (tmpWO - tmpBG)
#    img[withoutAtom == 0] = 1
#    img[img <= 0] = 1
#    img[img >= 1] = 1
#    img = -np.log(img)
#    return img



# # 读取tif文件
# tmpWI1 = "D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\pic\\tif\\WO230311_124135_695731.tif"
# tmpWI = Image.open(tmpWI1)
# # 将tif文件转换为numpy数组
# tif_array = np.array(tmpWI)
# # 保存为npy文件
# tmpWI = "D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\pic\\tif\\WO230311_124135_695731.npy"
# np.save(tmpWI, tif_array)

# # 读取tif文件
# tmpWI1 = "D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\pic\\tif\\BG230311_130036_713378.tif"
# tmpWI = Image.open(tmpWI1)
# # 将tif文件转换为numpy数组
# tif_array = np.array(tmpWI)
# # 保存为npy文件
# tmpWI = "D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\pic\\tif\\BG230311_130036_713378.npy"
# np.save(tmpWI, tif_array)

# ##################

# tmpWI = np.load('D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\pic\\tif\\WI230311_125746_303891.npy')
# tmpWO = np.load('D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\pic\\tif\\WO230311_124135_695731.npy')
# tmpBG = np.load('D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\pic\\tif\\BG230311_130036_713378.npy')

# def Calc_absImg(tmpWI, tmpWO, tmpBG):
#     withoutAtom = tmpWO - tmpBG
#     withoutAtomTest = copy.deepcopy(withoutAtom)  #用于判断何处数值为0
#     withoutAtom[withoutAtomTest == 0] = 1
#     img = (tmpWI - tmpBG) / (withoutAtom)
#     # use more efficient way to restrict image value
#     img[withoutAtomTest == 0] = 1
#     img[img <= 0] = 1
#     img[img >= 1] = 1
#     img = -np.log(img)
#     return img

# ODimage = Calc_absImg(tmpWI, tmpWO, tmpBG)
# npy_file = "D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\pic\\tif\\ODXIN.npy"
# np.save(npy_file, ODimage)

# # 裁剪需要的部分
# data_cropped = data_tif2npy[:960, :1280]

# # 将裁剪后的数据保存到新的npy文件
# np.save('D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\pic\\cropped_data.npy', data_cropped)

# # 查看图片矩阵
# print(depth_npy)
# print(data_tif2npy)
# print(data_int)
# print(data_cropped)
# # 查看格式
# print('depth_npy.dtype is ',depth_npy.dtype)
# print('data_tif2npy.dtype is ',data_tif2npy.dtype)
# # print(data_int.dtype)
# # print(data_cropped.dtype)
# # 查看尺寸
# print('depth_npy.shape is ',depth_npy.shape)
# print('data_tif2npy.shape is ',data_tif2npy.shape)
# # print(data_int.shape)
# # print(data_cropped.shape)

# # uint16
# uint8
# (2048, 2048)
# (960, 1280)




# filename = "end program.txt"
# try:
#     with open(filename, 'r') as f:
#         content = f.read().strip()
#         if content == '1':
#             print("1")
#         elif content == '0':
#             print("0")
#         else:
#             print("Invalid content in file.")
# except FileNotFoundError:
#     # If the file does not exist, create it and write "1" to it.
#     with open(filename, 'w') as f:
#         f.write("1")
#         print("1")
