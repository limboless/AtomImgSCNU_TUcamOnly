# from osgeo import gdal
# import os
# import numpy as np
# #import cv2

# # source_dir = "C:\\Users\\Administrator\\Desktop\\SAR_specific_models-master\\data\\FUSARship-1-512-tiff"
# # target_dir = "C:\\Users\\Administrator\\Desktop\\SAR_specific_models-master\\data\\slc_data"
# source_dir = "D:\\vscode\\becAtomImg SCNU20220425\\becAtomImg SCNU20220425\\AtomImgSCNU\\Debug\\tiffile"
# #source_dir = "D:\\vscode\\becAtomImg SCNU20220425\\becAtomImg SCNU20220425\\AtomImgSCNU\Debug\\tiffile\\Image1833825.tif"
# target_dir = "D:\\vscode\\becAtomImg SCNU20220425\\becAtomImg SCNU20220425\\AtomImgSCNU\\Debug\\tiffile"

# class_list = os.listdir(source_dir)
# for cla in class_list:
#     img_list = os.listdir(os.path.join(source_dir,cla))
#     for img in img_list:
#         # 拆分文件名
#         name = []
#         portion = os.path.splitext(img)  # 把文件名拆分为名字和后缀
#         if portion[1] == ".tif":
#             name = portion[0]
#             image_path = source_dir + "\\" + cla + "\\" + img
#             image_name = target_dir + "\\" + cla + "\\" + name + ".npy"
 
#             # 读取数组
#             image = gdal.Open(image_path)  # 读取栅格数据
#             #image = cv2.imread(image_path)
#             img_array = image.ReadAsArray()
 
#             # 创建文件夹
#             if not os.path.exists(target_dir + "\\" + cla):
#                 os.makedirs(target_dir + "\\" + cla)
 
#             K = np.array(img_array)
#             # 存储npy文件
#             np.save(image_name, K)
 

# import numpy as np
# from PIL import Image
# # 读取tif文件
# tif_file = "AtomImgSCNU\\Debug\\tiffile\\Image1833825.tif"
# tif_image = Image.open(tif_file)
# # 将tif文件转换为numpy数组
# tif_array = np.array(tif_image)
# # 保存为npy文件
# npy_file = "AtomImgSCNU\\Debug\\tiffile\\Image1833825.npy"
# np.save(npy_file, tif_array)



import numpy as np
from PIL import Image

def split_list(image):
    """
    design for Chameleon camera, take the  image to array
    :param image: list
    :return:
    """
    alist = []
    img_rows = image.getRows()
    img_cols = image.getCols()
    img_data = image.getData()[::-1]
    # print(img_data)
    # print(img_rows)
    # print(img_cols)
    for i in range(img_rows):
        alist.append(img_data[i*img_cols:(i+1)*img_cols])
    # print(np.array(alist))

    return np.array(alist)


# 读取tif文件
tif_file = "D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\test\\testimage\\WOtif.tif"
tif_image = Image.open(tif_file)
ceshi_image = split_list(tif_image)
print(ceshi_image)

# 保存为npy文件
# tmpWI = "D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\pic\\tif\\WO230311_124135_695731.npy"
# np.save(tmpWI, tif_array)




