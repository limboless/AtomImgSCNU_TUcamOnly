# from osgeo import gdal
# import os
# import numpy as np
 
# source_dir = "D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll\\Image1006242.tif"
# target_dir = "D:\\vscode\\becAtomImg SCNU20220425\\AtomImgSCNU\\Model\\Instruments\\Camera\\TUCamdll"
 
# class_list = os.listdir(source_dir)
# for cla in class_list:
#     img_list = os.listdir(os.path.join(source_dir,cla))
#     for img in img_list:
#         # 拆分文件名
#         name = []
#         portion = os.path.splitext(img)  # 把文件名拆分为名字和后缀
#         if portion[1] == ".tiff":
#             name = portion[0]
#             image_path = source_dir + "\\" + cla + "\\" + img
#             image_name = target_dir + "\\" + cla + "\\" + name + ".npy"
 
#             # 读取数组
#             image = gdal.Open(image_path)  # 读取栅格数据
#             img_array = image.ReadAsArray()
 
#             # 创建文件夹
#             if not os.path.exists(target_dir + "\\" + cla):
#                 os.makedirs(target_dir + "\\" + cla)
 
#             K = np.array(img_array)
#             # 存储npy文件
#             np.save(image_name, K)


import datetime
timestamp = datetime.datetime.now()
formatted_timestamp = str(timestamp).replace(":", "").replace("-", "").replace(".", "_")
print(formatted_timestamp)

