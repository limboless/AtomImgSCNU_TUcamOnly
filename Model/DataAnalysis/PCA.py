import numpy as np
from numpy import linalg as lig
from PIL import Image as image
import os

import matplotlib.pyplot as plt
from PIL import Image

'''read out all specific image filepath'''
def read_all_image_fp(cwd):
	filenamelist=[]
	get_dir= os.listdir(cwd)
	for i in get_dir:
		sub_dir=cwd+'\\'+i
		if os.path.isdir(sub_dir):
			read_all_image_fp(sub_dir)
		else:
			filenamelist.append(sub_dir)
	for h in filenamelist:
		print(h)
	return filenamelist
	

'''for fitting on the windows system(unnecessary)'''
def rawstring_obtainedq(text):
	raw_string=''
	raw_stringset=[]
	escape_dict={'\a':r'\a',
          '\b':r'\b',
         '\c':r'\c',
         '\f':r'\f',
         '\n':r'\n',
         '\r':r'\r',
         '\t':r'\t',
         '\v':r'\v',
         '\'':r'\'',
         '\"':r'\"',
         '\0':r'\0',
         '\1':r'\1',
         '\2':r'\2',
         '\3':r'\3',
         '\4':r'\4',
         '\5':r'\5',
         '\6':r'\6',
         '\7':r'\7',
         '\8':r'\8',
         '\9':r'\9'}
	for j in text:
		for i in j:
			try:
				raw_string+=escape_dict[i]
			except KeyError:
				raw_string+=i
		raw_stringset.append(raw_string)
		raw_string=''
	print(raw_stringset)
	return raw_stringset
	
'''for fitting on the windows system(unnecessary)'''
def rawstring_obtained(text):
	raw_string=''
	escape_dict={'\a':r'\a',
          '\b':r'\b',
         '\c':r'\c',
         '\f':r'\f',
         '\n':r'\n',
         '\r':r'\r',
         '\t':r'\t',
         '\v':r'\v',
         '\'':r'\'',
         '\"':r'\"',
         '\0':r'\0',
         '\1':r'\1',
         '\2':r'\2',
         '\3':r'\3',
         '\4':r'\4',
         '\5':r'\5',
         '\6':r'\6',
         '\7':r'\7',
         '\8':r'\8',
         '\9':r'\9'}
	for i in text:
		try:
			raw_string+=escape_dict[i]
		except KeyError:
			raw_string+=i
	print(raw_string)
	return raw_string

'''combine all file of the data set as a matrix'''
def figure_joint(filenamelist):
	size_of_file=len(filenamelist)
	count=0
	for i in filenamelist:
		im=image.open(i)
		imwoarr=np.array(im)
		[imheight,imwidth]=imwoarr.shape
		imwoarr=imwoarr.reshape(1,imheight*imwidth)
		if count==0:
			imwosetarr=np.zeros([size_of_file,imheight*imwidth],float)
			imwosetarr[count,:]=imwoarr
		else:
			imwosetarr[count,:]=imwoarr
		count+=1
		print(imwoarr)
		print('\n')
	print(imwosetarr)
	print('\n')
	print(imwosetarr.shape)
	return imwosetarr,imheight,imwidth,size_of_file
			
'''main part of principle component analysis, obtain the eigenvalues andthe eigenvectors for the principle component'''			
def pca_alg_eiva_eive(origin_array):
	'''work out the mean vector for every object'''
	mean_value=np.mean(origin_array,1)
	mean_value1=mean_value[:,np.newaxis]
	origin_array1=origin_array-mean_value1
	'''obtain the operator'''
	origin_mat=np.matrix(origin_array1)
	origin_mat_permute=np.matrix(origin_array1.T)
	suqare_example_operator=origin_mat*origin_mat_permute
	'''solve the eignvalue equation'''
	[eigenvalues,eigenvectors]=lig.eig(suqare_example_operator)
	print(eigenvalues)
	print('\n')
	print(eigenvectors)
	return eigenvalues,eigenvectors

'''obtain the principle image'''
def principle_component(eigenvalues,eigenvectors,imwosetarr1,imheight,imwidth):
	'''storing the eigenvalues and engenvextors as the key-value pairs and sorting them according to the eigenvalues'''
	eigendictionary={}
	count=0
	for i in eigenvectors.T:
		eigendictionary[eigenvalues[count]]=i
		count+=1
	sortedeigenvectors=[]
	for key in sorted(eigendictionary.keys()):
		sortedeigenvectors.append(eigendictionary[key])
	sortedeigenvalues=sorted(eigenvalues)
	'''obtaining the principle image'''
	count1=0
	length_of_eigenvalues=len(eigenvalues)
	principle_image=np.zeros([length_of_eigenvalues,imheight,imwidth],float)
	for j in sortedeigenvectors:
		prinvec=j[:,np.newaxis]
		prinvecmat=np.matrix(prinvec)
		imwosetarr=imwosetarr1.T
		imwosetmat=np.matrix(imwosetarr)
		prinimgmat=imwosetmat*prinvecmat.T
		prinimg1=np.array(prinimgmat)
		prinimg=prinimg1.reshape(imheight,imwidth)
		principle_image[count1,:,:]=prinimg
		count1+=1
	return principle_image
	
	

dirpath=input(r'Please offer the specific path of the folder in which the targget files are located and cached, if your computer systems please replace all \ with \\')
fp=read_all_image_fp(dirpath)
[imwoset,imh,imw,primsize]=figure_joint(fp)
[eiva,eive]=pca_alg_eiva_eive(imwoset)
prinimg=principle_component(eiva,eive,imwoset,imh,imw)
print(prinimg.shape)
print('\n')
print(prinimg)
for i in range(primsize):
	plt.imshow(prinimg[i,:,:])
	plt.show()

	
# fp=read_all_image_fp(r'F:\学習の資料\quantum\new laboratory\temporary group\my portfolio\Example\WO')
# [imwoset,imh,imw,primsize]=figure_joint(fp)
# [eiva,eive]=pca_alg_eiva_eive(imwoset)
# prinimg=principle_component(eiva,eive,imwoset,imh,imw)
# print(prinimg.shape)
# print('\n')
# print(prinimg)
# for i in range(primsize):
	# plt.imshow(prinimg[i,:,:])
	# plt.show()







