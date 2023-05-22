#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fringe Removal

"""



import numpy as np
from scipy.linalg import solve


'''Reduce matrix dimension,The first dimension represents the serial number of the image, and the second dimension is the elements of each image'''
def dimreduction(imgset):
    array3d = imgset.copy()
    n = array3d.shape[0]
    height = array3d.shape[1]
    width  = array3d.shape[2]
    array2d = array3d.reshape(n,height*width)
    return array2d

'''get the correlation matrix of the basis'''
def cormatrix(imgset):
    array3d = imgset.copy()
    array2d = dimreduction(array3d)
    ans = np.dot(array2d,array2d.T)
    return ans

'''extend basis of imgW0 and cut its edge'''
def extend(imgsetRaw,shiftMax,edge):
    array3d = imgsetRaw.copy()
    n0 = array3d.shape[0]
    height = array3d.shape[1]
    width = array3d.shape[2]
    mul = np.square(2*shiftMax+1)
    n = mul*n0
    imgset = np.zeros((n,height - 2*edge,width - 2*edge))
    temp = 0
    for i in list(range(0,n0)):
        img = np.squeeze(array3d[i,:,:])
        for j in list(range(-shiftMax,shiftMax+1)):
            for k in list(range(-shiftMax,shiftMax+1)):
                imgset[temp,:,:] = img[edge-j:height-edge-j,edge-k:width-edge-k]
                temp+=1
    return imgset

'''cut the edge of imgWi'''
def cutedge(imgsetRaw,edge):
    array3d = imgsetRaw.copy()
    n0 = array3d.shape[0]
    height = array3d.shape[1]
    width = array3d.shape[2]
    imgset = np.zeros((n0,height - 2*edge,width - 2*edge))
    temp = 0
    for i in list(range(0,n0)):
        img = np.squeeze(array3d[i,:,:])
        imgset[temp,:,:] = img[edge:height-edge,edge:width-edge]
        temp+=1
    return imgset


'''maskSize should be an even number,indicating the side length of the mask. The selected point is at the bottom right'''
def mask(imgset,xc,yc,maskSize):
    array3d = imgset.copy()
    n = array3d.shape[0]
    height = array3d.shape[1]
    width = array3d.shape[2]    
    xw = int(maskSize/2)
    yw = int(maskSize/2)
    mask = np.ones((height,width)) 
    mask[(yc-yw):(yc+yw),(xc-xw):(xc+xw)] = 0
    for i in list(range(0,n)):
        array3d[i,:,:]*mask
    return array3d 



'''get the projection matrix'''
def projectionM(imgWoset,imgWiset,i):#Starting from 1, operate on a single imgwi
    Wiarray2d = dimreduction(imgWiset)
    vector = np.squeeze(Wiarray2d[i,:])
    Woarray3d = imgWoset.copy()
    Woarray2d = dimreduction(Woarray3d)
    p = np.dot(Woarray2d,vector)
    return p


def FR(imgWoset,imgWiset,shiftMax,edge,Xc,Yc,maskSize):


    imgWoset = extend(imgWoset,shiftMax,edge)
    imgWiset = cutedge(imgWiset,edge)
    imgWoCset = 0*imgWiset
    
    imgWosetMask = mask(imgWoset,Xc,Yc,maskSize)
    imgWisetMask = mask(imgWiset,Xc,Yc,maskSize)
    
    n = imgWiset.shape[0]
    N = imgWoset.shape[0]
    
    C = cormatrix(imgWosetMask)

    
    for i in list(range(0,n)):
        P = projectionM(imgWosetMask, imgWisetMask, i)
        beta = solve(C,P)
        print(i)
    
        for j in list(range(0,N)):
            imgWoCset[i,:,:] += beta[j]*imgWoset[j,:,:]
            
 
    return imgWoCset

