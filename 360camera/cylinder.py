# -*- coding: utf-8 -*- 
"""
Created on Mon May 02 21:28:19 2016

@author: DK_stitching
"""

import cv2

import numpy as np

imageA = cv2.imread("C:/Users/1234/Desktop/embedded/DK_stitching/frame5.jpg")


img = imageA

out_img = np.multiply(img,0)

xdim = img.shape[1]
ydim = img.shape[0]

xc = xdim / 2
yc = ydim / 2


f = float(1150)
k1 = 0.1
k2 = 0.1



for y in range(1,ydim+1):
    for x in range(1,xdim+1):
        theta = (x - xc) / f
        h = (y - yc) / f
        
        xcap = np.sin(theta)
        ycap = h
        zcap = np.cos(theta)
        xn = xcap / zcap
        yn = ycap / zcap
        r = xn**2 + yn**2
        
        xd = xn * (1 + k1 * r + k2 * (r**2))
        yd = yn * (1 + k1 * r + k2 * (r**2))
        
        ximg = np.int(np.floor(f * xd + xc))
        yimg = np.int(np.floor(f * yd + yc))
        
        ximg_nf = (f * xd + xc)
        yimg_nf = (f * yd + yc)
        
        
        
        if ximg > 0 and ximg <= xdim and yimg > 0 and yimg <= ydim:
            #print y, yimg, x, ximg
            #break
            out_img[y-1,x-1,:] = img[yimg-1,ximg-1,:]
            
            


cv2.imwrite("frame5_c.jpg",out_img)
