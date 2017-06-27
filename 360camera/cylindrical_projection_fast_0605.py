import glob

import stitch
import seam as seamfunc
import numpy as np
import cv2
import findTransformMatrix as findTM
import cylindrical_projection as cyproj

import imutils
from timeit import default_timer as timer


num = 1

cap_0 = cv2.VideoCapture(0)
cap_0.set(3,1280)
cap_0.set(4,720)
cap_0.set(5,30)


cap_1 = cv2.VideoCapture(1)
cap_1.set(3,1280)
cap_1.set(4,720)
cap_1.set(5,30)

cap_2 = cv2.VideoCapture(2)
cap_2.set(3,1280)
cap_2.set(4,720)
cap_2.set(5,30)

cap_3 = cv2.VideoCapture(3)
cap_3.set(3,1280)
cap_3.set(4,720)
cap_3.set(5,30)

cap_4 = cv2.VideoCapture(4)
cap_4.set(3,1280)
cap_4.set(4,720)
cap_4.set(5,30)


cap_5 = cv2.VideoCapture(5)
cap_5.set(3,1280)
cap_5.set(4,720)
cap_5.set(5,30)


#cap.set(3,1980)
#cap.set(4,1080)
#cap.set(5,30)

#
#img_mapped0 = 0
#img_mapped1 = 0
#img_mapped2 = 0
#img_mapped3 = 0
#img_mapped4 = 0
#img_mapped5 = 0
#
#
#num = 1
#
#
#
#
#
#f = np.float(1150)
#k1 = 0.1
#k2 = 0.1
#
#
#
#theta = lambda y,x : (x - xc) / f
#h = lambda y,x : (y - yc) / f
#
#xcap = lambda y,x : np.sin(theta(y,x))
#ycap = lambda y,x : h(y,x)
#zcap = lambda y,x : np.cos(theta(y,x))
#xn = lambda y,x : xcap(y,x) / zcap(y,x)
#yn = lambda y,x : ycap(y,x) / zcap(y,x)
#r = lambda y,x : xn(y,x)**2 + yn(y,x)**2
#
#xd = lambda y,x : xn(y,x) * (1 + k1 * r(y,x) + k2 * (r(y,x)**2))
#yd = lambda y,x : yn(y,x) * (1 + k1 * r(y,x) + k2 * (r(y,x)**2))
#
#toXimg = lambda y,x : (np.floor(f * xd(y,x) + xc)).astype(int)
#toYimg = lambda y,x : (np.floor(f * yd(y,x) + yc)).astype(int)
#
#
#
#
#def makeMap(y, x):
#    #print x, y    
#    
#    ximg = toXimg(y,x)
#    yimg = toYimg(y,x)
#    return np.dstack((ximg, yimg)).astype(np.int16)
#
#
#def makeOriMap(y,x):
#    return np.dstack((x,y)).astype(np.int16)



while (True):

#	ret0, filename_0 = cap_0.read()  
#	ret1, filename_1 = cap_1.read()  
#	ret2, filename_2 = cap_2.read()  
#	ret3, filename_3 = cap_3.read()  
#	ret4, filename_4 = cap_4.read()  
#	ret5, filename_5 = cap_5.read()  
#	#filename = "1_right2.jpg"
#
#	file_list=[filename_0,filename_1,filename_2,filename_3,filename_4,filename_5]
#
#	imageA_0 = file_list[0]
#	imageA_1 = file_list[1]
#	imageA_2 = file_list[2]
#	imageA_3 = file_list[3]
#	imageA_4 = file_list[4]
#	imageA_5 = file_list[5]
#
#	#imageA = imutils.resize(cv2.imread(filename), width=1000)
#
#	img_0 = imageA_0
#	img_1 = imageA_1
#	img_2 = imageA_2
#	img_3 = imageA_3
#	img_4 = imageA_4
#	img_5 = imageA_5
#
#	xdim = img_0.shape[1]
#	ydim = img_0.shape[0]
#
#	xc = xdim / 2
#	yc = ydim / 2
#
#	out_img0 = np.multiply(img_0,0)
#	out_img1 = np.multiply(img_1,0)
#	out_img2 = np.multiply(img_2,0)
#	out_img3 = np.multiply(img_3,0)
#	out_img4 = np.multiply(img_4,0)
#	out_img5 = np.multiply(img_5,0)
#
#
#	map_ori_xy = np.fromfunction(makeOriMap, img_0.shape[:2], dtype=np.int16)
#
#	map_xy = np.fromfunction(makeMap, img_0.shape[:2], dtype=np.int16)
#
#	start = timer()
#	img_mapped0 = cv2.remap(img_0, (map_xy), None, False)
#
#
#	#cv2.imwrite('shot_%d_proj_crop0.jpg'%(num),img_mapped0)
#
#
#
#	map_ori_xy = np.fromfunction(makeOriMap, img_1.shape[:2], dtype=np.int16)
#
#	map_xy = np.fromfunction(makeMap, img_1.shape[:2], dtype=np.int16)
#
#	start = timer()
#	img_mapped1 = cv2.remap(img_1, (map_xy), None, False)
#
#
#	#cv2.imwrite('shot_%d_proj_crop1.jpg'%(num),img_mapped1)
#
#
#
#	map_ori_xy = np.fromfunction(makeOriMap, img_2.shape[:2], dtype=np.int16)
#
#	map_xy = np.fromfunction(makeMap, img_2.shape[:2], dtype=np.int16)
#
#	start = timer()
#	img_mapped2 = cv2.remap(img_2, (map_xy), None, False)
#
#
#	#cv2.imwrite('shot_%d_proj_crop2.jpg'%(num),img_mapped2)
#
#
#
#	map_ori_xy = np.fromfunction(makeOriMap, img_3.shape[:2], dtype=np.int16)
#
#	map_xy = np.fromfunction(makeMap, img_3.shape[:2], dtype=np.int16)
#
#	start = timer()
#	img_mapped3 = cv2.remap(img_3, (map_xy), None, False)
#
#
#	#cv2.imwrite('shot_%d_proj_crop3.jpg'%(num),img_mapped3)
#
#
#
#	map_ori_xy = np.fromfunction(makeOriMap, img_4.shape[:2], dtype=np.int16)
#
#	map_xy = np.fromfunction(makeMap, img_4.shape[:2], dtype=np.int16)
#
#	start = timer()
#	img_mapped4 = cv2.remap(img_4, (map_xy), None, False)
#
#
#	#cv2.imwrite('shot_%d_proj_crop4.jpg'%(num),img_mapped4)
#
#
#
#	map_ori_xy = np.fromfunction(makeOriMap, img_5.shape[:2], dtype=np.int16)
#
#	map_xy = np.fromfunction(makeMap, img_5.shape[:2], dtype=np.int16)
#
#	start = timer()
#	img_mapped5 = cv2.remap(img_5, (map_xy), None, False)



	#boundary dicided

      img_0 = cap_0[25:695, 0:1190]
      img_1 = cap_1[27:697, 110:990]
      img_2 = cap_2[30:700, 101:1160]
      img_3 = cap_3[25:695, 160:1169]
      img_4 = cap_4[27:697, 125:1070]
      img_5 = cap_5[29:699, 220:1280]

	im=[]
	im2=[]
	im3=[]
	im4=[]
	im5=[]

	for index in range(10):
	    
          img_b0 = cap_0[25:695, 1191+(index*3):1194+(index*3)]
          img_b1 = cap_1[27:697, 79+(index*3):82+(index*3)]
	    
          img_b11 = cap_1[27:697, 991+(index*10):1001+(index*10)]
          img_b2  = cap_2[30:700, 0+(index*10):10+(index*10)]
	    
          img_b22 = cap_2[30:700, 1161+(index*10):1171+(index*10)]
          img_b3  = cap_3[25:695, 59+(index*10):69+(index*10)]
	    
          img_b33 = cap_3[25:695, 1170+(index*5):1175+(index*5)]
          img_b4  = cap_4[27:697, 74+(index*5):79+(index*5)]
	    
          img_b44 = cap_4[27:697, 1071+(index*10):1081+(index*10)]
          img_b5  = cap_5[29:699, 119+(index*10):129+(index*10)]

	    out01=cv2.addWeighted(img_b0,1-(0.11*index),img_b1,0.11*index,0)
	    out12=cv2.addWeighted(img_b11,1-(0.11*index),img_b2,0.11*index,0)
	    out23=cv2.addWeighted(img_b22,1-(0.11*index),img_b3,0.11*index,0)
	    out34=cv2.addWeighted(img_b33,1-(0.11*index),img_b4,0.11*index,0)
	    out45=cv2.addWeighted(img_b44,1-(0.11*index),img_b5,0.11*index,0)
	    
	    im.append(out01)
	    im2.append(out12)
	    im3.append(out23)
	    im4.append(out34)
	    im5.append(out45)


	result_stit=np.hstack((img_0,im[0],im[1],im[2],im[3],im[4],im[5],im[6],im[7],im[8],im[9],img_1,im2[0],im2[1],im2[2],im2[3],im2[4],im2[5],im2[6],im2[7],im2[8],im2[9],img_2,im3[0],im3[1],im3[2],im3[3],im3[4],im3[5],im3[6],im3[7],im3[8],im3[9],img_3,im4[0],im4[1],im4[2],im4[3],im4[4],im4[5],im4[6],im4[7],im4[8],im4[9],img_4,im5[0],im5[1],im5[2],im5[3],im5[4],im5[5],im5[6],im5[7],im5[8],im5[9],img_5))
      
      resized=imutils.resize(result_stit,width=1280)
	#conbine horizontal

	#cv2.imwrite('05_29_result_cv2.png',result_stit)

	#result_crop=new_imk.crop((0,30,6136,618))
	#result_crop.save('28result_crop.png')



	#cv2.imshow('frame',after_mask_trans_image_array2.astype(crop_0.dtype))
	cv2.imshow('frame',resized)
 	
	if cv2.waitKey(1) & 0xFF == ord('z'):
		break


cv2.destroyAllWindows()
#leftStream.stop()

cv2.VideoCapture(0).release()
cv2.VideoCapture(1).release()
cv2.VideoCapture(2).release()
cv2.VideoCapture(3).release()
cv2.VideoCapture(4).release()
cv2.VideoCapture(5).release()
