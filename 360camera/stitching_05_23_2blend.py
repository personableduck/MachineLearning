# -*- coding: utf-8 -*-
"""
Created on Mon May 02 21:34:36 2016

@author: DK
"""

from PIL import Image

#move image1



#move image2

#newstep

img0 = Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/frame0.jpg')
img0 = img0.crop((0, 0, 1139, 745))
#img0.save("raw11.jpg")

img1 = Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/frame1.jpg')
#img1 = img.crop((0, 0, 1139, 745))
#img1.save("raw22.jpg")

img2= Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/frame2.jpg')

img3= Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/frame3.jpg')
img3 = img3.crop((0, 0, 1139, 745))

img4= Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/frame4.jpg')
img4 = img4.crop((0, 0, 1139, 800))

img5= Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/frame5.jpg')
img5 = img5.crop((0, 0, 1139, 745))
#boundary dicided

img_0 = img0.crop((0, 0, 949, 718))
img_1 = img1.crop((11, 13, 974, 731))
img_2 = img2.crop((11, 13, 1010, 731))
img_3 = img3.crop((11, 26, 1120, 744))
img_4 = img4.crop((11, 82, 949, 800))
img_5 = img5.crop((11, 30, 1139, 748))
#img_0.save('img_left.png')




im=[]
im2=[]
im3=[]
im4=[]
im5=[]

for index in range(10):

    img_b0 = img0.crop((950+(index*1), 0, 951+(index*1), 718))
    img_b1 = img1.crop((0+(index*1), 13, 1+(index*1), 731))
    
    img_b11 = img1.crop((975+(index*1), 13, 976+(index*1), 731))
    img_b2  = img2.crop((0+(index*1), 13, 1+(index*1), 731))
    
    img_b22 = img2.crop((1011+(index*1), 13, 1012+(index*1), 731))
    img_b3  = img3.crop((0+(index*1), 26, 1+(index*1), 744))
    
    img_b33 = img3.crop((1121+(index*1), 26, 1122+(index*1), 744))
    img_b4  = img4.crop((0+(index*1), 82, 1+(index*1), 800))
    
    img_b44 = img4.crop((950+(index*1), 82, 951+(index*1), 800))
    img_b5  = img5.crop((0+(index*1), 30, 1+(index*1), 748))
    

    out01=Image.blend(img_b0,img_b1, 0.11*index)
    out12=Image.blend(img_b11,img_b2, 0.11*index)
    out23=Image.blend(img_b22,img_b3, 0.11*index)
    out34=Image.blend(img_b33,img_b4, 0.11*index)
    out45=Image.blend(img_b44,img_b5, 0.11*index)
    
    im.append(out01)
    im2.append(out12)
    im3.append(out23)
    im4.append(out34)
    im5.append(out45)
    


#conbine horizontal

images = [img_0,im[0],im[1],im[2],im[3],im[4],im[5],im[6],im[7],im[8],im[9],img_1,im2[0],im2[1],im2[2],im2[3],im2[4],im2[5],im2[6],im2[7],im2[8],im2[9],img_2,im3[0],im3[1],im3[2],im3[3],im3[4],im3[5],im3[6],im3[7],im3[8],im3[9],img_3,im4[0],im4[1],im4[2],im4[3],im4[4],im4[5],im4[6],im4[7],im4[8],im4[9],img_4,im5[0],im5[1],im5[2],im5[3],im5[4],im5[5],im5[6],im5[7],im5[8],im5[9],img_5]
widths, heights = zip(*(i.size for i in images))


total_width = sum(widths)
max_height = max(heights)

new_imk = Image.new('RGB', (total_width, max_height))

x_offset = 0
for im in images:
  new_imk.paste(im, (x_offset,0))
  x_offset += im.size[0]

new_imk.save('05_28_result_num4.png')

result_crop=new_imk.crop((0,30,6136,618))
result_crop.save('28result_crop.png')