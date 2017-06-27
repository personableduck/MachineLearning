# -*- coding: utf-8 -*-
"""
Created on Mon May 02 21:34:36 2016

@author: DK
"""

from PIL import Image


img = Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/13_1.png')
img1 = img.crop((100, 0, 1180, 720))
img1.save("c_raw1.jpg")

img = Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/13_2.png')
img2 = img.crop((80, 0, 1200, 720))
img2.save("c_raw2.jpg")

img = Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/13_3.png')
img3 = img.crop((60, 0, 1220, 720))
img3.save("c_raw3.jpg")

#boundary dicided

img1= Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/c_raw1.jpg')
img2= Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/c_raw2.jpg')
img3= Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/c_raw3.jpg')

img_left = img1.crop((0, 0, 996, 720))
img_middle = img2.crop((110, 0, 1040, 720))
img_right = img3.crop((75, 0, 1160, 720))

img_left.save('img_left.png')
img_middle.save('img_middle.png')
img_right.save('img_right.png')

im=[]

for index in range(20):

    img_b1 = img1.crop((996+(index*2), 0, 998+(index*2), 720))
    img_b2 = img2.crop((70+(index*2), 0, 72+(index*2), 720))

    out=Image.blend(img_b1,img_b2, 0.053*index)
    im.append(out)

im2=[]

for index in range(20):

    img_b1 = img2.crop((1040+(index*2), 0, 1042+(index*2), 720))
    img_b2 = img3.crop((35+(index*2), 0, 37+(index*2), 720))

    out=Image.blend(img_b1,img_b2, 0.053*index)
    im2.append(out)



#combine part

images = [img_left,im[0],im[1],im[2],im[3],im[4],im[5],im[6],im[7],im[8],im[9],im[10],im[11],im[12],im[13],im[14],im[15],im[16],im[17],im[18],im[19],img_middle,im2[0],im2[1],im2[2],im2[3],im2[4],im2[5],im2[6],im2[7],im2[8],im2[9],im2[10],im2[11],im2[12],im2[13],im2[14],im2[15],im2[16],im2[17],im2[18],im2[19],img_right]
widths, heights = zip(*(i.size for i in images))

total_width = sum(widths)
max_height = max(heights)

new_im = Image.new('RGB', (total_width, max_height))

x_offset = 0
for im in images:
  new_im.paste(im, (x_offset,0))
  x_offset += im.size[0]

new_im.save('05_13_result_new5.png')