# -*- coding: utf-8 -*-
"""
Created on Mon May 02 21:34:36 2016

@author: DK
"""

from PIL import Image


img = Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/ab_1.png')
img1 = img.crop((75, 0, 565, 480))
img1.save("new1.jpg")

img = Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/ab_2.png')
img2 = img.crop((75, 0, 565, 480))
img2.save("new2.jpg")

#boundary dicided

img1= Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/new1.jpg')
img2= Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/new2.jpg')
img_left = img1.crop((0, 0, 360, 480))
img_right = img2.crop((40, 0, 490, 480))
img_left.save('img_left.png')
img_right.save('img_right.png')

im=[]

for index in range(20):

    img_b1 = img1.crop((360+(index*2), 0, 362+(index*2), 480))
    img_b2 = img2.crop((0+(index*2), 0, 2+(index*2), 480))

    out=Image.blend(img_b1,img_b2, 0.053*index)
    im.append(out)





#combine part

images = [img_left,im[0],im[1],im[2],im[3],im[4],im[5],im[6],im[7],im[8],im[9],im[10],im[11],im[12],im[13],im[14],im[15],im[16],im[17],im[18],im[19],img_right]
widths, heights = zip(*(i.size for i in images))

total_width = sum(widths)
max_height = max(heights)

new_im = Image.new('RGB', (total_width, max_height))

x_offset = 0
for im in images:
  new_im.paste(im, (x_offset,0))
  x_offset += im.size[0]

new_im.save('aa_result_new2.png')