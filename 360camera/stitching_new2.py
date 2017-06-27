# -*- coding: utf-8 -*-
"""
Created on Mon May 02 21:34:36 2016

@author: DK
"""

from PIL import Image

#move image

imgblack=Image.new("RGB",[1139,27],(0,0,0))
imgblack.save('black_image.png')

img_up= Image.open('C:/Users/1234/Downloads/6_shots/1.jpg')
img_down= imgblack

images = [img_up,img_down]
widths, heights = zip(*(i.size for i in images))

total_width = max(widths)
max_height = sum(heights)

new_im = Image.new('RGB', (total_width, max_height))

x_offset = 0
for im in images:
  new_im.paste(im, (x_offset,0))
  x_offset += im.size[0]

new_im.save('66_1.png')

#newstep

img = Image.open('C:/Users/1234/Downloads/6_shots/0.jpg')
img1 = img.crop((0, 0, 1139, 718))
img1.save("raw11.jpg")

img = new_im
img2 = img.crop((0, 0, 1139, 745))
img2.save("raw22.jpg")

#boundary dicided

img_left = img1.crop((0, 0, 949, 718))
img_right = img2.crop((10, 13, 1139, 731))
img_left.save('img_left.png')
img_right.save('img_right.png')



im=[]

for index in range(10):

    img_b1 = img1.crop((949+(index*1), 0, 950+(index*1), 718))
    img_b2 = img2.crop((0+(index*1), 13, 1+(index*1), 731))

    out=Image.blend(img_b1,img_b2, 0.11*index)
    im.append(out)



#conbine horizontal

images = [img_left,im[0],im[1],im[2],im[3],im[4],im[5],im[6],im[7],im[8],im[9],img_right]
widths, heights = zip(*(i.size for i in images))


total_width = sum(widths)
max_height = max(heights)

new_im2 = Image.new('RGB', (total_width, max_height))

x_offset = 0
for im in images:
  new_im2.paste(im, (x_offset,0))
  x_offset += im.size[0]

new_im2.save('05_21_result_new5.png')