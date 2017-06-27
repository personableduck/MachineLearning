# -*- coding: utf-8 -*-
"""
Created on Mon May 02 21:34:36 2016

@author: DK
"""

from PIL import Image

#move image1

imgblack=Image.new("RGB",[1139,27],(0,0,0))

img_up1= Image.open('C:/Users/1234/Downloads/6_shots/1.jpg')
img_down= imgblack

images = [img_up1,img_down]

widths, heights = zip(*(i.size for i in images))

total_width = max(widths)
max_height = sum(heights)

new_im = Image.new('RGB', (total_width, max_height))

x_offset = 0
for im in images:
  new_im.paste(im, (x_offset,0))
  x_offset += im.size[0]


#move image2
imgblack2=Image.new("RGB",[1139,47],(0,0,0))

img_up2= Image.open('C:/Users/1234/Downloads/6_shots/3.jpg')

images3 = [img_up2,imgblack2]

widths, heights = zip(*(i.size for i in images3))

total_width = max(widths)
max_height = sum(heights)

new_im3 = Image.new('RGB', (total_width, max_height))

x_offset = 0
for im in images3:
  new_im3.paste(im, (x_offset,0))
  x_offset += im.size[0]

new_im3.save("bb11.jpg")
#newstep

img0 = Image.open('C:/Users/1234/Downloads/6_shots/0.jpg')
#img0 = img.crop((0, 0, 1139, 718))
#img0.save("raw11.jpg")

img1 = new_im
#img1 = img.crop((0, 0, 1139, 745))
#img1.save("raw22.jpg")

img2= Image.open('C:/Users/1234/Downloads/6_shots/2.jpg')

img3= Image.open('C:/Users/1234/Downloads/6_shots/3.jpg')
#boundary dicided

img_0 = img0.crop((0, 0, 949, 718))
img_1 = img1.crop((10, 13, 974, 731))
img_2 = img2.crop((10, 13, 1010, 731))
img_3 = new_im3.crop((10, 26, 1139, 744))
img_4 = img3.crop((-10,0,1139,800))
#img_0.save('img_left.png')
img_4.save('img4.png')



im=[]
im2=[]
im3=[]

for index in range(10):

    img_b0 = img0.crop((949+(index*1), 0, 950+(index*1), 718))
    img_b1 = img1.crop((0+(index*1), 13, 1+(index*1), 731))
    
    img_b11 = img1.crop((974+(index*1), 13, 975+(index*1), 731))
    img_b2  = img2.crop((0+(index*1), 13, 1+(index*1), 731))
    
    img_b22 = img2.crop((1010+(index*1), 0, 1011+(index*1), 718))
    img_b3  = img3.crop((0+(index*1), 0, 1+(index*1), 718))

    out01=Image.blend(img_b0,img_b1, 0.11*index)
    out12=Image.blend(img_b11,img_b2, 0.11*index)
    out23=Image.blend(img_b22,img_b3, 0.11*index)
    
    im.append(out01)
    im2.append(out12)
    im3.append(out23)


#conbine horizontal

images = [img_0,im[0],im[1],im[2],im[3],im[4],im[5],im[6],im[7],im[8],im[9],img_1,im2[0],im2[1],im2[2],im2[3],im2[4],im2[5],im2[6],im2[7],im2[8],im2[9],img_2,im3[0],im3[1],im3[2],im3[3],im3[4],im3[5],im3[6],im3[7],im3[8],im3[9],img_3]
widths, heights = zip(*(i.size for i in images))


total_width = sum(widths)
max_height = max(heights)

new_imk = Image.new('RGB', (total_width, max_height))

x_offset = 0
for im in images:
  new_imk.paste(im, (x_offset,0))
  x_offset += im.size[0]

new_imk.save('05_21_result_num4.png')