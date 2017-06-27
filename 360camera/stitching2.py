# -*- coding: utf-8 -*-
"""
Created on Mon May 02 23:49:56 2016

@author: 1234
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May 02 21:34:36 2016

@author: DK
"""

from PIL import Image

img = Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/cc_1.png')
img1 = img.crop((80, 80, 880, 640))
img1.save("st001.jpg")

img = Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/cc_2.png')
img2 = img.crop((80, 80, 880, 640))
img2.save("st002.jpg")

img = Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/cc_1.png')
img1 = img.crop((80, 80, 796, 640))
img1.save("ct001.jpg")

img = Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/cc_2.png')
img2 = img.crop((340, 80, 880, 640))
img2.save("ct002.jpg")

#boundary dicided

img1= Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/cc_1.png')
img2= Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/cc_2.png')
img_left = img1.crop((80, 80, 771, 640))
img_right = img2.crop((370, 80, 880, 640))


im=[]

for index in range(11):

    img_b1 = img1.crop((771+(index*5), 80, 776+(index*5), 640))
    img_b2 = img2.crop((315+(index*5), 80, 320+(index*5), 640))

    out=Image.blend(img_b1,img_b2, 0.1*index)
    im.append(out)





#combine part

images = [img_left,im[0],im[1],im[2],im[3],im[4],im[5],im[6],im[7],im[8],im[9],im[10],img_right]
widths, heights = zip(*(i.size for i in images))

total_width = sum(widths)
max_height = max(heights)

new_im = Image.new('RGB', (total_width, max_height))

x_offset = 0
for im in images:
  new_im.paste(im, (x_offset,0))
  x_offset += im.size[0]

new_im.save('t_result_new1.png')