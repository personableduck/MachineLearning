# -*- coding: utf-8 -*-
"""
Created on Tue May 03 00:15:20 2016

@author: 1234
"""

img1= Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/result_new1.png')
img2= Image.open('C:/Users/1234/Desktop/embedded/DK_stitching/result_new2.png')
img_left = img1.crop((0, 0, 1100, 560))
img_right = img2.crop((645, 0, 1195, 560))
img_left.save("dd1.jpg")
img_right.save("dd2.jpg")

images = [img_left,img_right]
widths, heights = zip(*(i.size for i in images))

total_width = sum(widths)
max_height = max(heights)

new_im = Image.new('RGB', (total_width, max_height))

x_offset = 0
for im in images:
  new_im.paste(im, (x_offset,0))
  x_offset += im.size[0]

new_im.save('result_new1_2.png')