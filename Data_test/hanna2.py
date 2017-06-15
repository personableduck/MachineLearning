# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 21:44:28 2016

@author: 1234
"""

import csv
import matplotlib.pyplot as plt  

with open('C:/Users/1234/Downloads/Individual_Play_Data_KYWA.csv', 'rb') as f:
    reader = csv.reader(f)
    your_list = list(reader)
    
abstact_list=[]

for index in range (len(your_list)):
    if your_list[index][0] == 'fluorine1004@naver.com': #need to change email
            abstact_list.append(your_list[index])


new_order=sorted(abstact_list, key=lambda abstact_list: abstact_list[4])

delete_minus_score=[]

for index in range (len(new_order)):
    if float(new_order[index][1]) > 0:
        delete_minus_score.append(new_order[index])

picture_mat_x=[]

for index2 in range (len(delete_minus_score)):  
    picture_mat_x.append(delete_minus_score[index2][1])  #change number [1]=Score,[2]=palytime,[3]=gmaeId,[4]=gamestart,[5]=gameend

picture_mat_y=[]

for index2 in range (len(delete_minus_score)):  
    picture_mat_y.append(delete_minus_score[index2][2])  

maxvalue=[]

for index in range (len(picture_mat_y)):
    maxvalue.append(float(picture_mat_y[index]))


picture_mat_z=[]

count=0

for index2 in range (len(delete_minus_score)): 
    count=count+1
    picture_mat_z.append(count)

plt.figure(0)    
b=plt.scatter(picture_mat_z,picture_mat_y,c=picture_mat_z)  

plt.figure(1)   
c=plt.scatter(picture_mat_x,picture_mat_y,c=picture_mat_z) 
plt.ylim(0, max(maxvalue))
plt.xlabel('Score')
plt.ylabel('PlayTime')
plt.title('User:fluorine1004@naver.com')
plt.colorbar(b,orientation='horizontal')


picture_mat_x=[]
count=0

for index2 in range (len(delete_minus_score)): 
    count=count+1
    picture_mat_x.append(count)

picture_mat_y=[]

for index2 in range (len(delete_minus_score)):  
    picture_mat_y.append(delete_minus_score[index2][2])  

plt.figure(2) 
plt.plot(picture_mat_x,picture_mat_y,'ro-')  
plt.xlabel('Trial_#')
plt.ylabel('PlayTime')
plt.title('User:fluorine1004@naver.com')
plt.show()


picture_mat_x=[]
count=0

for index2 in range (len(delete_minus_score)): 
    count=count+1
    picture_mat_x.append(count)

picture_mat_y=[]

for index2 in range (len(delete_minus_score)):  
    picture_mat_y.append(delete_minus_score[index2][1])  

plt.figure(3) 
plt.plot(picture_mat_x,picture_mat_y,'bo-')  
plt.xlabel('Trial_#')
plt.ylabel('Score')
plt.title('User:fluorine1004@naver.com')
plt.show()

picture_mat_x=[]

count=0

for index2 in range (len(delete_minus_score)/5): 
    count=count+1
    picture_mat_x.append(count)

if len(delete_minus_score)%5 > 0:
    picture_mat_x.append(len(delete_minus_score)/5+1)

picture_mat_y=[]

for index2 in range (len(delete_minus_score)/5):  
    sum_average=0
    new_index=5*index2
    for index11 in range (0,5):
        sum_average=float(delete_minus_score[new_index+index11][2])+sum_average
    picture_mat_y.append(sum_average/5)
    
if len(delete_minus_score)%5 > 0:
    sum_average=0
    for index2 in range (len(delete_minus_score)%5):
        sum_average=float(delete_minus_score[(len(delete_minus_score)-len(delete_minus_score)%5)+index2][2])+sum_average
    picture_mat_y.append(sum_average/(len(delete_minus_score)%5))     

plt.figure(4) 
plt.plot(picture_mat_x,picture_mat_y,'ro-')  
plt.xlabel('Trial_#(average5)')
plt.ylabel('PlayTime(average5)')
plt.title('User:fluorine1004@naver.com')
plt.show()

picture_mat_x=[]

count=0

for index2 in range (len(delete_minus_score)/5): 
    count=count+1
    picture_mat_x.append(count)

if len(delete_minus_score)%5 > 0:
    picture_mat_x.append(len(delete_minus_score)/5+1)

picture_mat_y=[]

for index2 in range (len(delete_minus_score)/5):  
    sum_average=0
    new_index=5*index2
    for index11 in range (0,5):
        sum_average=float(delete_minus_score[new_index+index11][1])+sum_average
    picture_mat_y.append(sum_average/5)
    
if len(delete_minus_score)%5 > 0:
    sum_average=0
    for index2 in range (len(delete_minus_score)%5):
        sum_average=float(delete_minus_score[(len(delete_minus_score)-len(delete_minus_score)%5)+index2][1])+sum_average
    picture_mat_y.append(sum_average/(len(delete_minus_score)%5))    

plt.figure(5) 
plt.plot(picture_mat_x,picture_mat_y,'bo-')  
plt.xlabel('Trial_#(average5)')
plt.ylabel('Score(average5)')
plt.title('User:fluorine1004@naver.com')
plt.show()