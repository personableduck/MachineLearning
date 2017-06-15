# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 21:44:28 2016

@author: 1234
"""

import csv
import matplotlib.pyplot as plt  

with open('C:/Users/1234/Downloads/[Input]Individual_Play_Data_KYWA.csv', 'rb') as f:
    reader = csv.reader(f)
    your_list = list(reader)

with open('C:/Users/1234/Downloads/[Output]Individual_User_Maximum_Score_KYWA.csv', 'rb') as t:
    reader2 = csv.reader(t)
    user_list = list(reader2)

your_list_dlt=[]
for index in range (1,len(your_list)):
    your_list_dlt.append(your_list[index])

user_list_dlt=[]
for index in range (1,len(user_list)):
    user_list_dlt.append(user_list[index])

#delete_minuse_process
delete_minus_score=[]

for index in range (len(your_list_dlt)):
    if float(your_list_dlt[index][1]) > 0:
        delete_minus_score.append(your_list_dlt[index])
#-------------------
#score dividing by 25.52
#
#for index in range (len(delete_minus_score)):
#    delete_minus_score[index][1]=float(delete_minus_score[index][1])/25.52
#------------------------



#user_name_list process---
only_user_name=[]

for index in range (len(user_list_dlt)):
    only_user_name.append(user_list_dlt[index][0])
#------------------------------

#extraction pare----
abstact_list=[]
extracted_list=[]

for user_index in range (len(only_user_name)):
    for index in range (len(delete_minus_score)):
        if only_user_name[user_index] == delete_minus_score[index][0] : #need to change email
                abstact_list.append(delete_minus_score[index])          
    
    abstact_list=sorted(abstact_list, key=lambda abstact_list: abstact_list[4])
    extracted_list.append(abstact_list)
    abstact_list=[]          
#--------------------
    
#size filter---------over 5 time
size_filter=[]

for index in range (len(extracted_list)):
    if len(extracted_list[index]) >=5:
        size_filter.append(extracted_list[index])
    
#-------------------  


#boundary_point------------   
all_data=[]

all_data.append([99.41,2427])
all_data.append([97.41,1069])
all_data.append([95.06,132])
all_data.append([76.49,79])
all_data.append([68.10,57])
all_data.append([63.79,46])
all_data.append([54.94,39])
all_data.append([49.76,37])
all_data.append([33.00,32])
all_data.append([9.87,28])
all_data.append([3.61,27])
all_data.append([1.84,27])

#all_data0 = [99.41,2427]
#all_data1 = [97.41,1069]
#all_data2 = [95.06,132]
#all_data3 = [76.49,79]
#all_data4 = [68.10,57]
#all_data5 = [63.79,46]
#all_data6 = [54.94,39]
#all_data7 = [49.76,37]
#all_data8 = [33.00,32]
#all_data9 = [9.87,28]
#all_data10 = [3.61,27]
#all_data11 = [1.84,27]    

#--------------------------

#only 5 palying result--------------

only_five=[]

for index in range (len(size_filter)):
    length_end=len(size_filter[index])
    lenght_start= length_end - 5
    for fillter_index in range(lenght_start,length_end):
        only_five.append(size_filter[index][fillter_index])

#--------create excel file-----

import pyexcel as pe

sheet = pe.Sheet(only_five)
sheet.save_as("one.csv")

#---------------------------------
result_list=[]

for index_five in range (len(only_five)):
    compare_boundary=[]
    distance=0
    minimum_dist=0
    finding_num=0
    short_score=0
    short_time=0 
    
    for index_data in range (len(all_data)):
        distance=(float(only_five[index_five][1])-float(all_data[index_data][0]))**2+(float(only_five[index_five][2])-float(all_data[index_data][1]))**2 
        compare_boundary.append(distance)
    
    minimum_dist=min(compare_boundary)
    for index_check in range (len(compare_boundary)):
        if minimum_dist == compare_boundary[index_check]:
            finding_num=index_check
    
    short_score=abs(float(only_five[index_five][1])-float(all_data[finding_num][0]))
    short_time=abs(float(only_five[index_five][2])-float(all_data[finding_num][1]))       
    
    result_list.append(only_five[index_five]+[short_score]+[short_time])       

sheet = pe.Sheet(result_list)
sheet.save_as("result.csv")

#count game play time-----

count_list=[]

for index in range (len(size_filter)):
    count_list.append([size_filter[index][0][0]]+[len(size_filter[index])])
    
sheet = pe.Sheet(count_list)
sheet.save_as("count_list.csv")

#get average------------

average_list=[]

for index in range (len(count_list)):
    index_real=5*index
    average_score=0
    average_time=0
    
    for index_average in range(index_real,index_real+5):
        average_score=average_score+result_list[index_average][6]
        average_time=average_time+result_list[index_average][7]
    
    average_score=average_score/5
    average_time=average_time/5
    
    average_list.append(count_list[index]+[average_score]+[average_time])    
    
sheet = pe.Sheet(average_list)
sheet.save_as("average_list.csv")