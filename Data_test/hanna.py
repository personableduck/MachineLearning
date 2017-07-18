import csv
import matplotlib.pyplot as plt  
import matplotlib

with open('C:/Users/1234/Downloads/Individual_Play_Data_KYWA.csv', 'rb') as f:
    reader = csv.reader(f)
    your_list = list(reader)
    
abstact_list=[]

for index in range (len(your_list)):
    if your_list[index][0] == 'fluorine1004@naver.com': #need to change email
            abstact_list.append(your_list[index])

picture_mat_x=[]

for index2 in range (len(abstact_list)):  
    picture_mat_x.append(abstact_list[index2][1])  #change number [1]=Score,[2]=palytime,[3]=gmaeId,[4]=gamestart,[5]=gameend

picture_mat_y=[]

for index2 in range (len(abstact_list)):  
    picture_mat_y.append(abstact_list[index2][2])  

plt.figure(1)   
plt.plot(picture_mat_x,picture_mat_y,'go-')  
plt.xlabel('Score')
plt.ylabel('PlayTime')
plt.title('User:fluorine1004@naver.com')
plt.show()

new_order=sorted(abstact_list, key=lambda abstact_list: abstact_list[4])

picture_mat_x=[]
count=0

for index2 in range (len(new_order)): 
    count=count+1
    picture_mat_x.append(count)

picture_mat_y=[]

for index2 in range (len(new_order)):  
    picture_mat_y.append(new_order[index2][2])  

plt.figure(2) 
plt.plot(picture_mat_x,picture_mat_y)  
plt.xlabel('Timepass')
plt.ylabel('PlayTime')
plt.title('User:fluorine1004@naver.com')
plt.show()

new_order=sorted(abstact_list, key=lambda abstact_list: abstact_list[4])

picture_mat_x=[]
count=0

for index2 in range (len(new_order)): 
    count=count+1
    picture_mat_x.append(count)

picture_mat_y=[]

for index2 in range (len(new_order)):  
    picture_mat_y.append(new_order[index2][1])  

plt.figure(3) 
plt.plot(picture_mat_x,picture_mat_y)  
plt.xlabel('Timepass')
plt.ylabel('Score')
plt.title('User:fluorine1004@naver.com')
plt.show()

dived1=[]

for index in range (len(picture_mat_x)):
    if picture_mat_x[index] >= 50:
        dived1.append(picture_mat_x[index])
  

picture_mat_x=[]

for index2 in range (len(new_order)):  
    picture_mat_x.append(new_order[index2][1])  #change number [1]=Score,[2]=palytime,[3]=gmaeId,[4]=gamestart,[5]=gameend  
    
    
plt.figure(4)   
c=plt.scatter(picture_mat_x,picture_mat_y,c=picture_mat_x,vmin=-40,vmax=100) 

plt.colorbar(c,orientation='horizontal') 
