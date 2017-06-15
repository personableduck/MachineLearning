coor_from{1,iter-1} = copy_cf{1,iter-1};
coor_to{1,iter-1} = copy_ct{1,iter-1};

number_coor_to_matrix= size(coor_to{1,iter-1});
number_coor_to= number_coor_to_matrix(1,1);

for routine_coor_to= 1:number_coor_to

obj_coor((iter-1)*2,2)=0; %initialization // iter is image number
    
obj_coor(1,1) = coor_to{1,iter-1}(routine_coor_to,1); %input first and second figures from coor_from matrix
obj_coor(1,2) = coor_to{1,iter-1}(routine_coor_to,2);
obj_coor(2,1) = coor_from{1,iter-1}(routine_coor_to,1);
obj_coor(2,2) = coor_from{1,iter-1}(routine_coor_to,2);

object_matrix_number=0;

for number_routine_obj_coor= 2:(iter-1)

object_matrix_number=object_matrix_number+1;

number_coor_to_matrix_change= size(coor_to{1,iter-number_routine_obj_coor});
number_coor_to_matrix_change_value= number_coor_to_matrix_change(1,1);
    
for routine_number_for_coor_to= 1:number_coor_to_matrix_change_value
    
if((negative_threshold < (obj_coor((object_matrix_number*2),1) - coor_to{1,iter-number_routine_obj_coor}(routine_number_for_coor_to,1)) & (obj_coor((object_matrix_number*2),1) - coor_to{1,iter-number_routine_obj_coor}(routine_number_for_coor_to,1)) < positive_threshold) & (negative_threshold < (obj_coor((object_matrix_number*2),2) - coor_to{1,iter-number_routine_obj_coor}(routine_number_for_coor_to,2)) & (obj_coor((object_matrix_number*2),2) - coor_to{1,iter-number_routine_obj_coor}(routine_number_for_coor_to,2)) < positive_threshold))
    
    obj_coor((2*number_routine_obj_coor)-1,1) = coor_to{1,iter-number_routine_obj_coor}(routine_number_for_coor_to,1);
    obj_coor((2*number_routine_obj_coor)-1,2) = coor_to{1,iter-number_routine_obj_coor}(routine_number_for_coor_to,2);
    obj_coor((2*number_routine_obj_coor),1) = coor_from{1,iter-number_routine_obj_coor}(routine_number_for_coor_to,1);
    obj_coor((2*number_routine_obj_coor),2) = coor_from{1,iter-number_routine_obj_coor}(routine_number_for_coor_to,2);
    
 
end 

for detect= 1: (iter-1)*2 % if object don't move between frame, we have to show don't moving.This code is to make up for that
    
    if obj_coor(detect,1)==0 & obj_coor(detect,2) ==0
        obj_coor(detect,1)=obj_coor((detect-1),1);
        obj_coor(detect,2)=obj_coor((detect-1),2);
    end
end

end
end




obj_mat_cor(routine_coor_to,1)={obj_coor};
obj_coor=0;

end


%%%%%%SPEED

disfps= 10; %fps
timess= 9141/1000;  %fps*iter;

dkcor= size(obj_mat_cor);
dkcormax= dkcor(1,1);

spsp=0;

for ababs = 1: dkcormax
spsp=spsp+1;
sumdis=0;

for objds = 1:((iter-1)*2-1)
    
    obdistances= pixelsize*(sqrt((obj_mat_cor{ababs,1}(objds+1,1)-obj_mat_cor{ababs,1}(objds,1))^2 + (obj_mat_cor{ababs,1}(objds+1,2)-obj_mat_cor{ababs,1}(objds,2))^2)); % um
    
    sumdis = obdistances + sumdis;
    
end

speedss = sumdis/timess; %um/s

obj_spd(spsp,1) = speedss;
end

hist(obj_spd,25)
xlabel('Speed(um/sec)')
ylabel('Counts')
median(obj_spd)
mean(obj_spd)
dkcormax %objectnumber