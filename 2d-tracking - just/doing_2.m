
fps=3; %fps means frame per second   

for divided_frame = 1:(iter-1) %dividing fram number 
    
maxsize_num_xyz=size(xyz{1,divided_frame}); %read size
read_size_xyz=maxsize_num_xyz(1,1);

xyzs=xyz{1,divided_frame}; %copy elements of xyz
xyzs(:,3)=0;

last_distance_comp={0}; %initialization matrix
count_fair=0; %initialization value
prevent_dup=0; %initialization value


for loop_xyz=1:(read_size_xyz)-1

        if( find(prevent_dup == loop_xyz)) % eliminate duplication       
            continue
            
        end
            
distance_threshold = (110/fps)/pixelsize; %******threshold object tracking. need to control. maximum ditance to distinguish as a same object
%unit is um. // I limited that the fasted sperm's speed is 110 um/sec.
%fps means frame per second
count_distance_num=0;
distance_mat=0;

for realative_loop=(loop_xyz+1):read_size_xyz
      
    count_distance_num=count_distance_num+1;
    
    xx = xyz{1,divided_frame}(realative_loop,1) - xyz{1,divided_frame}(loop_xyz,1); % comparing x coordinates
    yy = xyz{1,divided_frame}(realative_loop,2) - xyz{1,divided_frame}(loop_xyz,2); % comparing y coordinates
    
    distance_mat(1,count_distance_num) = sqrt((xx^2)+(yy^2));

end
    
    last_distance_comp(1,loop_xyz)={distance_mat};

    [mini_distance,find_order] = min(last_distance_comp{1,loop_xyz}); %find minimum distance
    
    if( mini_distance <= distance_threshold) % finding moving objects, match function
        
        if (min_intensities_xyz{1,divided_frame}(loop_xyz) < min_intensities_xyz{1,divided_frame}(loop_xyz+find_order)) % compare intensity to decide before one. making oreder between 2 frames
        
            xyzs(loop_xyz,3)= (divided_frame*2);
            xyzs((loop_xyz+find_order),3)= (divided_frame*2)-1; 
        
            count_fair=count_fair+1;
            prevent_dup(count_fair,1) = (loop_xyz+find_order);
        
        else       
            xyzs(loop_xyz,3)= (divided_frame*2)-1;
            xyzs((loop_xyz+find_order),3)= (divided_frame*2);
        
        
            count_fair=count_fair+1;
            prevent_dup(count_fair,1) = (loop_xyz+find_order);
            
        end
  
    end
    
end


%------ slice devide by frame

finding_frame_num=max(xyzs(:,3)); %read devided frame number

size_func_xyzs=size(xyzs);
read_size_xyzs=size_func_xyzs(1,1);


for devided_frame_num_loop= finding_frame_num-1: finding_frame_num
    
    count_increase=0; %increasing number to count number

    for put_xyzs_num= 1:read_size_xyzs
        
        if xyzs(put_xyzs_num,3) == devided_frame_num_loop
        count_increase=count_increase+1;
        devided_coordinates{1,devided_frame_num_loop}(count_increase,:) = xyzs(put_xyzs_num,:);
        
        end
        
        
        if  xyzs(put_xyzs_num,3) == 0
        count_increase=count_increase+1;
        devided_coordinates{1,devided_frame_num_loop}(count_increase,:) = xyzs(put_xyzs_num,:);
        devided_coordinates{1,devided_frame_num_loop}(count_increase,3) = devided_frame_num_loop;

        end
            
    end
    

end

end



%---------- combine each frame matrix 
sum_dev_coordinate = devided_coordinates{1,1};

for sumloop_mat = 2:(iter-1)*2
sum_dev_coordinate=[sum_dev_coordinate; devided_coordinates{1,sumloop_mat}];
end

%---------------

trajectories_imfomation=track(sum_dev_coordinate,distance_threshold); % use IDL track funtion, need to change figure

[temp,ord] = sort(trajectories_imfomation(:,1));
xyzs_ord = trajectories_imfomation(ord,:);

%----- devide again by object

object_number_tj=max(trajectories_imfomation(:,4));

size_obj_num=size(trajectories_imfomation);
real_size_objnum=size_obj_num(1,1);



for loop_objnum=1:object_number_tj
    
    count_inc_obj=0;
    
for count_moving_loop=1:real_size_objnum
    
    if trajectories_imfomation(count_moving_loop,4) == loop_objnum
        count_inc_obj=count_inc_obj+1;
        object_trajectories{1,loop_objnum}(count_inc_obj,:) = trajectories_imfomation(count_moving_loop,:);
    end
    
end   
end   

figure, imshow(abs(recon_stack(:,:,1)),[])
hold on

for ssizz = 1:object_number_tj
plot((object_trajectories{1,ssizz}(:,1)),(object_trajectories{1,ssizz}(:,2)),'g.')
hold on
plot((object_trajectories{1,ssizz}(1,1)),(object_trajectories{1,ssizz}(1,2)),'b.')
end

% figure, imshow(recon_stack_diff(:,:,1),[])
% hold on
% plot(devided_coordinates{1,1}(:,1),devided_coordinates{1,1}(:,2),'g.')
% hold on
% plot(devided_coordinates{1,2}(:,1),devided_coordinates{1,2}(:,2),'b.')

% [temp,ord] = sort(xyzs(:,3));
% xyzs_ord = xyzs(ord,:);