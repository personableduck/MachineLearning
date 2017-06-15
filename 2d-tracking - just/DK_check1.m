object_number_count=0;

max_frame_size=fix((iter-1)/3)-1;
%295
for bundle_frame=0:max_frame_size

fast_obj_start=4*bundle_frame+1; % if you need to find starting point on the middle of frame, then you need to change the value as a loop(for function)
maxsize_fast_xyz=size(xyz{1,fast_obj_start}); %read size  

if(fast_obj_start>1)
    for rest_tj_loop=1:object_number_count
    object_number_count = object_number_count+1;
    end
end


%293
for inside_loop=1:maxsize_fast_xyz(1,1)-1

    count_distance_num=0;
    distance_fast_mat=0;
        
    if(fast_obj_start>1)
        for inside2_loop=1:maxsize_fast_xyz(1,1)
        count_distance_num=count_distance_num+1;
        
        xx_1 = xyz{1,fast_obj_start}(inside2_loop,1) - fast_total_coor{1,rest_tj_loop}(last,1); % comparing x coordinates   
        yy_1 = xyz{1,fast_obj_start}(inside2_loop,2) - fast_total_coor{1,rest_tj_loop}(last,2); % comparing y coordinates 
        
        distance_fast_mat(count_distance_num,1) = sqrt((xx_1^2)+(yy_1^2));
        end
    end
       
    for inside_bhd=inside_loop+1:maxsize_fast_xyz(1,1)
        
        if(fast_obj_start>1)
            continue
        end

        count_distance_num=count_distance_num+1;
           
        xx_1 = xyz{1,fast_obj_start}(inside_bhd,1) - xyz{1,fast_obj_start}(inside_loop,1); % comparing x coordinates   
        yy_1 = xyz{1,fast_obj_start}(inside_bhd,2) - xyz{1,fast_obj_start}(inside_loop,2); % comparing y coordinates
             
        distance_fast_mat(count_distance_num,1) = sqrt((xx_1^2)+(yy_1^2));
            
    end
  
    [row_dist,col_dist,check_dist]=find(distance_fast_mat <= distance_threshold_fast);
        
    if (check_dist == 1) %first distance threshold check
        
        row_dist_size=size(row_dist); %read size  
        angle_fast=0;
        dist_thresh_sec=0;
      
        for try_loop = 1: row_dist_size(1,1)
            

            if (min_intensities_xyz{1,fast_obj_start}(inside_loop) > min_intensities_xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop))) % compare intensity to decide before one. smaller intensity means later one       
                %decide vector andgle
                %each magnitude
                angle_fast = (atan((xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),2)-xyz{1,fast_obj_start}(inside_loop,2))/(xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),1)-xyz{1,fast_obj_start}(inside_loop,1))))*180/pi; % unit is degree
                dist_thresh_sec = distance_fast_mat(row_dist(try_loop));
                
                %first xyz{1,fast_obj_start}(inside_loop,1), xyz{1,fast_obj_start}(inside_loop,2)
                %second xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),1),xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),2)
                
                sec_dist_size=size(xyz{1,fast_obj_start+1}); %read size 
                
                distance_secd_mat=0;
                
                for dist_cmp_loop = 1: sec_dist_size(1,1)
                    
                    
                    xx_2 = xyz{1,fast_obj_start+1}(dist_cmp_loop,1) - xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),1);
                    yy_2 = xyz{1,fast_obj_start+1}(dist_cmp_loop,2) - xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),2);
                
                    distance_secd_mat(dist_cmp_loop,1) = sqrt((xx_2^2)+(yy_2^2));
                    
                end
                
                [row_secd,col_secd,check_secd]=find(distance_secd_mat < dist_thresh_sec * (1.3) & distance_secd_mat > dist_thresh_sec * (0.7) & distance_secd_mat ~= dist_thresh_sec);
                
                if (check_secd == 1)
              
                    row_secd_size=size(row_secd); %read size  
                    angle_secd_fast=0;
                    dist_thresh_third=0;
                   
                    for find_secd_loop = 1: row_secd_size(1,1)
                        if(dist_thresh_sec*1.3 < sqrt( (xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1)-xyz{1,fast_obj_start}(inside_loop,1))^2+(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2)-xyz{1,fast_obj_start}(inside_loop,2))^2 ))
                        if(angle_fast+control_vecter_angle > atan((xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2)-xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),2))/((xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1))-xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),1)))*180/pi & angle_fast-control_vecter_angle < atan((xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2)-xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),2))/((xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1))-xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),1)))*180/pi) % unit is degree
             
                            angle_secd_fast = atan((xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2)-xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),2))/((xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1))-xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),1)))*180/pi;
                            dist_thresh_third = distance_secd_mat(row_secd(find_secd_loop));
                            % finish third check
                            % (xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1)),(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2))
                            
                            thrd_dist_size=size(xyz{1,fast_obj_start+2}); %read size 
                            distance_thrd_mat=0;
                            
                            for dist_thrd_loop = 1: thrd_dist_size(1,1)
                                
                                xx_3 = xyz{1,fast_obj_start+2}(dist_thrd_loop,1) - (xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1));
                                yy_3 = xyz{1,fast_obj_start+2}(dist_thrd_loop,2) - (xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2));

                                distance_thrd_mat(dist_thrd_loop,1) = sqrt((xx_3^2)+(yy_3^2));
                                
                            end
                                          
                            [row_third,col_third,check_third]=find(distance_thrd_mat <= dist_thresh_third * (1.3) & distance_thrd_mat >= dist_thresh_third * (0.7));
                
                            if (check_third == 1)
                               
                                row_third_size=size(row_third);
                                angle_rest_th=0;
                                dist_rest_th=0;
                                
                                
                                           
                                for find_rest_loop = 1: row_third_size(1,1)

                                    if(dist_thresh_third*1.3 < sqrt((xyz{1,fast_obj_start+2}(row_third(find_rest_loop),1)- xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),1))^2+(xyz{1,fast_obj_start+2}(row_third(find_rest_loop),2)- xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),2))^2))
                                    if(angle_secd_fast+control_vecter_angle > atan((xyz{1,fast_obj_start+2}(row_third(find_rest_loop),2)-(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2)))/(xyz{1,fast_obj_start+2}(row_third(find_rest_loop),1)-(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1))))*180/pi & angle_secd_fast-control_vecter_angle < atan((xyz{1,fast_obj_start+2}(row_third(find_rest_loop),2)-(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2)))/(xyz{1,fast_obj_start+2}(row_third(find_rest_loop),1)-(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1))))*180/pi)     
   
                                        angle_rest_th=atan((xyz{1,fast_obj_start+2}(row_third(find_rest_loop),2)-(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2)))/(xyz{1,fast_obj_start+2}(row_third(find_rest_loop),1)-(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1))))*180/pi;
                                        dist_rest_th=distance_thrd_mat(row_third(find_rest_loop));
                                        
                                        %fourth xyz{1,fast_obj_start+2}(row_third(find_third_loop),1),xyz{1,fast_obj_start+2}(row_third(find_third_loop),2)
                                        if(fast_obj_start==1)
                                        object_number_count = object_number_count+1;
                                        end
                                        

                                        fst_mv_oj_coor(fast_obj_start,1) = xyz{1,fast_obj_start}(inside_loop,1);                              
                                        fst_mv_oj_coor(fast_obj_start,2) = xyz{1,fast_obj_start}(inside_loop,2);          
                                        fst_mv_oj_coor(fast_obj_start,3) = 2; %fast moving value is 2
                                        fst_mv_oj_coor(fast_obj_start,4) = object_number_count;
                                   
                                        fst_mv_oj_coor(fast_obj_start+1,1) = xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),1);
                                        fst_mv_oj_coor(fast_obj_start+1,2) = xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),2);          
                                        fst_mv_oj_coor(fast_obj_start+1,3) = 2; %fast moving value is 2
                                        fst_mv_oj_coor(fast_obj_start+1,4) = object_number_count;
                                   
                                        fst_mv_oj_coor(fast_obj_start+2,1) = xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1);
                                        fst_mv_oj_coor(fast_obj_start+2,2) = xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2);          
                                        fst_mv_oj_coor(fast_obj_start+2,3) = 2; %fast moving value is 2
                                        fst_mv_oj_coor(fast_obj_start+2,4) = object_number_count;
                                        
                                        fst_mv_oj_coor(fast_obj_start+3,1) = xyz{1,fast_obj_start+2}(row_third(find_rest_loop),1);
                                        fst_mv_oj_coor(fast_obj_start+3,2) = xyz{1,fast_obj_start+2}(row_third(find_rest_loop),2);          
                                        fst_mv_oj_coor(fast_obj_start+3,3) = 2; %fast moving value is 2
                                        fst_mv_oj_coor(fast_obj_start+3,4) = object_number_count;
                                        
                                        last=fast_obj_start+3;

                                        fast_total_coor(1,object_number_count)={fst_mv_oj_coor};
                                       

                                    end
                                    end
                                end
                                
                            end
                            
                        end
                        end   
                    end
        
                end

                
            else
                %need to change order for first
                
                angle_fast = (atan((xyz{1,fast_obj_start}(inside_loop,2)-xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),2))/(xyz{1,fast_obj_start}(inside_loop,1)-xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),1))))*180/pi; % unit is degree
                dist_thresh_sec = distance_fast_mat(row_dist(try_loop));
                  
                %first xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),1),xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),2)
                %second xyz{1,fast_obj_start}(inside_loop,1), xyz{1,fast_obj_start}(inside_loop,2)
         
                sec_dist_size=size(xyz{1,fast_obj_start+1}); %read size 
                
                distance_secd_mat=0;
                
                for dist_cmp_loop = 1: sec_dist_size(1,1)
                    
                    
                    xx_2 = xyz{1,fast_obj_start+1}(dist_cmp_loop,1) - xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),1);
                    yy_2 = xyz{1,fast_obj_start+1}(dist_cmp_loop,2) - xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),2);
                
                    distance_secd_mat(dist_cmp_loop,1) = sqrt((xx_2^2)+(yy_2^2));
                    
                end
                
                [row_secd,col_secd,check_secd]=find(distance_secd_mat < dist_thresh_sec * (1.3) & distance_secd_mat > dist_thresh_sec * (0.7) & distance_secd_mat ~= dist_thresh_sec);
                
                if (check_secd == 1)
                 
                    row_secd_size=size(row_secd); %read size  
                    angle_secd_fast=0;
                    dist_thresh_third=0;
                 
                    for find_secd_loop = 1: row_secd_size(1,1)
                        if(dist_thresh_sec*1.3 < sqrt( (xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1)-xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),1))^2+(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2)-xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),2))^2 ))
                        if(angle_fast+control_vecter_angle > atan((xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2)-xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),2))/((xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1))-xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),1)))*180/pi & angle_fast-control_vecter_angle < atan((xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2)-xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),2))/((xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1))-xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),1)))*180/pi) % unit is degree
             
                            angle_secd_fast = atan((xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2)-xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),2))/((xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1))-xyz{1,fast_obj_start}((inside_loop+row_dist(try_loop)),1)))*180/pi;
                            dist_thresh_third = distance_secd_mat(row_secd(find_secd_loop));
                            % finish third check
                            % (xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1)),(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2))
                            
                            thrd_dist_size=size(xyz{1,fast_obj_start+2}); %read size 
                            distance_thrd_mat=0;
                            
                            for dist_thrd_loop = 1: thrd_dist_size(1,1)
                                
                                xx_3 = xyz{1,fast_obj_start+2}(dist_thrd_loop,1) - (xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1));
                                yy_3 = xyz{1,fast_obj_start+2}(dist_thrd_loop,2) - (xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2));

                                distance_thrd_mat(dist_thrd_loop,1) = sqrt((xx_3^2)+(yy_3^2));
                                
                            end
                            
                            [row_third,col_third,check_third]=find(distance_thrd_mat <= dist_thresh_third * (1.3) & distance_thrd_mat >= dist_thresh_third * (0.7));
                
                            if (check_third == 1)
                               
                                row_third_size=size(row_third);
                                angle_rest_th=0;
                                dist_rest_th=0;
                                
                                
                                       
                                for find_rest_loop = 1: row_third_size(1,1)
                                   
                                    if(dist_thresh_third*1.3 < sqrt( (xyz{1,fast_obj_start+2}(row_third(find_rest_loop),1)-xyz{1,fast_obj_start}(inside_loop,1))^2+(xyz{1,fast_obj_start+2}(row_third(find_rest_loop),2)-xyz{1,fast_obj_start}(inside_loop,2))^2 ))
                                    if(angle_secd_fast+control_vecter_angle > atan((xyz{1,fast_obj_start+2}(row_third(find_rest_loop),2)-(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2)))/(xyz{1,fast_obj_start+2}(row_third(find_rest_loop),1)-(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1))))*180/pi & angle_secd_fast-control_vecter_angle < atan((xyz{1,fast_obj_start+2}(row_third(find_rest_loop),2)-(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2)))/(xyz{1,fast_obj_start+2}(row_third(find_rest_loop),1)-(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1))))*180/pi)     

                                        angle_rest_th=atan((xyz{1,fast_obj_start+2}(row_third(find_rest_loop),2)-(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2)))/(xyz{1,fast_obj_start+2}(row_third(find_rest_loop),1)-(xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1))))*180/pi;
                                        dist_rest_th=distance_thrd_mat(row_third(find_rest_loop));
                                        
                                        %fourth xyz{1,fast_obj_start+2}(row_third(find_third_loop),1),xyz{1,fast_obj_start+2}(row_third(find_third_loop),2)

                                        if(fast_obj_start==1)
                                        object_number_count = object_number_count+1;
                                        end        
                                        
                                        
                                        fst_mv_oj_coor(fast_obj_start,1) = xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),1);                              
                                        fst_mv_oj_coor(fast_obj_start,2) = xyz{1,fast_obj_start}(inside_loop+row_dist(try_loop),2);          
                                        fst_mv_oj_coor(fast_obj_start,3) = 2; %fast moving value is 2
                                        fst_mv_oj_coor(fast_obj_start,4) = object_number_count;
                                   
                                        fst_mv_oj_coor(fast_obj_start+1,1) = xyz{1,fast_obj_start}(inside_loop,1);
                                        fst_mv_oj_coor(fast_obj_start+1,2) = xyz{1,fast_obj_start}(inside_loop,2);          
                                        fst_mv_oj_coor(fast_obj_start+1,3) = 2; %fast moving value is 2
                                        fst_mv_oj_coor(fast_obj_start+1,4) = object_number_count;
                                   
                                        fst_mv_oj_coor(fast_obj_start+2,1) = xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),1);
                                        fst_mv_oj_coor(fast_obj_start+2,2) = xyz{1,fast_obj_start+1}(row_secd(find_secd_loop),2);          
                                        fst_mv_oj_coor(fast_obj_start+2,3) = 2; %fast moving value is 2
                                        fst_mv_oj_coor(fast_obj_start+2,4) = object_number_count;
                                        
                                        fst_mv_oj_coor(fast_obj_start+3,1) = xyz{1,fast_obj_start+2}(row_third(find_rest_loop),1);
                                        fst_mv_oj_coor(fast_obj_start+3,2) = xyz{1,fast_obj_start+2}(row_third(find_rest_loop),2);         
                                        fst_mv_oj_coor(fast_obj_start+3,3) = 2; %fast moving value is 2
                                        fst_mv_oj_coor(fast_obj_start+3,4) = object_number_count;
                                        
                                        last=fast_obj_start+3;
                                        
                                        fast_total_coor(1,object_number_count)={fst_mv_oj_coor};                           

                                    end
                                    end
                                    
                                end %finish first find_rest loop
                               
                            end
                            
                        end
                        end   
                    end
        
                end
                
                
            end%intensity
            end
        end
end
    
end
% 
% for a=1:3
% a
% 
% if (a==2)
% continue
% end
% for aa=1:3
% 
% aa
% end
% end