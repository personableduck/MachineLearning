%****************************************skip frame***************************

%distance_threshold_skip=6;
distance_threshold_skip=(40/fps)/pixelsize;

max_frame_size=fix((iter-4)/2);

skip_frame_nb=2*max_frame_size+1

new_number_count=object_number_count;

for rest_start_searching= 1:max_frame_size

    fast_obj_start=2*rest_start_searching-1;

    if(fast_obj_start == 1)
        maxsize_fast_xyz=size(xyz_11); %read size
        change_inside_loop = maxsize_fast_xyz(1,1);
        start_number_skip=1;
    else
        change_inside_loop = new_number_count; 
        start_number_skip = object_number_count;
    end
    
    for inside_loop=start_number_skip:change_inside_loop
        stop_search=0;
        distance_fast_mat=0;    
        maxsize_fast2_xyz=size(xyz{1,fast_obj_start+2}); %read size    
        
        for fast_second_mt=1:maxsize_fast2_xyz(1,1) %for fast sperm
        
            if(fast_obj_start == 1)
                
                xx_1 = xyz{1,fast_obj_start+2}(fast_second_mt,1) - xyz_11(inside_loop,1); % comparing x coordinates   
                yy_1 = xyz{1,fast_obj_start+2}(fast_second_mt,2) - xyz_11(inside_loop,2); % comparing y coordinates
                
                distance_fast_mat(fast_second_mt,1) = sqrt((xx_1^2)+(yy_1^2));
            else
                estm_ft_size=size(fast_total_coor{1,inside_loop});
                if( estm_ft_size(1,1) >= fast_obj_start-1 )
                    
                    xx_1 = xyz{1,fast_obj_start+2}(fast_second_mt,1) - fast_total_coor{1,inside_loop}(fast_obj_start,1); % comparing x coordinates   
                    yy_1 = xyz{1,fast_obj_start+2}(fast_second_mt,2) - fast_total_coor{1,inside_loop}(fast_obj_start,2); % comparing y coordinates
                    
                else     
                    xx_1=nan;
                    yy_1=nan;
                end
                
                distance_fast_mat(fast_second_mt,1) = sqrt((xx_1^2)+(yy_1^2));
            end
        end

        [row_dist,col_dist,check_dist]=find(distance_fast_mat <= distance_threshold_skip);     
    
        if (check_dist == 1) %first distance threshold check
            row_dist_size=size(row_dist); %read size  
            angle_fast=0;
            dist_thresh_sec=0;
            priority_angle_order=0;
            
            if(fast_obj_start == 1)  
                diff_vector=0;
                posit_nega=nan;
                vector=nan;
                for try_loop = 1: row_dist_size(1,1)
                    priority_angle_order(try_loop,1)=try_loop;
                    priority_angle_order(try_loop,2)=try_loop;%give priority for smaller angle!!
                end
            else 
                for try_loop = 1: row_dist_size(1,1)
                    angle_fast_test = atan((xyz{1,fast_obj_start+2}(row_dist(try_loop),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))/(xyz{1,fast_obj_start+2}(row_dist(try_loop),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi; % unit is degree   
                    %check
                    priority_angle_order(try_loop,1)=abs(cont_angle_th(inside_loop,1)-angle_fast_test);
                    priority_angle_order(try_loop,2)=try_loop;%give priority for smaller angle!!
                    
                    if(abs(cont_angle_th(inside_loop,1))>=75)
                        priority_angle_order(try_loop,1)=abs(abs(cont_angle_th(inside_loop,1))-abs(angle_fast_test));
                    end

                    [temp_pr1,ord_pr1] = sort(priority_angle_order(:,1));
                    priority_angle_order = priority_angle_order(ord_pr1,:);
                end
            end
            
            for try_loop2 = 1: row_dist_size(1,1)
       
              if(stop_search == 0) 
                pass_numfil=0;
                
                if(fast_obj_start == 1)
                    %decide vector andgle
                    %each magnitude
                    angle_fast = atan((xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2)-xyz_11(inside_loop,2))/(xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)-xyz_11(inside_loop,1)))*180/pi; % unit is degree
                    dist_thresh_sec = distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)));
                
                    if(angle_fast <= 45 & angle_fast >= -45) %vetorc
                        vector=1; % x vetorc
                        diff_vector = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1) - xyz_11(inside_loop,1);
                        posit_nega=sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                    else
                        vector=0; % y vetorc
                        diff_vector = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2) - xyz_11(inside_loop,2); 
                        posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                    end
                    %first xyz{1,fast_obj_start}(inside_loop,1), xyz{1,fast_obj_start}(inside_loop,2)      
                    %second xyz(row_dist(try_loop),1),xyz(row_dist(try_loop),2)
                    
                    posit_nega2 = posit_nega;
                    pass_numfil=1;
                
                else

                    if(cont_angle_th(inside_loop,1) <= 45 & cont_angle_th(inside_loop,1) >= -45) %vetorc
                        vector=1; % x vetorc
                        diff_vector = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1) - fast_total_coor{1,inside_loop}(fast_obj_start,1);
                        posit_nega=sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                    else 
                        vector=0; % y vetorc
                        diff_vector = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2) - fast_total_coor{1,inside_loop}(fast_obj_start,2); 
                        posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                    end

                    diff_vector2=0;
                    posit_nega2=nan;  

                    if(vector == 1)         
                        diff_vector2= xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1);    
                        posit_nega2=sign(diff_vector2);
                    elseif (vector == 0)            
                        diff_vector2= xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2);
                        posit_nega2=sign(diff_vector2);          
                    end
                    
                    angle_fast_cmp=atan((xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))/(xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi;

                    if(abs(cont_angle_th(inside_loop,1))+control_vecter_angle*multilple_short_angle > abs(angle_fast_cmp) & abs(cont_angle_th(inside_loop,1))-control_vecter_angle*multilple_short_angle < abs(angle_fast_cmp)) % unit is degree
                        pass_numfil=1;
                    end                   
                end
                
                if(posit_nega == posit_nega2) 
                    pass_value=1;   
                    if(fast_obj_start == 1)                                
                        pass_value=0;   
                    else 
                        if((distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)))+cont_rest_th(inside_loop,1))*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start-2,1))^2+(xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start-2,2))^2 ))    
                            pass_value=0; 
                        else 
                            pass_value=1; 
                        end 
                    end
                    if(pass_value == 0) 
                        if(pass_numfil == 1) 
                            if(fast_obj_start > 1)      
                                angle_fast = atan((xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))/(xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi;   
                                dist_thresh_sec = distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)));                                
                                % finish third check
                            end
                            distance_repeat_rest=0; 
                            maxsize_fast_xyz2=size(xyz{1,fast_obj_start+4}); %read size 
                            for dist_cmp_loop = 1: maxsize_fast_xyz2(1,1)
                                xx_2 = xyz{1,fast_obj_start+4}(dist_cmp_loop,1) - xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1);                               
                                yy_2 = xyz{1,fast_obj_start+4}(dist_cmp_loop,2) - xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2);
                                                              
                                distance_repeat_rest(dist_cmp_loop,1) = sqrt((xx_2^2)+(yy_2^2)); 
                            end
                            [row_secd,col_secd,check_secd]=find(distance_repeat_rest < distance_threshold_skip);
                                
                            if (check_secd == 1)
                                row_secd_size=size(row_secd); %read size   
                                angle_secd_fast=0; 
                                dist_thresh_third=0;  
                                priority_angle_order2=0;
     
                                for find_secd_loop = 1: row_secd_size(1,1) 
                                    angle_secd_fast_test = atan((xyz{1,fast_obj_start+4}(row_secd(find_secd_loop),2)-xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2))/(xyz{1,fast_obj_start+4}(row_secd(find_secd_loop),1)-xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)))*180/pi;
   
                                    priority_angle_order2(find_secd_loop,1)=abs(angle_fast-angle_secd_fast_test);
                                    priority_angle_order2(find_secd_loop,2)=find_secd_loop;%give priority for smaller angle!!
                                         
                                    if(abs(angle_fast)>=75)    
                                        priority_angle_order2(find_secd_loop,1)=abs(abs(angle_fast)-abs(angle_secd_fast_test)); 
                                    end
                                    [temp_pr1,ord_pr1] = sort(priority_angle_order2(:,1));  
                                    priority_angle_order2 = priority_angle_order2(ord_pr1,:);
                                end

                                for find_secd_loop2 = 1: row_secd_size(1,1)  
                                    if(stop_search == 0)                                       
                                        diff_vector2=0;
                                        posit_nega2=nan;
                                        pass_numfil=0;
                                        
                                        if(vector == 1)
                                            diff_vector2= xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1);
                                            posit_nega2=sign(diff_vector2);
                                        elseif (vector == 0)
                                            diff_vector2= xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2);                            
                                            posit_nega2=sign(diff_vector2);
                                        end
 
                                        if(posit_nega == posit_nega2)
                                            pass_value=1;
                                            
                                            if(fast_obj_start == 1)
                                                if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz_11(inside_loop,1))^2+(xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz_11(inside_loop,2))^2 ))
                                                    pass_value=0;
                                                else
                                                    pass_value=1;
                                                end
                                                
                                            else
                                                if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1))^2+(xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))^2 ))
                                                    pass_value=0;
                                                else
                                                    pass_value=1;
                                                end
                                            end

                                            if(pass_value == 0)                                       
                                                angle_secd_fast_cmp=atan((xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2))/(xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)))*180/pi; 
               
                                                if(abs(angle_fast)+control_vecter_angle*multilple_short_angle > abs(angle_secd_fast_cmp) & abs(angle_fast)-control_vecter_angle*multilple_short_angle < abs(angle_secd_fast_cmp)) % unit is degree                       
                                                    pass_numfil=1;
                                                end

                                                if(pass_numfil == 1) 
                                                    angle_secd_fast = atan((xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2))/(xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)))*180/pi;   
                                                    dist_thresh_third = distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)));                           
                                                    % finish third
   
                                                    if(fast_obj_start == 1)
    
                                                        new_number_count = new_number_count+1;
                                                        fst_mv_oj_coor=0;
                                         
                                                        fst_mv_oj_coor(1,1) = xyz_11(inside_loop,1);                              
                                                        fst_mv_oj_coor(1,2) = xyz_11(inside_loop,2);          
                                                        fst_mv_oj_coor(1,3) = new_number_count;
                                       
                                                        fst_mv_oj_coor(2,1) = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1);
                                                        fst_mv_oj_coor(2,2) = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2);           
                                                        fst_mv_oj_coor(2,3) = new_number_count;
     
                                                        fst_mv_oj_coor(3,1) = xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),1); 
                                                        fst_mv_oj_coor(3,2) = xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),2);                        
                                                        fst_mv_oj_coor(3,3) = new_number_count;
                    
                                                        angle_information(new_number_count,1)=abs(abs(angle_fast-angle_secd_fast));
                                          
                                                        if( fast_obj_start == 1)                                                                        
                                                            fast_total_coor(1,new_number_count)={fst_mv_oj_coor};                                                           
                                                            cont_rest_th(new_number_count,1)=dist_thresh_sec;                                      
                                                            cont_angle_th(new_number_count,1)=angle_fast;                                                
                                                            stop_search=1; %prevent duplicate trajectories                             
                                                        end
                                                        
                                                        %------prevent duplication--------                   
                                                        maxsize_fstc=size(fast_total_coor);                                        
                                                        if(fast_obj_start > 1)
                                                        maxsize_cmp_fst=size(fst_mv_oj_coor);               
                                                        duplication_number=0;                                                   
                                                        for prevent_duplication_loop=1:maxsize_fstc(1,2)	
                                                            duplication_happen=0;
                                                            for fst_check_loop =1:maxsize_cmp_fst(1,1)
                                                                [row_dplc_x,col_dplc_x,check_dplc_x]=find(fast_total_coor{1,prevent_duplication_loop}(:,1) == fst_mv_oj_coor(fst_check_loop,1));
                                                                [row_dplc_y,col_dplc_y,check_dplc_y]=find(fast_total_coor{1,prevent_duplication_loop}(:,2) == fst_mv_oj_coor(fst_check_loop,2));                                                 
                                                                rr_dplc=size(row_dplc_x); 
                                                                for dplc_chk_loop = 1:rr_dplc                                                                                                           
                                                                    [row_find_dplc,col_find_dplc,check_find_dplc]= find(row_dplc_y == row_dplc_x(dplc_chk_loop));
                                                                    if(check_find_dplc ~= 0)              
                                                                        duplication_happen=duplication_happen+1;  
                                                                    end                                                
                                                                end
                                                            end
                                                            if(duplication_happen >= 2)        
                                                                duplication_number=prevent_duplication_loop;
                                                            end
                                                        end

                                                        if(duplication_number ~= 0)                             
                                                            if(angle_information(duplication_number,1) < angle_information(new_number_count,1))                              
                                                                new_number_count = new_number_count-1; %duplication happen                         
                                                            else
                                                                fast_total_coor{1,duplication_number} = {0};
                                                                fast_total_coor{1,duplication_number} = fst_mv_oj_coor; %smaller angle is right trajectory         
                                                                new_number_count = new_number_count-1;                        
                                                            end
                                                        elseif(duplication_number == 0)                           
                                                            fast_total_coor(1,new_number_count)={fst_mv_oj_coor};   
                                                            cont_rest_th(new_number_count,1)=dist_thresh_sec;                        
                                                            cont_angle_th(new_number_count,1)=angle_fast;                    
                                                            stop_search=1; %prevent duplicate trajectories                              
                                                        end
                                                        end
                                                        
                                                    else
                                                        
                                                        fast_total_coor{1,inside_loop}(rest_start_searching+1,1) = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1);
                                                        fast_total_coor{1,inside_loop}(rest_start_searching+1,2) = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2);                                                      
                                                        fast_total_coor{1,inside_loop}(rest_start_searching+1,3) = inside_loop;
                                                        
                                                        cont_rest_th(inside_loop,1)=dist_thresh_sec;                                                          
                                                        cont_angle_th(inside_loop,1)=angle_fast;                        
                                                        stop_search=1; %prevent duplicate trajectories                                         
                                                    end              
                                                end   
                                            end
                                        end  
                                    end 
                                end 
                            end   
                        end
                    end  
                end
              end
            end
        end
    end  
end
                                       