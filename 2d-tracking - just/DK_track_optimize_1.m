%****************************************main trajectory***************************
max_frame_size=fix((iter-2)/4);

frame_number=4*max_frame_size+1

object_number_count=0;
cont_rest_th=0;
cont_angle_th=0;
fast_total_coor={0};
angle_information=0;

for rest_start_searching= 1:max_frame_size

    prevent_rest_check={0};
    angle_information=0;
    fast_obj_start=4*rest_start_searching-3;
    if(fast_obj_start == 1)
        maxsize_fast_xyz=size(xyz{1,fast_obj_start}); %read size
        change_inside_loop = maxsize_fast_xyz(1,1);
    else
        change_inside_loop = object_number_count;
    end
    
    for inside_loop=1:change_inside_loop
        stop_search=0;
        distance_fast_mat=0;    
        if(fast_obj_start == 1)
            [row_f1,col_f1,check_f1]=find( xyz{1,fast_obj_start+1}(:,1) < xyz{1,fast_obj_start}(inside_loop,1) + distance_threshold_fast & xyz{1,fast_obj_start+1}(:,1) > xyz{1,fast_obj_start}(inside_loop,1) - distance_threshold_fast);   
        else
            estm_ft_size=size(fast_total_coor{1,inside_loop});
            if( estm_ft_size(1,1) >= fast_obj_start )
                [row_f1,col_f1,check_f1]=find( xyz{1,fast_obj_start+1}(:,1) < fast_total_coor{1,inside_loop}(fast_obj_start,1) + distance_threshold_fast & xyz{1,fast_obj_start+1}(:,1) > fast_total_coor{1,inside_loop}(fast_obj_start,1) - distance_threshold_fast);
            end
        end
        if(check_f1 == 1)
            maxsize_fast2_xyz=size(row_f1); %read size
            distance_count1=0;
            for fast_second_mt=1:maxsize_fast2_xyz(1,1) %for fast sperm
                if(fast_obj_start == 1)
                    xx_1 = xyz{1,fast_obj_start+1}(row_f1(fast_second_mt),1) - xyz_1{1,fast_obj_start}(inside_loop,1); % comparing x coordinates   
                    yy_1 = xyz{1,fast_obj_start+1}(row_f1(fast_second_mt),2) - xyz_1{1,fast_obj_start}(inside_loop,2); % comparing y coordinates
                    distance_count1=distance_count1+1;
                    distance_fast_mat(distance_count1,1) = sqrt((xx_1^2)+(yy_1^2));
                    distance_fast_mat(distance_count1,2) = row_f1(fast_second_mt);
                else
                    estm_ft_size=size(fast_total_coor{1,inside_loop});
                    if( estm_ft_size(1,1) >= fast_obj_start )
                        xx_1 = xyz{1,fast_obj_start+1}(row_f1(fast_second_mt),1) - fast_total_coor{1,inside_loop}(fast_obj_start,1); % comparing x coordinates   
                        yy_1 = xyz{1,fast_obj_start+1}(row_f1(fast_second_mt),2) - fast_total_coor{1,inside_loop}(fast_obj_start,2); % comparing y coordinates
                    else     
                        xx_1=nan;
                        yy_1=nan;
                    end
                    distance_count1=distance_count1+1;
                    distance_fast_mat(distance_count1,1) = sqrt((xx_1^2)+(yy_1^2));
                    distance_fast_mat(distance_count1,2) = row_f1(fast_second_mt);
                end
            end
            if(fast_obj_start == 1) 
                [row_dist,col_dist,check_dist]=find(distance_fast_mat(:,1) <= distance_threshold_fast);     
            else
                if( cont_rest_th(inside_loop,1) <= decide_slow)
                    [row_dist,col_dist,check_dist]=find(distance_fast_mat(:,1) <= distance_threshold_fast);            
                else
                    [row_dist,col_dist,check_dist]=find(distance_fast_mat(:,1) <= cont_rest_th(inside_loop,1)*threshold_diff_dislocation & distance_fast_mat(:,1) <= distance_threshold_fast);
                end
            end
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
                        priority_angle_order(try_loop,1)=distance_fast_mat(row_dist(try_loop));
                        if(priority_angle_order(try_loop,1) == 0)
                            priority_angle_order(try_loop,1) = priority_zero; 
                        end
                        priority_angle_order(try_loop,2)=try_loop;%give priority for smaller distance!!
                        [temp_pr1,ord_pr1] = sort(priority_angle_order(:,1));
                        priority_angle_order = priority_angle_order(ord_pr1,:);
                    end
                else 
                    for try_loop = 1: row_dist_size(1,1) 
                        angle_fast_test = atan((xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(try_loop),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))/(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(try_loop),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi; % unit is degree   
                        priority_angle_order(try_loop,1)=abs(abs(cont_angle_th(inside_loop,1))-abs(angle_fast_test))/priority_degree_rate + distance_fast_mat(row_dist(try_loop));                  
                        if( (xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(try_loop),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2)) == 0 & (xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(try_loop),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)) == 0 )
                            priority_angle_order(try_loop,1) = priority_zero; 
                        end
                        priority_angle_order(try_loop,2)=try_loop;%give priority for smaller angle!!
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
                            angle_fast = atan((xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-xyz{1,fast_obj_start}(inside_loop,2))/(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-xyz{1,fast_obj_start}(inside_loop,1)))*180/pi; % unit is degree
                            if( (xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-xyz{1,fast_obj_start}(inside_loop,2)) == 0 & (xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-xyz{1,fast_obj_start}(inside_loop,1)) == 0 )
                                angle_fast = 0;
                            end
                            dist_thresh_sec = distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)));
                            if(angle_fast <= 45 & angle_fast >= -45) %vetorc
                                vector=1; % x vetorc
                                diff_vector = xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1) - xyz{1,fast_obj_start}(inside_loop,1);
                                posit_nega=sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                            else
                                vector=0; % y vetorc
                                diff_vector = xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2) - xyz{1,fast_obj_start}(inside_loop,2); 
                                posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                            end
                            posit_nega2 = posit_nega;
                            pass_numfil=1;
                        else
                            if(cont_angle_th(inside_loop,1) <= 45 & cont_angle_th(inside_loop,1) >= -45) %vetorc
                                vector=1; % x vetorc
                                diff_vector = fast_total_coor{1,inside_loop}(fast_obj_start,1) - fast_total_coor{1,inside_loop}(fast_obj_start-1,1) ;
                                posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                            else 
                                vector=0; % y vetorc
                                diff_vector = fast_total_coor{1,inside_loop}(fast_obj_start,2) - fast_total_coor{1,inside_loop}(fast_obj_start-1,2) ;
                                posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                            end
                            diff_vector2=0;
                            posit_nega2=nan;
                            if(vector == 1)         
                                diff_vector2= xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1);    
                                posit_nega2=sign(diff_vector2);
                            elseif (vector == 0)            
                                diff_vector2= xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2);
                                posit_nega2=sign(diff_vector2);          
                            end                    
                            if(posit_nega == posit_nega2 | posit_nega2 == 0 | posit_nega == 0)
                                angle_fast_test=atan((xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))/(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi;
                                if( (xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2)) == 0 & (xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)) == 0)
                                    angle_fast_test = 0;
                                end
                                if( cont_rest_th(inside_loop,1) > decide_slow_angle)           
                                    if(abs(cont_angle_th(inside_loop,1)) > tangent_control | abs(angle_fast_test) > tangent_control)              
                                        if(abs(cont_angle_th(inside_loop,1))+control_vecter_angle > abs(atan((xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start-1,2))/(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start-1,1)))*180/pi) & abs(cont_angle_th(inside_loop,1))-control_vecter_angle < abs(atan((xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start-1,2))/(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start-1,1)))*180/pi)) % unit is degree
                                            pass_numfil=1;
                                        end
                                    else
                                        if(cont_angle_th(inside_loop,1)+control_vecter_angle > atan((xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start-1,2))/(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start-1,1)))*180/pi & cont_angle_th(inside_loop,1)-control_vecter_angle < atan((xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start-1,2))/(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start-1,1)))*180/pi) % unit is degree
                                        pass_numfil=1;
                                        end
                                    end    
                                else
                                    if(abs(cont_angle_th(inside_loop,1)) > tangent_control | abs(angle_fast_test) > tangent_control)              
                                        if(abs(cont_angle_th(inside_loop,1))+control_vecter_angle*multilple_short_angle > abs(atan((xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start-1,2))/(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start-1,1)))*180/pi) & abs(cont_angle_th(inside_loop,1))-control_vecter_angle*multilple_short_angle < abs(atan((xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start-1,2))/(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start-1,1)))*180/pi)) % unit is degree
                                            pass_numfil=1;
                                        end
                                    else
                                        if(cont_angle_th(inside_loop,1)+control_vecter_angle*multilple_short_angle > atan((xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start-1,2))/(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start-1,1)))*180/pi & cont_angle_th(inside_loop,1)-control_vecter_angle*multilple_short_angle < atan((xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start-1,2))/(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start-1,1)))*180/pi) % unit is degree
                                        pass_numfil=1;
                                        end
                                    end
                                end
                            end
                        end
                        increase_check=1;
                        if(fast_obj_start == 1)                         
                            increase_check=0; 
                        else
                            if((distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)))+cont_rest_th(inside_loop,1))*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start-1,1))^2+(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start-1,2))^2 ))          
                                increase_check=0; 
                            else   
                                increase_check=1;
                            end          
                        end
                        if(increase_check == 0)       
                            if(pass_numfil == 1) 
                                if(fast_obj_start > 1)
                                    angle_fast = atan((xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))/(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi;
                                    if( (xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2)) == 0 & (xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)) == 0)
                                        angle_fast = 0;
                                    end
                                    dist_thresh_sec = distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)));                                
                                end   
                                distance_repeat_rest=0;
                                [row_f2,col_f2,check_f2]=find( xyz{1,fast_obj_start+2}(:,1) < xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1) + distance_threshold_fast & xyz{1,fast_obj_start+2}(:,1) > xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1) - distance_threshold_fast);       

                                if(check_f2 == 1)
                                    maxsize_fast_xyz2=size(row_f2); %read size
                                    distance_count2=0;
                                    for dist_cmp_loop = 1: maxsize_fast_xyz2(1,1)
                                        xx_2 = xyz{1,fast_obj_start+2}(row_f2(dist_cmp_loop),1) - xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1);                              
                                        yy_2 = xyz{1,fast_obj_start+2}(row_f2(dist_cmp_loop),2) - xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2);                 
                                        distance_count2=distance_count2+1;
                                        distance_repeat_rest(distance_count2,1) = sqrt((xx_2^2)+(yy_2^2));
                                        distance_repeat_rest(distance_count2,2) = row_f2(dist_cmp_loop);
                                    end 
                                    if(dist_thresh_sec <= decide_slow)
                                        [row_secd,col_secd,check_secd]=find(distance_repeat_rest(:,1) < distance_threshold_fast);
                                    else
                                        [row_secd,col_secd,check_secd]=find(distance_repeat_rest(:,1) < dist_thresh_sec*threshold_diff_dislocation & distance_repeat_rest(:,1) < distance_threshold_fast);
                                    end
                                    if (check_secd == 1)
                                        row_secd_size=size(row_secd); %read size  
                                        angle_secd_fast=0;
                                        dist_thresh_third=0;
                                        priority_angle_order2=0;
                                        for find_secd_loop = 1: row_secd_size(1,1)
                                            angle_secd_fast_test = atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2))/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)))*180/pi;
                                            priority_angle_order2(find_secd_loop,1)=abs(abs(angle_fast)-abs(angle_secd_fast_test))/priority_degree_rate + distance_repeat_rest(row_secd(find_secd_loop));
                                            if( (xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)) == 0 & (xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)) == 0 )
                                                priority_angle_order2(find_secd_loop,1) = priority_zero; 
                                            end
                                            priority_angle_order2(find_secd_loop,2)=find_secd_loop;%give priority for smaller angle!!
                                            [temp_pr1,ord_pr1] = sort(priority_angle_order2(:,1));
                                            priority_angle_order2 = priority_angle_order2(ord_pr1,:);
                                        end
                                        for find_secd_loop2 = 1: row_secd_size(1,1)
                                            if(stop_search == 0)
                                                pass_numfil=0;
                                                diff_vector2=0;
                                                posit_nega2=nan;
                                                if( fast_obj_start > 1)
                                                    if(angle_fast <= 45 & angle_fast >= -45) %vetorc
                                                        vector=1; % x vetorc
                                                        diff_vector = xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1);
                                                        posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                    else 
                                                        vector=0; % y vetorc
                                                        diff_vector = xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2);
                                                        posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                    end
                                                end
                                                if(vector == 1)
                                                    diff_vector2= xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1);
                                                    posit_nega2=sign(diff_vector2);
                                                elseif (vector == 0)
                                                    diff_vector2= xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2);                            
                                                    posit_nega2=sign(diff_vector2);
                                                end
                                                if(posit_nega == posit_nega2 | posit_nega2 == 0 | posit_nega == 0)
                                                    increase_check=1;
                                                    if(fast_obj_start == 1)
                                                        if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start}(inside_loop,1))^2+(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start}(inside_loop,2))^2 ))
                                                            increase_check=0;
                                                        else
                                                            increase_check=1;
                                                        end
                                                    else
                                                        if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1))^2+(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))^2 ))
                                                            increase_check=0;
                                                        else
                                                            increase_check=1;
                                                        end
                                                    end
                                                    if(increase_check == 0)
                                                        angle_secd_fast_cmp=atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2))/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)))*180/pi; 
                                                        if( (xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)) == 0 & (xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)) == 0 )
                                                            angle_secd_fast_cmp=0;
                                                        end
                                                        if(fast_obj_start==1)
                                                            if(dist_thresh_sec > decide_slow_angle)
                                                                if(abs(angle_fast)>tangent_control | abs(angle_secd_fast_cmp) > tangent_control)                                                     
                                                                    if(abs(angle_fast)+control_vecter_angle > abs(atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start}(inside_loop,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start}(inside_loop,1)))*180/pi) & abs(angle_fast)-control_vecter_angle < abs(atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start}(inside_loop,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start}(inside_loop,1)))*180/pi)) % unit is degree                       
                                                                        pass_numfil=1;
                                                                    end
                                                                else
                                                                    if(angle_fast+control_vecter_angle > atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start}(inside_loop,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start}(inside_loop,1)))*180/pi & angle_fast-control_vecter_angle < atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start}(inside_loop,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start}(inside_loop,1)))*180/pi) % unit is degree  
                                                                        pass_numfil=1;
                                                                    end
                                                                end
                                                            else
                                                                if(abs(angle_fast)>tangent_control | abs(angle_secd_fast_cmp) > tangent_control)                                                     
                                                                    if(abs(angle_fast)+control_vecter_angle*multilple_short_angle > abs(atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start}(inside_loop,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start}(inside_loop,1)))*180/pi) & abs(angle_fast)-control_vecter_angle*multilple_short_angle < abs(atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start}(inside_loop,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start}(inside_loop,1)))*180/pi)) % unit is degree                       
                                                                        pass_numfil=1;
                                                                    end
                                                                else
                                                                    if(angle_fast+control_vecter_angle*multilple_short_angle > atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start}(inside_loop,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start}(inside_loop,1)))*180/pi & angle_fast-control_vecter_angle*multilple_short_angle < atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start}(inside_loop,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start}(inside_loop,1)))*180/pi) % unit is degree  
                                                                        pass_numfil=1;
                                                                    end
                                                                end
                                                            end
                                                        else
                                                            if(dist_thresh_sec > decide_slow_angle)
                                                                if(abs(angle_fast)>tangent_control | abs(angle_secd_fast_cmp) > tangent_control)                                                     
                                                                    if(abs(angle_fast)+control_vecter_angle > abs(atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi) & abs(angle_fast)-control_vecter_angle < abs(atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi)) % unit is degree                       
                                                                        pass_numfil=1;
                                                                    end
                                                                else
                                                                    if(angle_fast+control_vecter_angle > atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi & angle_fast-control_vecter_angle < atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi) % unit is degree  
                                                                        pass_numfil=1;
                                                                    end
                                                                end
                                                            else
                                                                if(abs(angle_fast)>tangent_control | abs(angle_secd_fast_cmp) > tangent_control)                                                     
                                                                    if(abs(angle_fast)+control_vecter_angle*multilple_short_angle > abs(atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi) & abs(angle_fast)-control_vecter_angle*multilple_short_angle < abs(atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi)) % unit is degree                       
                                                                        pass_numfil=1;
                                                                    end
                                                                else
                                                                    if(angle_fast+control_vecter_angle*multilple_short_angle > atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi & angle_fast-control_vecter_angle*multilple_short_angle < atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2) )/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi) % unit is degree  
                                                                        pass_numfil=1;
                                                                    end
                                                                end
                                                            end
                                                        end
                                                        if(pass_numfil == 1) 
                                                            angle_secd_fast = atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2))/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)))*180/pi;   
                                                            if( (xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2)) == 0 & (xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)) == 0 )
                                                                angle_secd_fast=0;
                                                            end
                                                            dist_thresh_third = distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)));                           
                                                            % finish third
                                                            distance_thrd_mat=0;
                                                            [row_f3,col_f3,check_f3]=find( xyz{1,fast_obj_start+3}(:,1) < xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1) + distance_threshold_fast & xyz{1,fast_obj_start+3}(:,1) > xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1) - distance_threshold_fast);       

                                                            if(check_f3 == 1)
                                                                maxsize_fast_xyz3=size(row_f3);
                                                                distance_count3=0;
                                                                for dist_thrd_loop = 1: maxsize_fast_xyz3(1,1)                                   
                                                                    xx_3 = xyz{1,fast_obj_start+3}(row_f3(dist_thrd_loop),1) - xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1);
                                                                    yy_3 = xyz{1,fast_obj_start+3}(row_f3(dist_thrd_loop),2) - xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2);
                                                                    distance_count3=distance_count3+1;
                                                                    distance_thrd_mat(distance_count3,1) = sqrt((xx_3^2)+(yy_3^2));
                                                                    distance_thrd_mat(distance_count3,2) = row_f3(dist_thrd_loop);
                                                                end
                                                                if(dist_thresh_third <= decide_slow)
                                                                    [row_third,col_third,check_third]=find(distance_thrd_mat(:,1) < distance_threshold_fast );
                                                                else
                                                                    [row_third,col_third,check_third]=find(distance_thrd_mat(:,1) < dist_thresh_third*threshold_diff_dislocation & distance_thrd_mat(:,1) <= distance_threshold_fast );
                                                                end
                                                                if (check_third == 1)
                                                                    row_third_size=size(row_third);
                                                                    angle_rest_th=0;
                                                                    dist_rest_th=0;                         
                                                                    priority_angle_order3=0;
                                                                    for find_rest_loop = 1: row_third_size(1,1)
                                                                        angle_rest_th_test=atan((xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2))/(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)))*180/pi;                                                           
                                                                        priority_angle_order3(find_rest_loop,1)=abs(abs(angle_secd_fast)-abs(angle_rest_th_test))/priority_degree_rate + distance_thrd_mat(row_third(find_rest_loop));
                                                                        if( (xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)) == 0 & (xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)) ==0 )
                                                                            priority_angle_order3(find_rest_loop,1) = priority_zero;
                                                                        end
                                                                        priority_angle_order3(find_rest_loop,2)=find_rest_loop;%give priority for smaller angle!!
                                                                        [temp_pr1,ord_pr1] = sort(priority_angle_order3(:,1));
                                                                        priority_angle_order3 = priority_angle_order3(ord_pr1,:);
                                                                    end
                                                                    for find_rest_loop2 = 1: row_third_size(1,1) 
                                                                        if(stop_search == 0)
                                                                            pass_numfil=0; 
                                                                            diff_vector2=0;
                                                                            posit_nega2=nan;
                                                                            if(angle_secd_fast <= 45 & angle_secd_fast >= -45) %vetorc
                                                                                vector=1; % x vetorc
                                                                                diff_vector = xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1);
                                                                                posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                            else 
                                                                                vector=0; % y vetorc
                                                                                diff_vector = xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2);
                                                                                posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                            end
                                                                            if(vector == 1)       
                                                                                diff_vector2= xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1) - xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1);
                                                                                posit_nega2=sign(diff_vector2);                                           
                                                                            elseif(vector == 0)                                            
                                                                                diff_vector2= xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) - xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2);                                            
                                                                                posit_nega2=sign(diff_vector2);
                                                                            end
                                                                            if(posit_nega == posit_nega2 | posit_nega2 == 0 | posit_nega == 0)
                                                                                increase_check=1;
                                                                                if((distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)))+dist_thresh_third)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1))^2+(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2))^2 ))
                                                                                    increase_check=0;                                                     
                                                                                else            
                                                                                    increase_check=1;
                                                                                end  
                                                                                if(increase_check == 0)  
                                                                                    angle_rest_th_cmp=atan((xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2))/(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)))*180/pi;
                                                                                    if( (xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)) == 0 & (xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1))==0)
                                                                                        angle_rest_th_cmp=0;
                                                                                    end
                                                                                    if( dist_thresh_third > decide_slow_angle)
                                                                                        if(abs(angle_secd_fast)>tangent_control | abs(angle_rest_th_cmp)>tangent_control)    
                                                                                            if(abs(angle_secd_fast)+control_vecter_angle > abs(atan((xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2) )/(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)))*180/pi) & abs(angle_secd_fast)-control_vecter_angle < abs(atan((xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2) )/(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)))*180/pi))  
                                                                                                pass_numfil=1;
                                                                                            end
                                                                                        else
                                                                                            if(angle_secd_fast+control_vecter_angle > atan((xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2) )/(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)))*180/pi & angle_secd_fast-control_vecter_angle < atan((xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2) )/(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)))*180/pi)  
                                                                                                pass_numfil=1;
                                                                                            end
                                                                                        end
                                                                                    else
                                                                                        if(abs(angle_secd_fast)>tangent_control | abs(angle_rest_th_cmp)>tangent_control)    
                                                                                            if(abs(angle_secd_fast)+control_vecter_angle*multilple_short_angle > abs(atan((xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2) )/(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)))*180/pi) & abs(angle_secd_fast)-control_vecter_angle*multilple_short_angle < abs(atan((xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2) )/(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)))*180/pi))  
                                                                                                pass_numfil=1;
                                                                                            end
                                                                                        else
                                                                                            if(angle_secd_fast+control_vecter_angle*multilple_short_angle > atan((xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2) )/(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)))*180/pi & angle_secd_fast-control_vecter_angle*multilple_short_angle < atan((xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2) )/(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)-xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1)))*180/pi)  
                                                                                                pass_numfil=1;
                                                                                            end
                                                                                        end 
                                                                                    end
                                                                                    if( pass_numfil == 1) 
                                                                                        angle_rest_th=atan((xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2))/(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)))*180/pi;                 
                                                                                        if ( (xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2)) == 0 & (xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)) == 0 )
                                                                                            angle_rest_th = 0;
                                                                                        end
                                                                                        dist_rest_th=distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)));
                                                                                        [row_f4,col_f4,check_f4]=find( xyz{1,fast_obj_start+4}(:,1) < xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1) + distance_threshold_fast & xyz{1,fast_obj_start+4}(:,1) > xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1) - distance_threshold_fast);       

                                                                                        if(check_f4 == 1)
                                                                                            distance_four_mat=0;
                                                                                            maxsize_fast_xyz4=size(row_f4);
                                                                                            distance_count4=0;
                                                                                            for dist_four_loop = 1: maxsize_fast_xyz4(1,1)
                                                                                                xx_4 = xyz{1,fast_obj_start+4}(row_f4(dist_four_loop),1) - xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1);    
                                                                                                yy_4 = xyz{1,fast_obj_start+4}(row_f4(dist_four_loop),2) - xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2);                                
                                                                                                distance_count4=distance_count4+1;
                                                                                                distance_four_mat(distance_count4,1) = sqrt((xx_4^2)+(yy_4^2));
                                                                                                distance_four_mat(distance_count4,2) = row_f4(dist_four_loop);
                                                                                            end
                                                                                            if(dist_rest_th <= decide_slow)     
                                                                                                [row_four,col_four,check_four]=find(distance_four_mat(:,1) < distance_threshold_fast);                              
                                                                                            else
                                                                                                [row_four,col_four,check_four]=find(distance_four_mat(:,1) < dist_rest_th*threshold_diff_dislocation & distance_four_mat(:,1) <= distance_threshold_fast);   
                                                                                            end
                                                                                            if (check_four == 1)
                                                                                                row_four_size=size(row_four);
                                                                                                angle_final_th=0;
                                                                                                dist_final_th=0; 
                                                                                                priority_angle_order4=0;
                                                                                                for find_four_loop = 1: row_four_size(1,1)
                                                                                                    angle_final_th_test=atan((xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2))/(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)))*180/pi;
                                                                                                    priority_angle_order4(find_four_loop,1)=abs(abs(angle_rest_th)-abs(angle_final_th_test))/priority_degree_rate + distance_four_mat(row_four(find_four_loop));
                                                                                                    if( (xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)) == 0 & (xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) == 0)
                                                                                                        priority_angle_order4(find_four_loop,1)=priority_zero;
                                                                                                    end
                                                                                                    priority_angle_order4(find_four_loop,2)=find_four_loop;%give priority for smaller angle!!
                                                                                                    [temp_pr1,ord_pr1] = sort(priority_angle_order4(:,1));
                                                                                                    priority_angle_order4 = priority_angle_order4(ord_pr1,:);
                                                                                                end
                                                                                                for find_four_loop2 = 1: row_four_size(1,1)                                                          
                                                                                                    if(stop_search == 0)
                                                                                                        pass_numfil=0;
                                                                                                        diff_vector2=0;
                                                                                                        posit_nega2=nan;
                                                                                                        if(angle_rest_th <= 45 & angle_rest_th >= -45) %vetorc
                                                                                                            vector=1; % x vetorc
                                                                                                            diff_vector = xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1);
                                                                                                            posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                        else 
                                                                                                            vector=0; % y vetorc
                                                                                                            diff_vector = xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2);
                                                                                                            posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                        end
                                                                                                        if(vector == 1)       
                                                                                                            diff_vector2= xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1) - xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1);      
                                                                                                            posit_nega2=sign(diff_vector2);        
                                                                                                        elseif(vector == 0)                                            
                                                                                                            diff_vector2= xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2) - xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2);                                            
                                                                                                            posit_nega2=sign(diff_vector2);
                                                                                                        end
                                                                                                        if(posit_nega == posit_nega2 | posit_nega2 == 0 | posit_nega == 0)
                                                                                                            increase_check=1;
                                                                                                            if((distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)))+dist_rest_th)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1))^2+(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2))^2 ))       
                                                                                                                increase_check=0;    
                                                                                                            else
                                                                                                                increase_check=1;      
                                                                                                            end
                                                                                                            if(increase_check == 0)                                                                               
                                                                                                                angle_final_th_cmp=atan((xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2))/(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)))*180/pi;
                                                                                                                if( (xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)) == 0 & (xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) ==0 )
                                                                                                                    angle_final_th_cmp=0;
                                                                                                                end
                                                                                                                if( dist_rest_th > decide_slow_angle)
                                                                                                                    if(abs(angle_rest_th)>tangent_control | abs(angle_final_th_cmp) >tangent_control)
                                                                                                                        if(abs(angle_rest_th)+control_vecter_angle > abs(atan((xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)))*180/pi) & abs(angle_rest_th)-control_vecter_angle < abs(atan((xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)))*180/pi))
                                                                                                                            pass_numfil=1;
                                                                                                                        end
                                                                                                                    else
                                                                                                                        if(angle_rest_th+control_vecter_angle > atan((xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)))*180/pi & angle_rest_th-control_vecter_angle < atan((xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)))*180/pi)
                                                                                                                            pass_numfil=1;
                                                                                                                        end
                                                                                                                    end
                                                                                                                else  
                                                                                                                     if(abs(angle_rest_th)>tangent_control | abs(angle_final_th_cmp) >tangent_control)
                                                                                                                        if(abs(angle_rest_th)+control_vecter_angle*multilple_short_angle > abs(atan((xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)))*180/pi) & abs(angle_rest_th)-control_vecter_angle*multilple_short_angle  < abs(atan((xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)))*180/pi))
                                                                                                                            pass_numfil=1;
                                                                                                                        end
                                                                                                                    else
                                                                                                                        if(angle_rest_th+control_vecter_angle*multilple_short_angle  > atan((xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)))*180/pi & angle_rest_th-control_vecter_angle*multilple_short_angle  < atan((xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1)))*180/pi)
                                                                                                                            pass_numfil=1;
                                                                                                                        end
                                                                                                                    end
                                                                                                                end
                                                                                                                if( pass_numfil == 1)
                                                                                                                    angle_final_th=atan((xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2))/(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)))*180/pi;                                    
                                                                                                                    if( (xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2)) == 0 & (xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1))==0)
                                                                                                                        angle_final_th=0;
                                                                                                                    end
                                                                                                                    dist_final_th=distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)));
                                                                                                                    if(fast_obj_start == 1)
                                                                                                                        object_number_count = object_number_count+1;
                                                                                                                        fst_mv_oj_coor=0;

                                                                                                                        fst_mv_oj_coor(1,1) = xyz{1,fast_obj_start}(inside_loop,1);                              
                                                                                                                        fst_mv_oj_coor(1,2) = xyz{1,fast_obj_start}(inside_loop,2);          
                                                                                                                        fst_mv_oj_coor(1,3) = object_number_count;

                                                                                                                        fst_mv_oj_coor(2,1) = xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1);
                                                                                                                        fst_mv_oj_coor(2,2) = xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2);       
                                                                                                                        fst_mv_oj_coor(2,3) = object_number_count;

                                                                                                                        fst_mv_oj_coor(3,1) = xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1);
                                                                                                                        fst_mv_oj_coor(3,2) = xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2);          
                                                                                                                        fst_mv_oj_coor(3,3) = object_number_count;

                                                                                                                        fst_mv_oj_coor(4,1) = xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1);
                                                                                                                        fst_mv_oj_coor(4,2) = xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2);          
                                                                                                                        fst_mv_oj_coor(4,3) = object_number_count;

                                                                                                                        fst_mv_oj_coor(5,1) = xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1);
                                                                                                                        fst_mv_oj_coor(5,2) = xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2);
                                                                                                                        fst_mv_oj_coor(5,3) = object_number_count;

                                                                                                                        if( abs(angle_final_th) > tangent_control)
                                                                                                                            angle_information(object_number_count,1)=(abs(abs(angle_fast)-abs(angle_secd_fast))+abs(abs(angle_secd_fast)-abs(angle_rest_th))+abs(abs(angle_rest_th)-abs(angle_final_th)))/priority_degree_rate + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) );
                                                                                                                        else
                                                                                                                            angle_information(object_number_count,1)=(abs(angle_fast-angle_secd_fast)+abs(angle_secd_fast-angle_rest_th)+abs(angle_rest_th-angle_final_th))/priority_degree_rate + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) );
                                                                                                                        end
                                                                                                                        distance_six_mat=0;
                                                                                                                        [row_f5,col_f5,check_f5]=find( xyz{1,fast_obj_start+4}(:,1) < xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1) + distance_threshold_fast & xyz{1,fast_obj_start+4}(:,1) > xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1) - distance_threshold_fast);       

                                                                                                                        if(check_f5 == 1)
                                                                                                                            stop_search2=0;
                                                                                                                            maxsize_fast_xyz6=size(row_f5);
                                                                                                                            distance_count6=0;
                                                                                                                            for dist_six_loop = 1: maxsize_fast_xyz6(1,1)                                   
                                                                                                                                xx_6 = xyz{1,fast_obj_start+4}(row_f5(dist_six_loop),1) - xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1);
                                                                                                                                yy_6 = xyz{1,fast_obj_start+4}(row_f5(dist_six_loop),2) - xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2);                                                                                                                                       
                                                                                                                                distance_count6=distance_count6+1;
                                                                                                                                distance_six_mat(distance_count6,1) = sqrt((xx_6^2)+(yy_6^2));
                                                                                                                                distance_six_mat(distance_count6,2) = row_f5(dist_six_loop);
                                                                                                                            end
                                                                                                                            if(dist_final_th <= decide_slow)
                                                                                                                                [row_six,col_six,check_six]=find(distance_six_mat(:,1) <= distance_threshold_fast & distance_six_mat(:,1) ~= 0);
                                                                                                                            else
                                                                                                                                [row_six,col_six,check_six]=find(distance_six_mat(:,1) <= dist_final_th*threshold_diff_dislocation & distance_six_mat(:,1) <= distance_threshold_fast & distance_six_mat(:,1) ~= 0);
                                                                                                                            end
                                                                                                                            if (check_six == 1)
                                                                                                                                row_six_size=size(row_six);
                                                                                                                                angle_six_th=0;
                                                                                                                                dist_six_th=0;                         
                                                                                                                                priority_angle_order6=0;
                                                                                                                                for find_six_loop = 1: row_six_size(1,1)
                                                                                                                                    angle_six_th_test=atan((xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),2)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2))/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),1)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)))*180/pi;                                                           
                                                                                                                                    priority_angle_order6(find_six_loop,1)=abs(abs(angle_final_th)-abs(angle_six_th_test))/priority_degree_rate + distance_six_mat(row_six(find_six_loop));
                                                                                                                                    if( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),2)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)) == 0 & (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),1)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)) )
                                                                                                                                        priority_angle_order6(find_six_loop,1)= priority_zero; 
                                                                                                                                    end
                                                                                                                                    priority_angle_order6(find_six_loop,2)=find_six_loop;%give priority for smaller angle!!
                                                                                                                                    [temp_pr1,ord_pr1] = sort(priority_angle_order6(:,1));
                                                                                                                                    priority_angle_order6 = priority_angle_order6(ord_pr1,:);
                                                                                                                                end
                                                                                                                                for find_six_loop2 = 1: row_six_size(1,1)
                                                                                                                                    pass_numfil=0;  
                                                                                                                                    if(stop_search2 == 0)
                                                                                                                                        diff_vector2=0;
                                                                                                                                        posit_nega2=nan;
                                                                                                                                        if(angle_final_th <= 45 & angle_final_th >= -45) %vetorc
                                                                                                                                            vector=1; % x vetorc
                                                                                                                                            diff_vector = xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1);
                                                                                                                                            posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                                        else 
                                                                                                                                            vector=0; % y vetorc
                                                                                                                                            diff_vector = xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2);
                                                                                                                                            posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                                        end
                                                                                                                                        if(vector == 1)       
                                                                                                                                            diff_vector2= xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1) - xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1);
                                                                                                                                            posit_nega2=sign(diff_vector2);                                           
                                                                                                                                        elseif(vector == 0)                                            
                                                                                                                                            diff_vector2= xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2) - xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2);
                                                                                                                                            posit_nega2=sign(diff_vector2);
                                                                                                                                        end
                                                                                                                                        if(posit_nega == posit_nega2 | posit_nega2 == 0 | posit_nega == 0)
                                                                                                                                            increase_check=1; 
                                                                                                                                            if((distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)))+dist_final_th)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1))^2+(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2))^2 ))
                                                                                                                                                increase_check=0;                                                     
                                                                                                                                            else            
                                                                                                                                                increase_check=1;
                                                                                                                                            end  
                                                                                                                                            if(increase_check == 0)
                                                                                                                                                angle_six_th_cmp= atan((xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2))/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)))*180/pi;
                                                                                                                                                if( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)) == 0 & (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)) == 0 )
                                                                                                                                                    angle_six_th_cmp=0;
                                                                                                                                                end
                                                                                                                                                if( dist_final_th > decide_slow_angle)
                                                                                                                                                    if(abs(angle_final_th)>tangent_control | abs(angle_six_th_cmp)>tangent_control)    
                                                                                                                                                        if(abs(angle_final_th)+control_vecter_angle > abs(atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi) & abs(angle_final_th)-control_vecter_angle < abs(atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi))  
                                                                                                                                                            pass_numfil=1;
                                                                                                                                                        end
                                                                                                                                                    else
                                                                                                                                                        if(angle_final_th+control_vecter_angle > atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi & angle_final_th-control_vecter_angle < atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi)  
                                                                                                                                                            pass_numfil=1;
                                                                                                                                                        end
                                                                                                                                                    end
                                                                                                                                                else
                                                                                                                                                    if(abs(angle_final_th)>tangent_control | abs(angle_six_th_cmp)>tangent_control)
                                                                                                                                                        if(abs(angle_final_th)+control_vecter_angle*multilple_short_angle > abs(atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi) & abs(angle_final_th)-control_vecter_angle*multilple_short_angle < abs(atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi))  
                                                                                                                                                            pass_numfil=1;
                                                                                                                                                        end
                                                                                                                                                    else
                                                                                                                                                        if(angle_final_th+control_vecter_angle*multilple_short_angle > atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi & angle_final_th-control_vecter_angle*multilple_short_angle < atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi)  
                                                                                                                                                            pass_numfil=1;
                                                                                                                                                        end
                                                                                                                                                    end
                                                                                                                                                end
                                                                                                                                                if( pass_numfil == 1)
                                                                                                                                                    fst_mv_oj_coor(6,1) = xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1);
                                                                                                                                                    fst_mv_oj_coor(6,2) = xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2);
                                                                                                                                                    fst_mv_oj_coor(6,3) = object_number_count;
                                                                                                                                                    stop_search2=1;
                                                                                                                                                end
                                                                                                                                            end
                                                                                                                                        end
                                                                                                                                    end
                                                                                                                                end
                                                                                                                            end                                                                                                                                
                                                                                                                        end
                                                                                                                        %-----prevent duplication----
                                                                                                                        if( fast_total_coor{1,1} == 0)
                                                                                                                            fast_total_coor(1,object_number_count)={fst_mv_oj_coor};
                                                                                                                            stop_search=1; %prevent duplicate trajectories
                                                                                                                            cont_rest_th(object_number_count,1)=dist_final_th;
                                                                                                                            cont_angle_th(object_number_count,1)=angle_final_th;
                                                                                                                            if(abs(angle_final_th) > tangent_control)
                                                                                                                                cont_angle_cmp(object_number_count,1)=abs(angle_rest_th)-abs(angle_final_th);
                                                                                                                            else
                                                                                                                                cont_angle_cmp(object_number_count,1)=angle_rest_th-angle_final_th;
                                                                                                                            end
                                                                                                                        else
                                                                                                                            %------prevent duplication--------
                                                                                                                            maxsize_fstc=size(fast_total_coor);     
                                                                                                                            maxsize_cmp_fst=size(fst_mv_oj_coor);
                                                                                                                            duplication_number=0;
                                                                                                                            duplecount=0;
                                                                                                                            duplehappen=0;
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
                                                                                                                                    duplecount=duplecount+1;
                                                                                                                                    duplication_number(duplecount,1)=prevent_duplication_loop;
                                                                                                                                end
                                                                                                                            end
                                                                                                                            if(duplecount ~= 0)
                                                                                                                                for check_duple_loop=1:duplecount
                                                                                                                                    if(angle_information(duplication_number(check_duple_loop,1),1) < angle_information(object_number_count,1))
                                                                                                                                        duplehappen=1;
                                                                                                                                    end
                                                                                                                                end
                                                                                                                                if( duplehappen == 1)
                                                                                                                                    object_number_count= object_number_count-1; %duplication happen
                                                                                                                                else 
                                                                                                                                    for check_duple_loop=1:duplecount
                                                                                                                                        fast_total_coor(1,duplication_number(check_duple_loop,1)) = {0};
                                                                                                                                        cont_rest_th(object_number_count,1)=dist_final_th;
                                                                                                                                        cont_angle_th(object_number_count,1)=angle_final_th;
                                                                                                                                        if(abs(angle_final_th) > tangent_control)
                                                                                                                                            cont_angle_cmp(object_number_count,1)=abs(angle_rest_th)-abs(angle_final_th);
                                                                                                                                        else
                                                                                                                                            cont_angle_cmp(object_number_count,1)=angle_rest_th-angle_final_th;
                                                                                                                                        end      
                                                                                                                                        cont_rest_th(duplication_number(check_duple_loop,1),1)=cont_rest_th(object_number_count,1);
                                                                                                                                        cont_angle_th(duplication_number(check_duple_loop,1),1)=cont_angle_th(object_number_count,1);  
                                                                                                                                        cont_angle_cmp(duplication_number(check_duple_loop,1),1)=cont_angle_cmp(object_number_count,1);
                                                                                                                                        object_number_count = object_number_count-1;
                                                                                                                                        stop_search=1;
                                                                                                                                        maxsize_fst_check=size(fst_mv_oj_coor);
                                                                                                                                        for thirdnumberchange=1:maxsize_fst_check(1,1)
                                                                                                                                            fst_mv_oj_coor(thirdnumberchange,3) = duplication_number(check_duple_loop,1);
                                                                                                                                        end
                                                                                                                                        fast_total_coor(1,duplication_number(check_duple_loop,1)) = {fst_mv_oj_coor}; %smaller angle is right trajectory
                                                                                                                                    end
                                                                                                                                end
                                                                                                                            else
                                                                                                                                fast_total_coor(1,object_number_count)={fst_mv_oj_coor};
                                                                                                                                stop_search=1; %prevent duplicate trajectories
                                                                                                                                cont_rest_th(object_number_count,1)=dist_final_th;
                                                                                                                                cont_angle_th(object_number_count,1)=angle_final_th;
                                                                                                                                if(abs(angle_secd_fast) > tangent_control)
                                                                                                                                    cont_angle_cmp(object_number_count,1)=abs(angle_fast)-abs(angle_secd_fast);
                                                                                                                                else
                                                                                                                                    cont_angle_cmp(object_number_count,1)=angle_fast-angle_secd_fast;
                                                                                                                                end
                                                                                                                            end
                                                                                                                        end
                                                                                                                    else
                                                                                                                        fst_mv_oj_coor=0;

                                                                                                                        fst_mv_oj_coor(1,1) = xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1);
                                                                                                                        fst_mv_oj_coor(1,2) = xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2);       
                                                                                                                        fst_mv_oj_coor(1,3) = inside_loop;

                                                                                                                        fst_mv_oj_coor(2,1) = xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1);
                                                                                                                        fst_mv_oj_coor(2,2) = xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2);          
                                                                                                                        fst_mv_oj_coor(2,3) = inside_loop;

                                                                                                                        fst_mv_oj_coor(3,1) = xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1);
                                                                                                                        fst_mv_oj_coor(3,2) = xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2); 
                                                                                                                        fst_mv_oj_coor(3,3) = inside_loop;

                                                                                                                        fst_mv_oj_coor(4,1) = xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1);
                                                                                                                        fst_mv_oj_coor(4,2) = xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2);
                                                                                                                        fst_mv_oj_coor(4,3) = inside_loop; 

                                                                                                                        if( abs(angle_final_th) > tangent_control)
                                                                                                                            angle_information(inside_loop,1)=(abs(abs(angle_fast)-abs(angle_secd_fast))+abs(abs(angle_secd_fast)-abs(angle_rest_th))+abs(abs(angle_rest_th)-abs(angle_final_th)))/priority_degree_rate + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) );
                                                                                                                        else
                                                                                                                            angle_information(inside_loop,1)=(abs(angle_fast-angle_secd_fast)+abs(angle_secd_fast-angle_rest_th)+abs(angle_rest_th-angle_final_th))/priority_degree_rate + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) );
                                                                                                                        end

                                                                                                                        [row_f5,col_f5,check_f5]=find( xyz{1,fast_obj_start+4}(:,1) < xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1) + distance_threshold_fast & xyz{1,fast_obj_start+4}(:,1) > xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1) - distance_threshold_fast);       
                                                                                                                        if(check_f5 == 1)
                                                                                                                            distance_six_mat=0;
                                                                                                                            stop_search2=0;
                                                                                                                            if( dist_final_th > decide_slow)                                                                                                            
                                                                                                                                maxsize_fast_xyz6=size(row_f5);
                                                                                                                                distance_count6=0;
                                                                                                                                for dist_six_loop = 1: maxsize_fast_xyz6(1,1)                                   
                                                                                                                                    xx_6 = xyz{1,fast_obj_start+4}(row_f5(dist_six_loop),1) - xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1);
                                                                                                                                    yy_6 = xyz{1,fast_obj_start+4}(row_f5(dist_six_loop),2) - xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2);                                                                                                                                       
                                                                                                                                    distance_count6=distance_count6+1;
                                                                                                                                    distance_six_mat(distance_count6,1) = sqrt((xx_6^2)+(yy_6^2));
                                                                                                                                    distance_six_mat(distance_count6,2) = row_f5(dist_six_loop);
                                                                                                                                end
                                                                                                                                if(dist_final_th <= decide_slow)
                                                                                                                                    [row_six,col_six,check_six]=find(distance_six_mat(:,1) <= distance_threshold_fast & distance_six_mat(:,1) ~= 0);
                                                                                                                                else
                                                                                                                                    [row_six,col_six,check_six]=find(distance_six_mat(:,1) <= dist_final_th*threshold_diff_dislocation & distance_six_mat(:,1) <= distance_threshold_fast & distance_six_mat(:,1) ~= 0);
                                                                                                                                end

                                                                                                                                if (check_six == 1)
                                                                                                                                    row_six_size=size(row_six);
                                                                                                                                    angle_six_th=0;
                                                                                                                                    dist_six_th=0;                         
                                                                                                                                    priority_angle_order6=0;
                                                                                                                                    for find_six_loop = 1: row_six_size(1,1)
                                                                                                                                        angle_six_th_test=atan((xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),2)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2))/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),1)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)))*180/pi;                                                           
                                                                                                                                        priority_angle_order6(find_six_loop,1)=abs(abs(angle_final_th)-abs(angle_six_th_test))/priority_degree_rate + distance_six_mat(row_six(find_six_loop));
                                                                                                                                        if( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),2)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)) == 0 & (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),1) - xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)) )
                                                                                                                                            priority_angle_order6(find_six_loop,1)=0;
                                                                                                                                        end
                                                                                                                                        priority_angle_order6(find_six_loop,2)=find_six_loop;%give priority for smaller angle!!
                                                                                                                                        [temp_pr1,ord_pr1] = sort(priority_angle_order6(:,1));
                                                                                                                                        priority_angle_order6 = priority_angle_order6(ord_pr1,:);
                                                                                                                                    end
                                                                                                                                    for find_six_loop2 = 1: row_six_size(1,1)
                                                                                                                                        pass_numfil=0;  
                                                                                                                                        if(stop_search2 == 0)
                                                                                                                                            diff_vector2=0;
                                                                                                                                            posit_nega2=nan;
                                                                                                                                            if(angle_final_th <= 45 & angle_final_th >= -45) %vetorc
                                                                                                                                                vector=1; % x vetorc
                                                                                                                                                diff_vector = xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1);
                                                                                                                                                posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                                            else 
                                                                                                                                                vector=0; % y vetorc
                                                                                                                                                diff_vector = xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2);
                                                                                                                                                posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                                            end
                                                                                                                                            if(vector == 1)       
                                                                                                                                                diff_vector2= xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1) - xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1);
                                                                                                                                                posit_nega2=sign(diff_vector2);                                           
                                                                                                                                            elseif(vector == 0)                                            
                                                                                                                                                diff_vector2= xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2) - xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2);
                                                                                                                                                posit_nega2=sign(diff_vector2);
                                                                                                                                            end
                                                                                                                                            if(posit_nega == posit_nega2 | posit_nega2 == 0 | posit_nega == 0)
                                                                                                                                                increase_check=1; 
                                                                                                                                                if((distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)))+dist_final_th)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1))^2+(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2))^2 ))
                                                                                                                                                    increase_check=0;                                                     
                                                                                                                                                else            
                                                                                                                                                    increase_check=1;
                                                                                                                                                end
                                                                                                                                                if(increase_check == 0)
                                                                                                                                                    angle_six_th_cmp= atan((xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2))/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)))*180/pi;
                                                                                                                                                    if( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2)) == 0 & (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)) == 0 )
                                                                                                                                                        angle_six_th_cmp=0;
                                                                                                                                                    end
                                                                                                                                                    if( dist_final_th > decide_slow_angle)
                                                                                                                                                        if(abs(angle_final_th)>tangent_control | abs(angle_six_th_cmp)>tangent_control)    
                                                                                                                                                            if(abs(angle_final_th)+control_vecter_angle > abs(atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi) & abs(angle_final_th)-control_vecter_angle < abs(atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi))  
                                                                                                                                                                pass_numfil=1;
                                                                                                                                                            end
                                                                                                                                                        else
                                                                                                                                                            if(angle_final_th+control_vecter_angle > atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi & angle_final_th-control_vecter_angle < atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi)  
                                                                                                                                                                pass_numfil=1;
                                                                                                                                                            end
                                                                                                                                                        end
                                                                                                                                                    else
                                                                                                                                                        if(abs(angle_final_th)>tangent_control | abs(angle_six_th_cmp)>tangent_control)
                                                                                                                                                            if(abs(angle_final_th)+control_vecter_angle*multilple_short_angle > abs(atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi) & abs(angle_final_th)-control_vecter_angle*multilple_short_angle < abs(atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi))  
                                                                                                                                                                pass_numfil=1;
                                                                                                                                                            end
                                                                                                                                                        else
                                                                                                                                                            if(angle_final_th+control_vecter_angle*multilple_short_angle > atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi & angle_final_th-control_vecter_angle*multilple_short_angle < atan( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2) )/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1)) )*180/pi)  
                                                                                                                                                                pass_numfil=1;
                                                                                                                                                            end
                                                                                                                                                        end
                                                                                                                                                    end
                                                                                                                                                    if( pass_numfil == 1)
                                                                                                                                                        fst_mv_oj_coor(5,1) = xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1);
                                                                                                                                                        fst_mv_oj_coor(5,2) = xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2);
                                                                                                                                                        stop_search2=1;
                                                                                                                                                    end
                                                                                                                                                end
                                                                                                                                            end
                                                                                                                                        end
                                                                                                                                    end
                                                                                                                                end
                                                                                                                            end
                                                                                                                        end
                                                                                                                        %-----prevent duplication----
                                                                                                                        if( prevent_rest_check{1,1} == 0)

                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+1,1) = fst_mv_oj_coor(1,1);
                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+1,2) = fst_mv_oj_coor(1,2);       
                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+1,3) = inside_loop;

                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+2,1) = fst_mv_oj_coor(2,1);
                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+2,2) = fst_mv_oj_coor(2,2);         
                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+2,3) = inside_loop;

                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+3,1) = fst_mv_oj_coor(3,1);
                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+3,2) = fst_mv_oj_coor(3,2);       
                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+3,3) = inside_loop;

                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+4,1) = fst_mv_oj_coor(4,1);
                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+4,2) = fst_mv_oj_coor(4,2);
                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+4,3) = inside_loop;  

                                                                                                                            prevent_rest_check{1,inside_loop}(1,1) = fast_total_coor{1,inside_loop}(fast_obj_start,1);
                                                                                                                            prevent_rest_check{1,inside_loop}(1,2) = fast_total_coor{1,inside_loop}(fast_obj_start,2); 
                                                                                                                            prevent_rest_check{1,inside_loop}(2,1) = fst_mv_oj_coor(1,1);
                                                                                                                            prevent_rest_check{1,inside_loop}(2,2) = fst_mv_oj_coor(1,2); 
                                                                                                                            prevent_rest_check{1,inside_loop}(3,1) = fst_mv_oj_coor(2,1);
                                                                                                                            prevent_rest_check{1,inside_loop}(3,2) = fst_mv_oj_coor(2,2); 
                                                                                                                            prevent_rest_check{1,inside_loop}(4,1) = fst_mv_oj_coor(3,1);
                                                                                                                            prevent_rest_check{1,inside_loop}(4,2) = fst_mv_oj_coor(3,2);
                                                                                                                            prevent_rest_check{1,inside_loop}(5,1) = fst_mv_oj_coor(4,1);
                                                                                                                            prevent_rest_check{1,inside_loop}(5,2) = fst_mv_oj_coor(4,2);
                                                                                                                            if(stop_search2==1)
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+5,1) = fst_mv_oj_coor(5,1);
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+5,2) = fst_mv_oj_coor(5,2);
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+5,3) = inside_loop;
                                                                                                                                prevent_rest_check{1,inside_loop}(6,1)= fst_mv_oj_coor(5,1);
                                                                                                                                prevent_rest_check{1,inside_loop}(6,2)= fst_mv_oj_coor(5,2);
                                                                                                                            end
                                                                                                                            stop_search=1; %prevent duplicate trajectories
                                                                                                                            cont_rest_th(inside_loop,1)=dist_final_th;
                                                                                                                            cont_angle_th(inside_loop,1)=angle_final_th;
                                                                                                                            if(abs(angle_final_th) > tangent_control)
                                                                                                                                cont_angle_cmp(inside_loop,1)=abs(angle_rest_th)-abs(angle_final_th);
                                                                                                                            else
                                                                                                                                cont_angle_cmp(inside_loop,1)=angle_rest_th-angle_final_th;
                                                                                                                            end
                                                                                                                        else
                                                                                                                            %------prevent duplication--------
                                                                                                                            maxsize_fstc=size(prevent_rest_check);     
                                                                                                                            maxsize_cmp_fst=size(fst_mv_oj_coor);
                                                                                                                            duplication_number=0;
                                                                                                                            duplecount=0;
                                                                                                                            duplehappen=0;
                                                                                                                            for prevent_duplication_loop=1:maxsize_fstc(1,2)	
                                                                                                                                duplication_happen=0;
                                                                                                                                for fst_check_loop =1:maxsize_cmp_fst(1,1)
                                                                                                                                    [row_dplc_x,col_dplc_x,check_dplc_x]=find(prevent_rest_check{1,prevent_duplication_loop}(:,1) == fst_mv_oj_coor(fst_check_loop,1));
                                                                                                                                    [row_dplc_y,col_dplc_y,check_dplc_y]=find(prevent_rest_check{1,prevent_duplication_loop}(:,2) == fst_mv_oj_coor(fst_check_loop,2));
                                                                                                                                    rr_dplc=size(row_dplc_x);
                                                                                                                                    for dplc_chk_loop = 1:rr_dplc
                                                                                                                                        [row_find_dplc,col_find_dplc,check_find_dplc]= find(row_dplc_y == row_dplc_x(dplc_chk_loop));
                                                                                                                                        if(check_find_dplc ~= 0)
                                                                                                                                            duplication_happen=duplication_happen+1;
                                                                                                                                        end
                                                                                                                                    end
                                                                                                                                end
                                                                                                                                if(duplication_happen >= 2) 
                                                                                                                                    duplecount=duplecount+1;
                                                                                                                                    duplication_number(duplecount,1)=prevent_duplication_loop;
                                                                                                                                end
                                                                                                                            end
                                                                                                                            if(duplecount ~= 0)
                                                                                                                                for check_duple_loop=1:duplecount
                                                                                                                                    if(angle_information(duplication_number(check_duple_loop,1),1) < angle_information(inside_loop,1))
                                                                                                                                        duplehappen=1;
                                                                                                                                    end
                                                                                                                                end
                                                                                                                                if( duplehappen == 1)
                                                                                                                                    %duplication happen
                                                                                                                                else 
                                                                                                                                    for check_duple_loop=1:duplecount
                                                                                                                                        cont_rest_th(inside_loop,1)=dist_final_th;
                                                                                                                                        cont_angle_th(inside_loop,1)=angle_final_th;
                                                                                                                                        if(abs(angle_final_th) > tangent_control)
                                                                                                                                            cont_angle_cmp(inside_loop,1)=abs(angle_rest_th)-abs(angle_final_th);
                                                                                                                                        else
                                                                                                                                            cont_angle_cmp(inside_loop,1)=angle_rest_th-angle_final_th;
                                                                                                                                        end      
                                                                                                                                        cont_rest_th(duplication_number(check_duple_loop,1),1)=cont_rest_th(inside_loop,1);
                                                                                                                                        cont_angle_th(duplication_number(check_duple_loop,1),1)=cont_angle_th(inside_loop,1);  
                                                                                                                                        cont_angle_cmp(duplication_number(check_duple_loop,1),1)=cont_angle_cmp(inside_loop,1);
                                                                                                                                        stop_search=1;
                                                                                                                                        maxsize_fst_check=size(fst_mv_oj_coor);
                                                                                                                                        for thirdnumberchange=1:maxsize_fst_check(1,1)
                                                                                                                                            fst_mv_oj_coor(thirdnumberchange,3) = duplication_number(check_duple_loop,1);
                                                                                                                                        end
                                                                                                                                        %smaller angle is right trajectory
                                                                                                                                        fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+1,1) = fst_mv_oj_coor(1,1);
                                                                                                                                        fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+1,2) = fst_mv_oj_coor(1,2);       
                                                                                                                                        fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+1,3) = inside_loop;

                                                                                                                                        fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+2,1) = fst_mv_oj_coor(2,1);
                                                                                                                                        fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+2,2) = fst_mv_oj_coor(2,2);         
                                                                                                                                        fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+2,3) = inside_loop;

                                                                                                                                        fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+3,1) = fst_mv_oj_coor(3,1);
                                                                                                                                        fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+3,2) = fst_mv_oj_coor(3,2);       
                                                                                                                                        fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+3,3) = inside_loop;

                                                                                                                                        fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+4,1) = fst_mv_oj_coor(4,1);
                                                                                                                                        fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+4,2) = fst_mv_oj_coor(4,2);
                                                                                                                                        fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+4,3) = inside_loop;  

                                                                                                                                        prevent_rest_check{1,inside_loop}(1,1) = nan;
                                                                                                                                        prevent_rest_check{1,inside_loop}(1,2) = nan;

                                                                                                                                        prevent_rest_check{1,duplication_number(check_duple_loop,1)}(1,1) = fast_total_coor{1,inside_loop}(fast_obj_start,1);
                                                                                                                                        prevent_rest_check{1,duplication_number(check_duple_loop,1)}(1,2) = fast_total_coor{1,inside_loop}(fast_obj_start,2);
                                                                                                                                        prevent_rest_check{1,duplication_number(check_duple_loop,1)}(2,1) = fst_mv_oj_coor(1,1);
                                                                                                                                        prevent_rest_check{1,duplication_number(check_duple_loop,1)}(2,2) = fst_mv_oj_coor(1,2);  
                                                                                                                                        prevent_rest_check{1,duplication_number(check_duple_loop,1)}(3,1) = fst_mv_oj_coor(2,1);
                                                                                                                                        prevent_rest_check{1,duplication_number(check_duple_loop,1)}(3,2) = fst_mv_oj_coor(2,2);  
                                                                                                                                        prevent_rest_check{1,duplication_number(check_duple_loop,1)}(4,1) = fst_mv_oj_coor(3,1);
                                                                                                                                        prevent_rest_check{1,duplication_number(check_duple_loop,1)}(4,2) = fst_mv_oj_coor(3,2);  
                                                                                                                                        prevent_rest_check{1,duplication_number(check_duple_loop,1)}(5,1) = fst_mv_oj_coor(4,1);
                                                                                                                                        prevent_rest_check{1,duplication_number(check_duple_loop,1)}(5,2) = fst_mv_oj_coor(4,2);

                                                                                                                                        if(stop_search2==1)
                                                                                                                                            fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+5,1) = fst_mv_oj_coor(5,1);
                                                                                                                                            fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+5,2) = fst_mv_oj_coor(5,2);
                                                                                                                                            fast_total_coor{1,duplication_number(check_duple_loop,1)}(fast_obj_start+5,3) = inside_loop;

                                                                                                                                            prevent_rest_check{1,duplication_number(check_duple_loop,1)}(6,1)= fst_mv_oj_coor(5,1);
                                                                                                                                            prevent_rest_check{1,duplication_number(check_duple_loop,1)}(6,2)= fst_mv_oj_coor(5,2);
                                                                                                                                        end

                                                                                                                                    end
                                                                                                                                end
                                                                                                                            else
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+1,1) = fst_mv_oj_coor(1,1);
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+1,2) = fst_mv_oj_coor(1,2);       
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+1,3) = inside_loop;

                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+2,1) = fst_mv_oj_coor(2,1);
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+2,2) = fst_mv_oj_coor(2,2);         
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+2,3) = inside_loop;

                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+3,1) = fst_mv_oj_coor(3,1);
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+3,2) = fst_mv_oj_coor(3,2);       
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+3,3) = inside_loop;

                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+4,1) = fst_mv_oj_coor(4,1);
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+4,2) = fst_mv_oj_coor(4,2);
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start+4,3) = inside_loop;  

                                                                                                                                prevent_rest_check{1,inside_loop}(1,1) = fast_total_coor{1,inside_loop}(fast_obj_start,1);
                                                                                                                                prevent_rest_check{1,inside_loop}(1,2) = fast_total_coor{1,inside_loop}(fast_obj_start,2);
                                                                                                                                prevent_rest_check{1,inside_loop}(2,1) = fst_mv_oj_coor(1,1);
                                                                                                                                prevent_rest_check{1,inside_loop}(2,2) = fst_mv_oj_coor(1,2);  
                                                                                                                                prevent_rest_check{1,inside_loop}(3,1) = fst_mv_oj_coor(2,1);
                                                                                                                                prevent_rest_check{1,inside_loop}(3,2) = fst_mv_oj_coor(2,2);  
                                                                                                                                prevent_rest_check{1,inside_loop}(4,1) = fst_mv_oj_coor(3,1);
                                                                                                                                prevent_rest_check{1,inside_loop}(4,2) = fst_mv_oj_coor(3,2);  
                                                                                                                                prevent_rest_check{1,inside_loop}(5,1) = fst_mv_oj_coor(4,1);
                                                                                                                                prevent_rest_check{1,inside_loop}(5,2) = fst_mv_oj_coor(4,2);

                                                                                                                                if(stop_search2==1)
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start+5,1) = fst_mv_oj_coor(5,1);
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start+5,2) = fst_mv_oj_coor(5,2);
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start+5,3) = inside_loop;

                                                                                                                                    prevent_rest_check{1,inside_loop}(6,1)= fst_mv_oj_coor(5,1);
                                                                                                                                    prevent_rest_check{1,inside_loop}(6,2)= fst_mv_oj_coor(5,2);
                                                                                                                                end
                                                                                                                                stop_search=1; %prevent duplicate trajectories
                                                                                                                                cont_rest_th(inside_loop,1)=dist_final_th;
                                                                                                                                cont_angle_th(inside_loop,1)=angle_final_th;
                                                                                                                                if(abs(angle_final_th) > tangent_control)
                                                                                                                                    cont_angle_cmp(inside_loop,1)=abs(angle_rest_th)-abs(angle_final_th);
                                                                                                                                else
                                                                                                                                    cont_angle_cmp(inside_loop,1)=angle_rest_th-angle_final_th;
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
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end %for find_secd_loop
                                    end% for check_secd 
                                end
                            end
                        end                
                    end
                end %end for try_loop
            end %end for check_dist    
        end
        if( fast_obj_start > 1)
           if( stop_search == 0 )
               prevent_rest_check{1,inside_loop}(1,1) = nan;
               prevent_rest_check{1,inside_loop}(1,2) = nan;
           end
       end
    end
end
