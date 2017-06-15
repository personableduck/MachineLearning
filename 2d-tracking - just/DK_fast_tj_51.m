%New trajectories5 by DUCK-HA HWANG

%------solve diff holo dislocation, change xyz coordinate information.

constant_distance_rate=0.6; % endure 40% difference
threshold_diff_dislocation= 1+constant_distance_rate; %*****check it!! important value (1+1*40%) for combining close coordinates
control_vecter_angle = 40; %angle control
multilple_short_angle = 1.25; %find wiggling moving for slow one
increase_distance_threshold=0.8; %make sure about linear moving for fast
combine_limitation=1; % combine xyz noise distance
tangent_control=70;%tangent absolute value 60 degree
priority_degree_rate=30; %30

fps=3; %fps means frame per second
distance_threshold_fast = (140/fps)/pixelsize; %******threshold object tracking. need to control. maximum ditance to distinguish as a same object
%unit is um. // I limited that the fasted sperm's speed is 120 um/sec.
%fps means frame per second

%distance_threshold_skip=4.5;
distance_threshold_skip=(30/fps)/pixelsize; %<----30 before

decide_slow=(40/fps)/pixelsize;
decide_slow_angle=(60/fps)/pixelsize;%criteria for slow moving
%decide_slow=11; %distancce that ovelap happen on the differ_holo_stack

% %-----------combine xyz
% 
% for combine_compare_loop=1:iter-2
%     
%     behind_compare=combine_compare_loop+1;
%     
%     maxsize_cmb_cp=size(xyz{1,combine_compare_loop}); %read size
%     
%     maxsize_cmb_cp2=size(xyz{1,behind_compare}); %read size
% 
%     
%     for inft_loop=1:maxsize_cmb_cp(1,1)
%         for bhid_loop=1:maxsize_cmb_cp2(1,1)
% 
%             xx_shrt = xyz{1,combine_compare_loop}(inft_loop,1) - xyz{1,behind_compare}(bhid_loop,1); % comparing x coordinates
%             yy_shrt = xyz{1,combine_compare_loop}(inft_loop,2) - xyz{1,behind_compare}(bhid_loop,2); % comparing y coordinates
%         
%             distance_shrt_limit= sqrt((xx_shrt^2)+(yy_shrt^2));
% 
%             if (distance_shrt_limit <= combine_limitation)         
%                 xyz{1,behind_compare}(bhid_loop,1) = xyz{1,combine_compare_loop}(inft_loop,1);
%                 xyz{1,behind_compare}(bhid_loop,2) = xyz{1,combine_compare_loop}(inft_loop,2);
%                 
%             end
%         end
%     end
%     
% end

% %------change xyz
% fast_obj_start=1; %make a loop
% 
% %----------make a combined matrix
% 
% sum3_xyz=xyz{1,fast_obj_start};
% 
% for sum_fst_loop=fast_obj_start+1:iter-1
% 
%     comp_xyz_new=0;
%     comp_xyz_new2=0;
%     
%     comp_xyz_new=xyz{1,sum_fst_loop-1};
%     comp_xyz_new2=xyz{1,sum_fst_loop};
% 
%     size_cmp_n=size(comp_xyz_new);
%     
%     for loop_cmp_del=1:size_cmp_n(1,1)
%     
%         [row_del_x]=find(comp_xyz_new2(:,1) == comp_xyz_new(loop_cmp_del,1));
%         [row_del_y]=find(comp_xyz_new2(:,2) == comp_xyz_new(loop_cmp_del,2));
%     
%         rr_ss=size(row_del_x);
%         for check_loop_dlt=1:rr_ss
%             [row_del_f,col_del_f,check_del_f]=find(row_del_y == row_del_x(check_loop_dlt));
%             
%             if(check_del_f ~= 0)
%                 comp_xyz_new2(row_del_x(check_loop_dlt),:)=nan;
%             end
%         end
%         
%     end
%     
%     sum3_xyz=[sum3_xyz; comp_xyz_new2];
% end
%%--------------------------------make starting point
% xyz_1=0;
% xyz_1=xyz;
% for xyz_change_loop=1:iter-2
%     
%     comp_xyz_new=0;
%     comp_xyz_new2=0;
%     rsh_number=0;
%     researching_special=0;
% 
%     comp_xyz_new=xyz{1,xyz_change_loop};
%     comp_xyz_new2=xyz{1,xyz_change_loop+1};
% 
%     size_cmp_n=size(comp_xyz_new2);
% 
%     for loop_cmp_del=1:size_cmp_n(1,1)
% 
%         [row_del_x]=find(comp_xyz_new(:,1) == comp_xyz_new2(loop_cmp_del,1));
%         [row_del_y]=find(comp_xyz_new(:,2) == comp_xyz_new2(loop_cmp_del,2));
% 
%         rr_ss=size(row_del_x);
%         for check_loop_dlt=1:rr_ss
%             [row_del_f,col_del_f,check_del_f]=find(row_del_y == row_del_x(check_loop_dlt));
% 
%             if(check_del_f ~= 0)
%                 rsh_number=rsh_number+1;
%                 researching_special(rsh_number,1)=xyz{1,xyz_change_loop}(row_del_x(check_loop_dlt),1);
%                 researching_special(rsh_number,2)=xyz{1,xyz_change_loop}(row_del_x(check_loop_dlt),2);
%                 researching_special(rsh_number,3)=row_del_x(check_loop_dlt);
%             end
%         end
%     end
% 
%     maxsize_rsch_xyz1=size(researching_special); %read size
% 
%     for find_start_loop=1:maxsize_rsch_xyz1(1,1)-1
%         [row_search_strt,col_search_strt,check_search_strt]=find( researching_special(:,1) < researching_special(find_start_loop,1) + distance_threshold_fast & researching_special(:,1) > researching_special(find_start_loop,1) - distance_threshold_fast); 
%         maxsize_check_search_strt=size(check_search_strt);
%         if(maxsize_check_search_strt(1,1) > 1)
%             counting_strt=0;
%             distance_strt=0;
%             for row_search_strt_loop=1:maxsize_check_search_strt(1,1)
%                 srt_xx=researching_special((row_search_strt(row_search_strt_loop)),1)-researching_special(find_start_loop,1);
%                 srt_yy=researching_special((row_search_strt(row_search_strt_loop)),2)-researching_special(find_start_loop,2);
%                 if(srt_xx == 0)
%                     srt_xx=nan;
%                 end
%                 if(srt_yy == 0)
%                     srt_yy=nan;
%                 end
%                 counting_strt=counting_strt+1;
%                 distance_strt(counting_strt,1)=sqrt(srt_xx^2+srt_yy^2);
%                 distance_strt(counting_strt,2)=row_search_strt_loop;
%             end
%             [row_distrt,col_distrt,check_distrt]= find(distance_strt(:,1) <= distance_threshold_fast);
%             if(check_distrt == 1)           
%                 left_only=nan;
%                 left_one_moving=0;
%                 left_one_moving(1,1)=min_intensities_xyz{1,xyz_change_loop}(researching_special(find_start_loop,3));
%                 left_one_moving(1,2)=researching_special(find_start_loop,3);
%                 maxsize_row_distrt=size(row_distrt);
%                 for row_distrt_loop=1:maxsize_row_distrt(1,1)
%                     left_one_moving(row_distrt_loop+1,1)=min_intensities_xyz{1,xyz_change_loop}(researching_special(row_search_strt(distance_strt(row_distrt(row_distrt_loop),2)),3));
%                     left_one_moving(row_distrt_loop+1,2)=researching_special(row_search_strt(distance_strt(row_distrt(row_distrt_loop),2)),3);
%                 end
%                 left_only=max(left_one_moving(:,1));
%                 maxsize_left_one_moving=size(left_one_moving);
%                 for left_one_moving_loop=1:maxsize_left_one_moving(1,1)
%                     if(left_one_moving(left_one_moving_loop,1) ~= left_only)
%                         xyz_1{1,xyz_change_loop}(left_one_moving(left_one_moving_loop,2),:) = nan;
%                     end
%                 end
%             else
%                 [row_search_strt2,col_search_strt2,check_search_strt2]=find( xyz{1,xyz_change_loop}(:,1) < xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1) + distance_threshold_fast & xyz{1,xyz_change_loop}(:,1) > xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1) - distance_threshold_fast);
%                 maxsize_check_search_strt2=size(check_search_strt2);
%                 if(maxsize_check_search_strt2(1,1) > 1)
%                     counting_strt=0;
%                     distance_strt=0;
%                     check_distrt=0;
%                     row_distrt=0;
%                     col_distrt=0;
%                     for row_search_strt_loop=1:maxsize_check_search_strt2(1,1)
%                         srt_xx=xyz{1,xyz_change_loop}(row_search_strt2(row_search_strt_loop),1)-xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1);
%                         srt_yy=xyz{1,xyz_change_loop}(row_search_strt2(row_search_strt_loop),2)-xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),2);
%                         if(srt_xx == 0)
%                             srt_xx=nan;
%                         end
%                         if(srt_yy == 0)
%                             srt_yy=nan;
%                         end
%                         counting_strt=counting_strt+1;
%                         distance_strt(counting_strt,1)=sqrt(srt_xx^2+srt_yy^2);
%                         distance_strt(counting_strt,2)=row_search_strt_loop;
%                     end
%                     [row_distrt,col_distrt,check_distrt]= find(distance_strt(:,1) <= distance_threshold_fast);
%                     if(check_distrt == 1)
%                         xyz_1{1,xyz_change_loop}(researching_special(find_start_loop,3),:)=nan;
%                     end
%                 end
%             end
%         else
%             [row_search_strt3,col_search_strt3,check_search_strt3]=find( xyz{1,xyz_change_loop}(:,1) < xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1) + distance_threshold_fast & xyz{1,xyz_change_loop}(:,1) > xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1) - distance_threshold_fast);
%             maxsize_check_search_strt3=size(check_search_strt3);
%             if(maxsize_check_search_strt3(1,1) > 1)
%                 counting_strt=0;
%                 distance_strt=0;
%                 check_distrt=0;
%                 row_distrt=0;
%                 col_distrt=0;
%                 for row_search_strt_loop=1:maxsize_check_search_strt3(1,1)
%                     srt_xx=xyz{1,xyz_change_loop}(row_search_strt3(row_search_strt_loop),1)-xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1);
%                     srt_yy=xyz{1,xyz_change_loop}(row_search_strt3(row_search_strt_loop),2)-xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),2);
%                     if(srt_xx == 0)
%                         srt_xx=nan;
%                     end
%                     if(srt_yy == 0)
%                         srt_yy=nan;
%                     end
%                     counting_strt=counting_strt+1;
%                     distance_strt(counting_strt,1)=sqrt(srt_xx^2+srt_yy^2);
%                     distance_strt(counting_strt,2)=row_search_strt_loop;
%                 end
%                 [row_distrt,col_distrt,check_distrt]= find(distance_strt(:,1) <= distance_threshold_fast);
%                 if(check_distrt == 1)
%                     xyz_1{1,xyz_change_loop}(researching_special(find_start_loop,3),:)=nan;
%                 end
%             end
%         end
%     end
% end

%****************************************main trajectory***************************

max_frame_size=fix((iter-4)/2);

frame_number=2*max_frame_size+3

object_number_count=0;
cont_rest_th=0;
cont_angle_th=0;
fast_total_coor={0};
angle_information=0;

for rest_start_searching= 1:max_frame_size

    fast_obj_start=2*rest_start_searching-1;
   
    if(fast_obj_start == 1)
        maxsize_fast_xyz=size(xyz{1,fast_obj_start}); %read size
        change_inside_loop = maxsize_fast_xyz(1,1);
    else
        change_inside_loop = object_number_count;
    end
    
    for inside_loop=1:change_inside_loop
        stop_search=0;
        distance_fast_mat=0;    
        maxsize_fast2_xyz=size(xyz{1,fast_obj_start+1}); %read size    
        
        for fast_second_mt=1:maxsize_fast2_xyz(1,1) %for fast sperm
            if(fast_obj_start == 1)
                xx_1 = xyz{1,fast_obj_start+1}(fast_second_mt,1) - xyz_1{1,fast_obj_start}(inside_loop,1); % comparing x coordinates   
                yy_1 = xyz{1,fast_obj_start+1}(fast_second_mt,2) - xyz_1{1,fast_obj_start}(inside_loop,2); % comparing y coordinates
                
                distance_fast_mat(fast_second_mt,1) = sqrt((xx_1^2)+(yy_1^2));
            else
                estm_ft_size=size(fast_total_coor{1,inside_loop});
                if( estm_ft_size(1,1) >= fast_obj_start )
                    xx_1 = xyz{1,fast_obj_start+1}(fast_second_mt,1) - fast_total_coor{1,inside_loop}(fast_obj_start,1); % comparing x coordinates   
                    yy_1 = xyz{1,fast_obj_start+1}(fast_second_mt,2) - fast_total_coor{1,inside_loop}(fast_obj_start,2); % comparing y coordinates
                else     
                    xx_1=nan;
                    yy_1=nan;
                end                
                distance_fast_mat(fast_second_mt,1) = sqrt((xx_1^2)+(yy_1^2));
            end
        end
        if(fast_obj_start == 1) 
            [row_dist,col_dist,check_dist]=find(distance_fast_mat <= distance_threshold_fast);     
        else
            if( cont_rest_th(inside_loop,1) <= decide_slow)
                [row_dist,col_dist,check_dist]=find(distance_fast_mat <= distance_threshold_fast);            
            else
                [row_dist,col_dist,check_dist]=find(distance_fast_mat < cont_rest_th(inside_loop,1)*threshold_diff_dislocation & distance_fast_mat <= distance_threshold_fast);
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
                    priority_angle_order(try_loop,2)=try_loop;%give priority for smaller distance!!
                    [temp_pr1,ord_pr1] = sort(priority_angle_order(:,1));
                    priority_angle_order = priority_angle_order(ord_pr1,:);
                end
            else 
                for try_loop = 1: row_dist_size(1,1) 
                    angle_fast_test = atan((xyz{1,fast_obj_start+1}(row_dist(try_loop),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))/(xyz{1,fast_obj_start+1}(row_dist(try_loop),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi; % unit is degree   
                    priority_angle_order(try_loop,1)=abs(cont_angle_th(inside_loop,1)-angle_fast_test)/priority_degree_rate + distance_fast_mat(row_dist(try_loop));
                    priority_angle_order(try_loop,2)=try_loop;%give priority for smaller angle!!
                    if(abs(cont_angle_th(inside_loop,1))>=tangent_control)
                        priority_angle_order(try_loop,1)=abs(abs(cont_angle_th(inside_loop,1))-abs(angle_fast_test))/priority_degree_rate + distance_fast_mat(row_dist(try_loop));
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
                        angle_fast = atan((xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2)-xyz{1,fast_obj_start}(inside_loop,2))/(xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)-xyz{1,fast_obj_start}(inside_loop,1)))*180/pi; % unit is degree
                        dist_thresh_sec = distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)));
                        if(angle_fast <= 45 & angle_fast >= -45) %vetorc
                            vector=1; % x vetorc
                            diff_vector = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1) - xyz{1,fast_obj_start}(inside_loop,1);
                            posit_nega=sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                        else
                            vector=0; % y vetorc
                            diff_vector = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2) - xyz{1,fast_obj_start}(inside_loop,2); 
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
                            diff_vector2= xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1);    
                            posit_nega2=sign(diff_vector2);
                        elseif (vector == 0)            
                            diff_vector2= xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2);
                            posit_nega2=sign(diff_vector2);          
                        end
                        angle_fast_test=atan((xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))/(xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi;
                        if( cont_rest_th(inside_loop,1) > decide_slow_angle)           
                            if(abs(cont_angle_th(inside_loop,1)) > tangent_control | abs(angle_fast_test) > tangent_control)              
                                if(abs(cont_angle_th(inside_loop,1))+control_vecter_angle > abs(angle_fast_test) & abs(cont_angle_th(inside_loop,1))-control_vecter_angle < abs(angle_fast_test)) % unit is degree
                                    pass_numfil=1;
                                end
                            else
                                if(cont_angle_th(inside_loop,1)+control_vecter_angle > angle_fast_test & cont_angle_th(inside_loop,1)-control_vecter_angle < angle_fast_test) % unit is degree
                                pass_numfil=1;
                                end
                            end    
                        else
                            if(abs(cont_angle_th(inside_loop,1))+control_vecter_angle*multilple_short_angle > abs(angle_fast_test) & abs(cont_angle_th(inside_loop,1))-control_vecter_angle*multilple_short_angle < abs(angle_fast_test)) % unit is degree
                                pass_numfil=1;
                            end
                        end
                    end        
                    if(posit_nega == posit_nega2) 
                        pass_value=1;
                        increase_check=1;
                        if(fast_obj_start == 1)                         
                            pass_value=0; 
                        else
                            if((distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)))+cont_rest_th(inside_loop,1))*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start-1,1))^2+(xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start-1,2))^2 ))          
                                increase_check=0; 
                            else   
                                increase_check=1;
                            end

                            if(increase_check == 0)                
                                if(abs(cont_angle_th(inside_loop,1)) > tangent_control)
                                    if( abs(cont_angle_cmp(inside_loop,1)) + control_vecter_angle > abs(abs(cont_angle_th(inside_loop,1)) - abs(angle_fast_test)) & abs(cont_angle_cmp(inside_loop,1)) - control_vecter_angle < abs(abs(cont_angle_th(inside_loop,1)) - abs(angle_fast_test)) )
                                        pass_value=0; 
                                    else
                                        pass_value=1;    
                                    end
                                else
                                    if( abs(cont_angle_cmp(inside_loop,1)) + control_vecter_angle > abs(cont_angle_th(inside_loop,1) - angle_fast_test) & abs(cont_angle_cmp(inside_loop,1)) - control_vecter_angle < abs(cont_angle_th(inside_loop,1) - angle_fast_test) )
                                        pass_value=0; 
                                    else
                                        pass_value=1;    
                                    end
                                end
                            end
                        end
                        if(pass_value == 0)       
                            if(pass_numfil == 1)
                                if(fast_obj_start > 1) 
                                    angle_fast = atan((xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))/(xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi;
                                    dist_thresh_sec = distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)));                                
                                    % finish third check
                                end   
                                distance_repeat_rest=0;
                                maxsize_fast_xyz2=size(xyz{1,fast_obj_start+2}); %read size 

                                for dist_cmp_loop = 1: maxsize_fast_xyz2(1,1)
                                    xx_2 = xyz{1,fast_obj_start+2}(dist_cmp_loop,1) - xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1);                              
                                    yy_2 = xyz{1,fast_obj_start+2}(dist_cmp_loop,2) - xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2);

                                    distance_repeat_rest(dist_cmp_loop,1) = sqrt((xx_2^2)+(yy_2^2));
                                end
                                if(dist_thresh_sec <= decide_slow)
                                    [row_secd,col_secd,check_secd]=find(distance_repeat_rest < distance_threshold_fast);
                                else
                                    [row_secd,col_secd,check_secd]=find(distance_repeat_rest < dist_thresh_sec*threshold_diff_dislocation & distance_repeat_rest <= distance_threshold_fast);
                                end
                                if (check_secd == 1)
                                    row_secd_size=size(row_secd); %read size  
                                    angle_secd_fast=0;
                                    dist_thresh_third=0;
                                    priority_angle_order2=0;

                                    for find_secd_loop = 1: row_secd_size(1,1)
                                        angle_secd_fast_test = atan((xyz{1,fast_obj_start+2}(row_secd(find_secd_loop),2)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2))/(xyz{1,fast_obj_start+2}(row_secd(find_secd_loop),1)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)))*180/pi;
                                        priority_angle_order2(find_secd_loop,1)=abs(angle_fast-angle_secd_fast_test)/priority_degree_rate + distance_repeat_rest(row_secd(find_secd_loop));
                                        priority_angle_order2(find_secd_loop,2)=find_secd_loop;%give priority for smaller angle!!
                                        if(abs(angle_fast)>=tangent_control)
                                            priority_angle_order2(find_secd_loop,1)=abs(abs(angle_fast)-abs(angle_secd_fast_test))/priority_degree_rate + distance_repeat_rest(row_secd(find_secd_loop));
                                        end
                                        [temp_pr1,ord_pr1] = sort(priority_angle_order2(:,1));
                                        priority_angle_order2 = priority_angle_order2(ord_pr1,:);
                                    end

                                    for find_secd_loop2 = 1: row_secd_size(1,1)
                                        pass_numfil=0;
                                        if(stop_search == 0)
                                            diff_vector2=0;
                                            posit_nega2=nan;
                                        
                                            if( fast_obj_start > 1)
                                                if(angle_fast <= 45 & angle_fast >= -45) %vetorc
                                                    vector=1; % x vetorc
                                                    diff_vector = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1);
                                                    posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                else 
                                                    vector=0; % y vetorc
                                                    diff_vector = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2);
                                                    posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                end
                                            end

                                            if(vector == 1)
                                                diff_vector2= xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1);
                                                posit_nega2=sign(diff_vector2);
                                            elseif (vector == 0)
                                                diff_vector2= xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2);                            
                                                posit_nega2=sign(diff_vector2);
                                            end

                                            if(posit_nega == posit_nega2)
                                                pass_value=1;
                                                increase_check=1;

                                                if(fast_obj_start == 1)
                                                    if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz{1,fast_obj_start}(inside_loop,1))^2+(xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz{1,fast_obj_start}(inside_loop,2))^2 ))
                                                        pass_value=0;
                                                    else
                                                        pass_value=1;
                                                    end
                                                else
                                                    if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1))^2+(xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))^2 ))
                                                        increase_check=0;
                                                    else
                                                        increase_check=1;
                                                    end
                                                    if(increase_check == 0)
                                                        angle_secd_fast_test = atan((xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2))/(xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)))*180/pi;
                                                        if(abs(angle_fast) > tangent_control)
                                                            if( abs(abs(cont_angle_th(inside_loop,1))-abs(angle_fast)) + control_vecter_angle > abs(abs(angle_fast) - abs(angle_secd_fast_test)) & abs(abs(cont_angle_th(inside_loop,1))-abs(angle_fast)) - control_vecter_angle < abs(abs(angle_fast) - abs(angle_secd_fast_test)) )
                                                                pass_value=0;
                                                            else
                                                                pass_value=1;
                                                            end 
                                                        else
                                                            if( abs(cont_angle_th(inside_loop,1)-angle_fast) + control_vecter_angle > abs(angle_fast - angle_secd_fast_test) & abs(cont_angle_th(inside_loop,1)-angle_fast) - control_vecter_angle < abs(angle_fast - angle_secd_fast_test) )
                                                                pass_value=0;
                                                            else
                                                                pass_value=1;
                                                            end    
                                                        end
                                                    end
                                                end
                                                if(pass_value == 0)
                                                    angle_secd_fast_cmp=atan((xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2))/(xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)))*180/pi; 

                                                    if(dist_thresh_sec > decide_slow_angle)
                                                        if(abs(angle_fast)>tangent_control | abs(angle_secd_fast_cmp) > tangent_control)                                                     
                                                            if(abs(angle_fast)+control_vecter_angle > abs(angle_secd_fast_cmp) & abs(angle_fast)-control_vecter_angle < abs(angle_secd_fast_cmp)) % unit is degree                       
                                                                pass_numfil=1;
                                                            end

                                                        else

                                                            if(angle_fast+control_vecter_angle > angle_secd_fast_cmp & angle_fast-control_vecter_angle < angle_secd_fast_cmp) % unit is degree  
                                                                pass_numfil=1;
                                                            end
                                                        end
                                                    else
                                                        if(abs(angle_fast)+control_vecter_angle*multilple_short_angle > abs(angle_secd_fast_cmp) & abs(angle_fast)-control_vecter_angle*multilple_short_angle < abs(angle_secd_fast_cmp)) % unit is degree                       
                                                            pass_numfil=1;
                                                        end
                                                    end

                                                    if(pass_numfil == 1) 
                                                        angle_secd_fast = atan((xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2))/(xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)))*180/pi;   
                                                        dist_thresh_third = distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)));                           
                                                        % finish third
                                                        distance_thrd_mat=0;
                                                        maxsize_fast_xyz3=size(xyz{1,fast_obj_start+3});

                                                        for dist_thrd_loop = 1: maxsize_fast_xyz3(1,1)                                   
                                                            xx_3 = xyz{1,fast_obj_start+3}(dist_thrd_loop,1) - xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1);
                                                            yy_3 = xyz{1,fast_obj_start+3}(dist_thrd_loop,2) - xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2);
                                                            distance_thrd_mat(dist_thrd_loop,1) = sqrt((xx_3^2)+(yy_3^2));                           
                                                        end

                                                        if(dist_thresh_third <= decide_slow)
                                                            [row_third,col_third,check_third]=find(distance_thrd_mat < distance_threshold_fast);
                                                        else
                                                            [row_third,col_third,check_third]=find(distance_thrd_mat < dist_thresh_third*threshold_diff_dislocation & distance_thrd_mat <= distance_threshold_fast);
                                                        end

                                                        if (check_third == 1)

                                                            row_third_size=size(row_third);
                                                            angle_rest_th=0;
                                                            dist_rest_th=0;                         
                                                            priority_angle_order3=0;

                                                            for find_rest_loop = 1: row_third_size(1,1)
                                                                angle_rest_th_test=atan((xyz{1,fast_obj_start+3}(row_third(find_rest_loop),2)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2))/(xyz{1,fast_obj_start+3}(row_third(find_rest_loop),1)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)))*180/pi;                                                           
                                                                priority_angle_order3(find_rest_loop,1)=abs(angle_secd_fast-angle_rest_th_test)/priority_degree_rate + distance_thrd_mat(row_third(find_rest_loop));
                                                                priority_angle_order3(find_rest_loop,2)=find_rest_loop;%give priority for smaller angle!!
                                                                if(abs(angle_secd_fast)>=tangent_control)
                                                                    priority_angle_order3(find_rest_loop,1)=abs(abs(angle_secd_fast)-abs(angle_rest_th_test))/priority_degree_rate + distance_thrd_mat(row_third(find_rest_loop));
                                                                end
                                                                [temp_pr1,ord_pr1] = sort(priority_angle_order3(:,1));
                                                                priority_angle_order3 = priority_angle_order3(ord_pr1,:);
                                                            end

                                                            for find_rest_loop2 = 1: row_third_size(1,1)
                                                                pass_numfil=0;  
                                                                if(stop_search == 0)
                                                                    diff_vector2=0;
                                                                    posit_nega2=nan;
                                                                    if(angle_secd_fast <= 45 & angle_secd_fast >= -45) %vetorc
                                                                        vector=1; % x vetorc
                                                                        diff_vector = xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1);
                                                                        posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                    else 
                                                                        vector=0; % y vetorc
                                                                        diff_vector = xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2);
                                                                        posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                    end
                                                                    if(vector == 1)       
                                                                        diff_vector2= xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1) - xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1);
                                                                        posit_nega2=sign(diff_vector2);                                           
                                                                    elseif(vector == 0)                                            
                                                                        diff_vector2= xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2) - xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2);                                            
                                                                        posit_nega2=sign(diff_vector2);
                                                                    end

                                                                    if(posit_nega == posit_nega2)
                                                                        pass_value=1;
                                                                        increase_check=1;
                                                                        if((distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)))+dist_thresh_third)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1))^2+(xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2))^2 ))
                                                                            increase_check=0;                                                     
                                                                        else            
                                                                            increase_check=1;
                                                                        end  
                                                                        if( increase_check == 0)
                                                                            angle_rest_th_test=atan((xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2))/(xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)))*180/pi;                                                                 
                                                                            if( abs(angle_secd_fast) > tangent_control)
                                                                                if( abs(abs(angle_fast)-abs(angle_secd_fast)) + control_vecter_angle > abs(abs(angle_secd_fast)-abs(angle_rest_th_test)) & abs(abs(angle_fast)-abs(angle_secd_fast)) + control_vecter_angle > abs(abs(angle_secd_fast)-abs(angle_rest_th_test)) )
                                                                                    pass_value=0;
                                                                                else
                                                                                    pass_value=1;  
                                                                                end
                                                                            else  
                                                                                if( abs(angle_fast-angle_secd_fast) + control_vecter_angle > abs(angle_secd_fast-angle_rest_th_test) & abs(angle_fast-angle_secd_fast) + control_vecter_angle > abs(angle_secd_fast-angle_rest_th_test) )
                                                                                    pass_value=0;
                                                                                else
                                                                                    pass_value=1;  
                                                                                end
                                                                            end
                                                                        end
                                                                        if(pass_value == 0)                                                      
                                                                            if( dist_thresh_third > decide_slow_angle)
                                                                                if(abs(angle_secd_fast)>tangent_control | abs(angle_rest_th_test)>tangent_control)    
                                                                                    if(abs(angle_secd_fast)+control_vecter_angle > abs(angle_rest_th_test) & abs(angle_secd_fast)-control_vecter_angle < abs(angle_rest_th_test))  
                                                                                        pass_numfil=1;
                                                                                    end
                                                                                else
                                                                                    if(angle_secd_fast+control_vecter_angle > angle_rest_th_test & angle_secd_fast-control_vecter_angle < angle_rest_th_test)  
                                                                                        pass_numfil=1;
                                                                                    end
                                                                                end
                                                                            else
                                                                                if(abs(angle_secd_fast)+control_vecter_angle*multilple_short_angle > abs(angle_rest_th_test) & abs(angle_secd_fast)-control_vecter_angle*multilple_short_angle < abs(angle_rest_th_test))  
                                                                                    pass_numfil=1;
                                                                                end  
                                                                            end

                                                                            if( pass_numfil == 1)
                                                                                angle_rest_th=atan((xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2))/(xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)))*180/pi;                 
                                                                                dist_rest_th=distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)));
                                                                                distance_four_mat=0;
                                                                                maxsize_fast_xyz4=size(xyz{1,fast_obj_start+4});
                                                                                for dist_four_loop = 1: maxsize_fast_xyz4(1,1)
                                                                                    xx_4 = xyz{1,fast_obj_start+4}(dist_four_loop,1) - xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1);    
                                                                                    yy_4 = xyz{1,fast_obj_start+4}(dist_four_loop,2) - xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2);     
                                                                                    distance_four_mat(dist_four_loop,1) = sqrt((xx_4^2)+(yy_4^2));                           
                                                                                end
                                                                                if(dist_rest_th <= decide_slow)     
                                                                                    [row_four,col_four,check_four]=find(distance_four_mat < distance_threshold_fast);                              
                                                                                else
                                                                                    [row_four,col_four,check_four]=find(distance_four_mat < dist_rest_th*threshold_diff_dislocation & distance_four_mat <= distance_threshold_fast);   
                                                                                end

                                                                                if (check_four == 1)
                                                                                    row_four_size=size(row_four);
                                                                                    angle_final_th=0;
                                                                                    dist_final_th=0; 
                                                                                    priority_angle_order4=0;

                                                                                    for find_four_loop = 1: row_four_size(1,1)
                                                                                        angle_final_th_test=atan((xyz{1,fast_obj_start+4}(row_four(find_four_loop),2)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2))/(xyz{1,fast_obj_start+4}(row_four(find_four_loop),1)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)))*180/pi;
                                                                                        priority_angle_order4(find_four_loop,1)=abs(angle_rest_th-angle_final_th_test)/priority_degree_rate + distance_four_mat(row_four(find_four_loop));
                                                                                        priority_angle_order4(find_four_loop,2)=find_four_loop;%give priority for smaller angle!!
                                                                                        if(abs(angle_rest_th)>=tangent_control)
                                                                                            priority_angle_order4(find_four_loop,1)=abs(abs(angle_rest_th)-abs(angle_final_th_test))/priority_degree_rate + distance_four_mat(row_four(find_four_loop));
                                                                                        end
                                                                                        [temp_pr1,ord_pr1] = sort(priority_angle_order4(:,1));
                                                                                        priority_angle_order4 = priority_angle_order4(ord_pr1,:);
                                                                                    end

                                                                                    for find_four_loop2 = 1: row_four_size(1,1)                                                          
                                                                                        pass_numfil=0;
                                                                                        if(stop_search == 0)  
                                                                                            diff_vector2=0;
                                                                                            posit_nega2=nan;
                                                                                            if(angle_rest_th <= 45 & angle_rest_th >= -45) %vetorc
                                                                                                vector=1; % x vetorc
                                                                                                diff_vector = xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1);
                                                                                                posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                            else 
                                                                                                vector=0; % y vetorc
                                                                                                diff_vector = xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2);
                                                                                                posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                            end
                                                                                            if(vector == 1)       
                                                                                                diff_vector2= xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1) - xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1);      
                                                                                                posit_nega2=sign(diff_vector2);        
                                                                                            elseif(vector == 0)                                            
                                                                                                diff_vector2= xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2) - xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2);                                            
                                                                                                posit_nega2=sign(diff_vector2);
                                                                                            end

                                                                                            if(posit_nega == posit_nega2)
                                                                                                pass_value=1;
                                                                                                increase_check=1;
                                                                                                if((distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)))+dist_rest_th)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1))^2+(xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2))^2 ))       
                                                                                                    increase_check=0;    
                                                                                                else
                                                                                                    increase_check=1;      
                                                                                                end
                                                                                                if( increase_check == 0)                                                                                      
                                                                                                    angle_final_th_test=atan((xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2))/(xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)))*180/pi;
                                                                                                    if(abs(angle_rest_th) > tangent_control)
                                                                                                        if( abs(abs(angle_secd_fast)-abs(angle_rest_th)) + control_vecter_angle > abs(abs(angle_rest_th)-abs(angle_final_th_test)) & abs(abs(angle_secd_fast)-abs(angle_rest_th)) - control_vecter_angle < abs(abs(angle_rest_th)-abs(angle_final_th_test)) )
                                                                                                            pass_value=0;    
                                                                                                        else
                                                                                                            pass_value=1;      
                                                                                                        end
                                                                                                    else 
                                                                                                        if( abs(angle_secd_fast-angle_rest_th) + control_vecter_angle > abs(angle_rest_th-angle_final_th_test) & abs(angle_secd_fast-angle_rest_th) - control_vecter_angle < abs(angle_rest_th-angle_final_th_test) )
                                                                                                            pass_value=0;    
                                                                                                        else
                                                                                                            pass_value=1;      
                                                                                                        end
                                                                                                    end
                                                                                                end
                                                                                                if(pass_value == 0)                                                                                             
                                                                                                    angle_final_th_cmp=atan((xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2))/(xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)))*180/pi;
                                                                                                    if( dist_rest_th > decide_slow_angle)
                                                                                                        if(abs(angle_rest_th)>tangent_control | abs(angle_final_th_cmp) >tangent_control)
                                                                                                            if(abs(angle_rest_th)+control_vecter_angle > abs(angle_final_th_cmp) & abs(angle_rest_th)-control_vecter_angle < abs(angle_final_th_cmp))
                                                                                                                pass_numfil=1;
                                                                                                            end
                                                                                                        else
                                                                                                            if(angle_rest_th+control_vecter_angle > angle_final_th_cmp & angle_rest_th-control_vecter_angle < angle_final_th_cmp)
                                                                                                                pass_numfil=1;
                                                                                                            end
                                                                                                        end
                                                                                                    else  
                                                                                                        if(abs(angle_rest_th)+control_vecter_angle*multilple_short_angle > abs(angle_final_th_cmp) & abs(angle_rest_th)-control_vecter_angle*multilple_short_angle < abs(angle_final_th_cmp))
                                                                                                            pass_numfil=1;
                                                                                                        end
                                                                                                    end

                                                                                                    if( pass_numfil == 1)
                                                                                                        angle_final_th=atan((xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2))/(xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)))*180/pi;                                    
                                                                                                        dist_final_th=distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)));

                                                                                                        if(fast_obj_start == 1)
                                                                                                            object_number_count = object_number_count+1;
                                                                                                            fst_mv_oj_coor=0;

                                                                                                            fst_mv_oj_coor(1,1) = xyz{1,fast_obj_start}(inside_loop,1);                              
                                                                                                            fst_mv_oj_coor(1,2) = xyz{1,fast_obj_start}(inside_loop,2);          
                                                                                                            fst_mv_oj_coor(1,3) = object_number_count;

                                                                                                            fst_mv_oj_coor(2,1) = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1);
                                                                                                            fst_mv_oj_coor(2,2) = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2);       
                                                                                                            fst_mv_oj_coor(2,3) = object_number_count;

                                                                                                            fst_mv_oj_coor(3,1) = xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1);
                                                                                                            fst_mv_oj_coor(3,2) = xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2);          
                                                                                                            fst_mv_oj_coor(3,3) = object_number_count;

                                                                                                            fst_mv_oj_coor(4,1) = xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1);
                                                                                                            fst_mv_oj_coor(4,2) = xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2);          
                                                                                                            fst_mv_oj_coor(4,3) = object_number_count;

                                                                                                            fst_mv_oj_coor(5,1) = xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1);
                                                                                                            fst_mv_oj_coor(5,2) = xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2);
                                                                                                            fst_mv_oj_coor(5,3) = object_number_count;

                                                                                                            if( abs(angle_final_th) > tangent_control)
                                                                                                                angle_information(object_number_count,1)=(abs(abs(angle_fast)-abs(angle_secd_fast))+abs(abs(angle_secd_fast)-abs(angle_rest_th))+abs(abs(angle_rest_th)-abs(angle_final_th)))/priority_degree_rate + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) );
                                                                                                            else
                                                                                                                angle_information(object_number_count,1)=(abs(angle_fast-angle_secd_fast)+abs(angle_secd_fast-angle_rest_th)+abs(angle_rest_th-angle_final_th))/priority_degree_rate + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) );
                                                                                                            end
                                                                                                            distance_six_mat=0;
                                                                                                            maxsize_fast_xyz6=size(xyz{1,fast_obj_start+4});
                                                                                                            for dist_six_loop = 1: maxsize_fast_xyz6(1,1)                                   
                                                                                                                xx_6 = xyz{1,fast_obj_start+4}(dist_six_loop,1) - xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1);
                                                                                                                yy_6 = xyz{1,fast_obj_start+4}(dist_six_loop,2) - xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2);
                                                                                                                distance_six_mat(dist_six_loop,1) = sqrt((xx_6^2)+(yy_6^2));                           
                                                                                                            end

                                                                                                            if(dist_final_th <= decide_slow)
                                                                                                                [row_six,col_six,check_six]=find(distance_six_mat < distance_threshold_fast);
                                                                                                            else
                                                                                                                [row_six,col_six,check_six]=find(distance_six_mat < dist_final_th*threshold_diff_dislocation & distance_six_mat <= distance_threshold_fast);
                                                                                                            end

                                                                                                            if (check_six == 1)
                                                                                                                row_six_size=size(row_six);
                                                                                                                angle_six_th=0;
                                                                                                                dist_six_th=0;                         
                                                                                                                priority_angle_order6=0;
                                                                                                                for find_six_loop = 1: row_six_size(1,1)
                                                                                                                    angle_six_th_test=atan((xyz{1,fast_obj_start+4}(row_six(find_six_loop),2)-xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2))/(xyz{1,fast_obj_start+4}(row_six(find_six_loop),1)-xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2)))*180/pi;                                                           
                                                                                                                    priority_angle_order6(find_six_loop,1)=abs(angle_final_th-angle_six_th_test)/priority_degree_rate + distance_six_mat(row_six(find_six_loop));
                                                                                                                    priority_angle_order6(find_six_loop,2)=find_six_loop;%give priority for smaller angle!!
                                                                                                                    if(abs(angle_final_th)>=tangent_control)
                                                                                                                        priority_angle_order6(find_six_loop,1)=abs(abs(angle_final_th)-abs(angle_six_th_test))/priority_degree_rate + distance_six_mat(row_six(find_six_loop));
                                                                                                                    end
                                                                                                                    [temp_pr1,ord_pr1] = sort(priority_angle_order6(:,1));
                                                                                                                    priority_angle_order6 = priority_angle_order6(ord_pr1,:);
                                                                                                                end
                                                                                                                for find_six_loop2 = 1: row_six_size(1,1)
                                                                                                                    pass_numfil=0;  
                                                                                                                    if(stop_search == 0)
                                                                                                                        diff_vector2=0;
                                                                                                                        posit_nega2=nan;
                                                                                                                        if(angle_final_th <= 45 & angle_final_th >= -45) %vetorc
                                                                                                                            vector=1; % x vetorc
                                                                                                                            diff_vector = xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1);
                                                                                                                            posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                        else 
                                                                                                                            vector=0; % y vetorc
                                                                                                                            diff_vector = xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2);
                                                                                                                            posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                        end
                                                                                                                        
                                                                                                                        if(vector == 1)       
                                                                                                                            diff_vector2= xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),1) - xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1);
                                                                                                                            posit_nega2=sign(diff_vector2);                                           
                                                                                                                        elseif(vector == 0)                                            
                                                                                                                            diff_vector2= xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),2) - xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2);
                                                                                                                            posit_nega2=sign(diff_vector2);
                                                                                                                        end
                                                                                                                        if(posit_nega == posit_nega2)
                                                                                                                            pass_value=1;
                                                                                                                            increase_check=1; 
                                                                                                                            if((distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)))+dist_final_th)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),1)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1))^2+(xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),2)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2))^2 ))
                                                                                                                                increase_check=0;                                                     
                                                                                                                            else            
                                                                                                                                increase_check=1;
                                                                                                                            end  
                                                                                                                            if( increase_check == 0)
                                                                                                                                angle_six_th_test=atan((xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),2)-xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2))/(xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),1)-xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1)))*180/pi;                                                                 
                                                                                                                                if( abs(angle_secd_fast) > tangent_control)
                                                                                                                                    if( abs(abs(angle_rest_th)-abs(angle_final_th)) + control_vecter_angle > abs(abs(angle_final_th)-abs(angle_six_th_test)) & abs(abs(angle_rest_th)-abs(angle_final_th)) - control_vecter_angle < abs(abs(angle_final_th)-abs(angle_six_th_test)) )
                                                                                                                                        pass_value=0;
                                                                                                                                    else
                                                                                                                                        pass_value=1;  
                                                                                                                                    end
                                                                                                                                else  
                                                                                                                                    if( abs(angle_rest_th-angle_final_th) + control_vecter_angle > abs(angle_final_th-angle_six_th_test) & abs(angle_rest_th-angle_final_th) + control_vecter_angle > abs(angle_final_th-angle_six_th_test) )
                                                                                                                                        pass_value=0;
                                                                                                                                    else
                                                                                                                                        pass_value=1;  
                                                                                                                                    end
                                                                                                                                end
                                                                                                                            end
                                                                                                                            if(pass_value == 0)                                                     
                                                                                                                                if( dist_final_th > decide_slow_angle)
                                                                                                                                    if(abs(angle_final_th)>tangent_control | abs(angle_six_th_test)>tangent_control)    
                                                                                                                                        if(abs(angle_final_th)+control_vecter_angle > abs(angle_six_th_test) & abs(angle_final_th)-control_vecter_angle < abs(angle_six_th_test))  
                                                                                                                                            pass_numfil=1;
                                                                                                                                        end
                                                                                                                                    else
                                                                                                                                        if(angle_final_th+control_vecter_angle > angle_six_th_test & angle_final_th-control_vecter_angle < angle_six_th_test)  
                                                                                                                                            pass_numfil=1;
                                                                                                                                        end
                                                                                                                                    end
                                                                                                                                else
                                                                                                                                    if(abs(angle_final_th)+control_vecter_angle*multilple_short_angle > abs(angle_six_th_test) & abs(angle_final_th)-control_vecter_angle*multilple_short_angle < abs(angle_six_th_test))  
                                                                                                                                        pass_numfil=1;
                                                                                                                                    end  
                                                                                                                                end

                                                                                                                                if( pass_numfil == 1)
                                                                                                                                    fst_mv_oj_coor(6,1) = xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),1);
                                                                                                                                    fst_mv_oj_coor(6,2) = xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),2);
                                                                                                                                    fst_mv_oj_coor(6,3) = object_number_count;
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
                                                                                                                cont_rest_th(object_number_count,1)=dist_thresh_third;
                                                                                                                cont_angle_th(object_number_count,1)=angle_secd_fast;
                                                                                                                if(abs(angle_secd_fast) > tangent_control)
                                                                                                                    cont_angle_cmp(object_number_count,1)=abs(angle_fast)-abs(angle_secd_fast);
                                                                                                                else
                                                                                                                    cont_angle_cmp(object_number_count,1)=angle_fast-angle_secd_fast;
                                                                                                                end
                                                                                                            else                
                                                                                                            %------prevent duplication--------
                                                                                                                maxsize_fstc=size(fast_total_coor);     
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
                                                                                                                    if(angle_information(duplication_number,1) < angle_information(object_number_count,1))
                                                                                                                        object_number_count= object_number_count-1; %duplication happen
                                                                                                                    else
                                                                                                                        fast_total_coor(1,duplication_number) = {0};
                                                                                                                        cont_rest_th(object_number_count,1)=dist_thresh_third;
                                                                                                                        cont_angle_th(object_number_count,1)=angle_secd_fast;
                                                                                                                        if(abs(angle_secd_fast) > tangent_control)
                                                                                                                            cont_angle_cmp(object_number_count,1)=abs(angle_fast)-abs(angle_secd_fast);
                                                                                                                        else
                                                                                                                            cont_angle_cmp(object_number_count,1)=angle_fast-angle_secd_fast;
                                                                                                                        end      
                                                                                                                        cont_rest_th(duplication_number,1)=cont_rest_th(object_number_count,1);
                                                                                                                        cont_angle_th(duplication_number,1)=cont_angle_th(object_number_count,1);  
                                                                                                                        cont_angle_cmp(duplication_number,1)=cont_angle_cmp(object_number_count,1);
                                                                                                                        object_number_count = object_number_count-1;
                                                                                                                        stop_search=1;

                                                                                                                        fst_mv_oj_coor(1,3) = duplication_number;
                                                                                                                        fst_mv_oj_coor(2,3) = duplication_number;
                                                                                                                        fst_mv_oj_coor(3,3) = duplication_number;
                                                                                                                        fst_mv_oj_coor(4,3) = duplication_number;
                                                                                                                        fast_total_coor(1,duplication_number) = {fst_mv_oj_coor}; %smaller angle is right trajectory
                                                                                                                    end
                                                                                                                else
                                                                                                                    fast_total_coor(1,object_number_count)={fst_mv_oj_coor};
                                                                                                                    stop_search=1; %prevent duplicate trajectories
                                                                                                                    cont_rest_th(object_number_count,1)=dist_thresh_third;
                                                                                                                    cont_angle_th(object_number_count,1)=angle_secd_fast;
                                                                                                                    if(abs(angle_secd_fast) > tangent_control)
                                                                                                                        cont_angle_cmp(object_number_count,1)=abs(angle_fast)-abs(angle_secd_fast);
                                                                                                                    else
                                                                                                                        cont_angle_cmp(object_number_count,1)=angle_fast-angle_secd_fast;
                                                                                                                    end
                                                                                                                end
                                                                                                            end
                                                                                                        else            
                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+1,1) = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1);
                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+1,2) = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2);       
                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+1,3) = inside_loop;

                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+2,1) = xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1);
                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+2,2) = xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2);          
                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+2,3) = inside_loop;

                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+3,1) = xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1);
                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+3,2) = xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2);          
                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+3,3) = inside_loop;

                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+4,1) = xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1);
                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+4,2) = xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2);
                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+4,3) = inside_loop;  

                                                                                                            distance_six_mat=0;
                                                                                                            maxsize_fast_xyz6=size(xyz{1,fast_obj_start+4});
                                                                                                            for dist_six_loop = 1: maxsize_fast_xyz6(1,1)                                   
                                                                                                                xx_6 = xyz{1,fast_obj_start+4}(dist_six_loop,1) - xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1);
                                                                                                                yy_6 = xyz{1,fast_obj_start+4}(dist_six_loop,2) - xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2);
                                                                                                                distance_six_mat(dist_six_loop,1) = sqrt((xx_6^2)+(yy_6^2));                           
                                                                                                            end

                                                                                                            if(dist_final_th <= decide_slow)
                                                                                                                [row_six,col_six,check_six]=find(distance_six_mat < distance_threshold_fast);
                                                                                                            else
                                                                                                                [row_six,col_six,check_six]=find(distance_six_mat < dist_final_th*threshold_diff_dislocation & distance_six_mat <= distance_threshold_fast);
                                                                                                            end

                                                                                                            if (check_six == 1)
                                                                                                                row_six_size=size(row_six);
                                                                                                                angle_six_th=0;
                                                                                                                dist_six_th=0;                         
                                                                                                                priority_angle_order6=0;
                                                                                                                for find_six_loop = 1: row_six_size(1,1)
                                                                                                                    angle_six_th_test=atan((xyz{1,fast_obj_start+4}(row_six(find_six_loop),2)-xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2))/(xyz{1,fast_obj_start+4}(row_six(find_six_loop),1)-xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2)))*180/pi;                                                           
                                                                                                                    priority_angle_order6(find_six_loop,1)=abs(angle_final_th-angle_six_th_test)/priority_degree_rate + distance_six_mat(row_six(find_six_loop));
                                                                                                                    priority_angle_order6(find_six_loop,2)=find_six_loop;%give priority for smaller angle!!
                                                                                                                    if(abs(angle_final_th)>=tangent_control)
                                                                                                                        priority_angle_order6(find_six_loop,1)=abs(abs(angle_final_th)-abs(angle_six_th_test))/priority_degree_rate + distance_six_mat(row_six(find_six_loop));
                                                                                                                    end
                                                                                                                    [temp_pr1,ord_pr1] = sort(priority_angle_order6(:,1));
                                                                                                                    priority_angle_order6 = priority_angle_order6(ord_pr1,:);
                                                                                                                end
                                                                                                                for find_six_loop2 = 1: row_six_size(1,1)
                                                                                                                    pass_numfil=0;  
                                                                                                                    if(stop_search == 0)
                                                                                                                        diff_vector2=0;
                                                                                                                        posit_nega2=nan;
                                                                                                                        if(angle_final_th <= 45 & angle_final_th >= -45) %vetorc
                                                                                                                            vector=1; % x vetorc
                                                                                                                            diff_vector = xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1);
                                                                                                                            posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                        else 
                                                                                                                            vector=0; % y vetorc
                                                                                                                            diff_vector = xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2);
                                                                                                                            posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                        end
                                                                                                                        if(vector == 1)       
                                                                                                                            diff_vector2= xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),1) - xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1);
                                                                                                                            posit_nega2=sign(diff_vector2);                                           
                                                                                                                        elseif(vector == 0)                                            
                                                                                                                            diff_vector2= xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),2) - xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2);
                                                                                                                            posit_nega2=sign(diff_vector2);
                                                                                                                        end
                                                                                                                        if(posit_nega == posit_nega2)
                                                                                                                            pass_value=1;
                                                                                                                            increase_check=1; 
                                                                                                                            if((distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)))+dist_final_th)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),1)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1))^2+(xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),2)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2))^2 ))
                                                                                                                                increase_check=0;                                                     
                                                                                                                            else            
                                                                                                                                increase_check=1;
                                                                                                                            end  
                                                                                                                            if( increase_check == 0)
                                                                                                                                angle_six_th_test=atan((xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),2)-xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2))/(xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),1)-xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1)))*180/pi;                                                                 
                                                                                                                                if( abs(angle_secd_fast) > tangent_control)
                                                                                                                                    if( abs(abs(angle_rest_th)-abs(angle_final_th)) + control_vecter_angle > abs(abs(angle_final_th)-abs(angle_six_th_test)) & abs(abs(angle_rest_th)-abs(angle_final_th)) - control_vecter_angle < abs(abs(angle_final_th)-abs(angle_six_th_test)) )
                                                                                                                                        pass_value=0;
                                                                                                                                    else
                                                                                                                                        pass_value=1;  
                                                                                                                                    end
                                                                                                                                else  
                                                                                                                                    if( abs(angle_rest_th-angle_final_th) + control_vecter_angle > abs(angle_final_th-angle_six_th_test) & abs(angle_rest_th-angle_final_th) + control_vecter_angle > abs(angle_final_th-angle_six_th_test) )
                                                                                                                                        pass_value=0;
                                                                                                                                    else
                                                                                                                                        pass_value=1;  
                                                                                                                                    end
                                                                                                                                end
                                                                                                                            end
                                                                                                                            if(pass_value == 0)              %%%%%%%%%%%%%%%%%%%%%%%%%%                                        
                                                                                                                                if( dist_final_th > decide_slow_angle)
                                                                                                                                    if(abs(angle_final_th)>tangent_control | abs(angle_six_th_test) >tangent_control)    
                                                                                                                                        if(abs(angle_final_th)+control_vecter_angle > abs(angle_six_th_test) & abs(angle_final_th)-control_vecter_angle < abs(angle_six_th_test))  
                                                                                                                                            pass_numfil=1;
                                                                                                                                        end
                                                                                                                                    else
                                                                                                                                        if(angle_final_th+control_vecter_angle > angle_six_th_test & angle_final_th-control_vecter_angle < angle_six_th_test)  
                                                                                                                                            pass_numfil=1;
                                                                                                                                        end
                                                                                                                                    end
                                                                                                                                else
                                                                                                                                    if(abs(angle_final_th)+control_vecter_angle*multilple_short_angle > abs(angle_six_th_test) & abs(angle_final_th)-control_vecter_angle*multilple_short_angle < abs(angle_six_th_test))  
                                                                                                                                        pass_numfil=1;
                                                                                                                                    end  
                                                                                                                                end

                                                                                                                                if( pass_numfil == 1)
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start+4,1) = xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),1);
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start+4,2) = xyz{1,fast_obj_start+4}(row_six(priority_angle_order6(find_six_loop2,2)),2);
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start+4,3) = inside_loop;
                                                                                                                                end
                                                                                                                            end
                                                                                                                        end
                                                                                                                    end
                                                                                                                end
                                                                                                            end
                                                                                                            stop_search=1; %prevent duplicate trajectories
                                                                                                            cont_rest_th(inside_loop,1)=dist_thresh_third;
                                                                                                            cont_angle_th(inside_loop,1)=angle_secd_fast;                   
                                                                                                            if(abs(angle_secd_fast) > tangent_control)
                                                                                                                cont_angle_cmp(inside_loop,1)=abs(angle_fast)-abs(angle_secd_fast);
                                                                                                            else
                                                                                                                cont_angle_cmp(inside_loop,1)=angle_fast-angle_secd_fast;
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
end

%------------------------slow--------------------------------------------------
prevent_duplication_strt2=0;
maxsize_fst_coor= size(fast_total_coor);

for delete_strt_loop= 1:maxsize_fst_coor(1,2)
prevent_duplication_strt2(delete_strt_loop,1)=fast_total_coor{1,delete_strt_loop}(1,1);
prevent_duplication_strt2(delete_strt_loop,2)=fast_total_coor{1,delete_strt_loop}(1,2);
end

xyz_11=0;
comp_prevent_strt=0;
rsh_number=0;
    
xyz_11=xyz{1,1};
comp_prevent_strt=prevent_duplication_strt2;

size_cmp_prevent_1=size(comp_prevent_strt);
    
for loop_cmp_del=1:size_cmp_prevent_1(1,1)
    
    [row_del_x]=find(xyz_11(:,1) == comp_prevent_strt(loop_cmp_del,1));
    [row_del_y]=find(xyz_11(:,2) == comp_prevent_strt(loop_cmp_del,2));
    
    rr_ss=size(row_del_x);
    for check_loop_dlt=1:rr_ss
        [row_del_f,col_del_f,check_del_f]=find(row_del_y == row_del_x(check_loop_dlt));
            
        if(check_del_f ~= 0)
            xyz_11(row_del_x(check_loop_dlt),:)=nan;
        end
    end
end

%****************************************skip frame***************************


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
        start_number_skip = object_number_count+1;
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
                if( estm_ft_size(1,1) > rest_start_searching )
                    
                    xx_1 = xyz{1,fast_obj_start+2}(fast_second_mt,1) - fast_total_coor{1,inside_loop}(rest_start_searching,1); % comparing x coordinates   
                    yy_1 = xyz{1,fast_obj_start+2}(fast_second_mt,2) - fast_total_coor{1,inside_loop}(rest_start_searching,2); % comparing y coordinates
                    
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
                    angle_fast_test = atan((xyz{1,fast_obj_start+2}(row_dist(try_loop),2)-fast_total_coor{1,inside_loop}(rest_start_searching,2))/(xyz{1,fast_obj_start+2}(row_dist(try_loop),1)-fast_total_coor{1,inside_loop}(rest_start_searching,1)))*180/pi; % unit is degree   
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
                        diff_vector = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1) - fast_total_coor{1,inside_loop}(rest_start_searching,1);
                        posit_nega=sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                    else 
                        vector=0; % y vetorc
                        diff_vector = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2) - fast_total_coor{1,inside_loop}(rest_start_searching,2); 
                        posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                    end

                    diff_vector2=0;
                    posit_nega2=nan;  

                    if(vector == 1)         
                        diff_vector2= xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(rest_start_searching,1);    
                        posit_nega2=sign(diff_vector2);
                    elseif (vector == 0)            
                        diff_vector2= xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(rest_start_searching,2);
                        posit_nega2=sign(diff_vector2);          
                    end
                    
                    angle_fast_test=atan((xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(rest_start_searching,2))/(xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(rest_start_searching,1)))*180/pi;

                    if(abs(cont_angle_th(inside_loop,1))+control_vecter_angle*multilple_short_angle > abs(angle_fast_test) & abs(cont_angle_th(inside_loop,1))-control_vecter_angle*multilple_short_angle < abs(angle_fast_test)) % unit is degree
                        pass_numfil=1;
                    end                   
                end
                
                if(posit_nega == posit_nega2) 
                    pass_value=1;   
                    if(fast_obj_start == 1)                                
                        pass_value=0;   
                    else 
                        if((distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)))+cont_rest_th(inside_loop,1))*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(rest_start_searching-1,1))^2+(xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(rest_start_searching-1,2))^2 ))    
                            pass_value=0; 
                        else 
                            pass_value=1; 
                        end 
                    end
                    if(pass_value == 0) 
                        if(pass_numfil == 1) 
                            if(fast_obj_start > 1)      
                                angle_fast = atan((xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(rest_start_searching,2))/(xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(rest_start_searching,1)))*180/pi;   
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
                            [row_secd,col_secd,check_secd]=find(distance_repeat_rest <= distance_threshold_skip);
                                
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
                                                if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-fast_total_coor{1,inside_loop}(rest_start_searching,1))^2+(xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-fast_total_coor{1,inside_loop}(rest_start_searching,2))^2 ))
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
  
                                                        %------prevent duplication--------                   
                                                        maxsize_fstc=size(fast_total_coor);                                        
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
                                                            new_number_count = new_number_count-1; %duplication happen                         
                                                        else                           
                                                            fast_total_coor(1,new_number_count)={fst_mv_oj_coor};   
                                                            cont_rest_th(new_number_count,1)=dist_thresh_sec;                        
                                                            cont_angle_th(new_number_count,1)=angle_fast;                    
                                                            stop_search=1; %prevent duplicate trajectories                              
                                                        end
                                                        
                                                    else
                                                        
                                                        fast_total_coor{1,inside_loop}(rest_start_searching+1,1) = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1);
                                                        fast_total_coor{1,inside_loop}(rest_start_searching+1,2) = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2);                                                      
                                                        fast_total_coor{1,inside_loop}(rest_start_searching+1,3) = inside_loop;

                                                        fast_total_coor{1,inside_loop}(rest_start_searching+2,1) = xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),1);
                                                        fast_total_coor{1,inside_loop}(rest_start_searching+2,2) = xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),2);                                 
                                                        fast_total_coor{1,inside_loop}(rest_start_searching+2,3) = inside_loop;

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

%*******figure*****

figure
imshow(abs(recon_stack(:,:,1)),[])
hold on

dddeee= size(fast_total_coor);

for pictureN=1:dddeee(1,2)
hold on
plot((fast_total_coor{1,pictureN}(:,1)),(fast_total_coor{1,pictureN}(:,2)),'g')
hold on
plot((fast_total_coor{1,pictureN}(1,1)),(fast_total_coor{1,pictureN}(1,2)),'b.')
end

figure
imshow(abs(recon_stack(:,:,1)),[])
hold on
plot(xyz{1,1}(:,1),xyz{1,1}(:,2),'g.')
set(gca, 'XLim', [800, 900], 'YLim', [100,200])
axis on

figure
imshow(abs(recon_stack(:,:,1)),[])
hold on
plot(sum3_xyz(:,1),sum3_xyz(:,2),'g.')
set(gca, 'XLim', [650, 680], 'YLim', [20,70])
axis on

figure, imshow(recon_stack_diff(:,:,1),[])
hold on
plot(xyz{1,1}(:,1),xyz{1,1}(:,2),'g.')
set(gca, 'XLim', [2420, 2520], 'YLim', [550,650])
axis on

set(gca, 'XLim', [900,960], 'YLim', [0,60])
axis on
