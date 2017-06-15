%New trajectories5 by DUCK-HA HWANG

%------solve diff holo dislocation, change xyz coordinate information.

constant_distance_rate=0.8; % endure 80% difference
threshold_diff_dislocation= 1+constant_distance_rate; %*****check it!! important value (1+1*60%) for combining close coordinates
control_vecter_angle = 30; %angle control 45
multilple_short_angle = 1.5; %find wiggling moving for slow one 1.78
increase_distance_threshold=0.81; %make sure about linear moving for fast
combine_limitation=1; % combine xyz noise distance
tangent_control=70;%tangent absolute value 60 degree
priority_degree_rate=9; %9
compare_angle_time=1;
priority_zero=20; %instead of zero value put this value

fps=4; %fps means frame per second
distance_threshold_fast = (140/fps)/pixelsize; %******threshold object tracking. need to control. maximum ditance to distinguish as a same object
%unit is um. // I limited that the fasted sperm's speed is 140 um/sec.
%fps means frame per second

%distance_threshold_skip=4.5;
distance_threshold_skip=(30/fps)/pixelsize; %<----30 before

decide_slow=(40/fps)/pixelsize;
decide_slow_angle=(60/fps)/pixelsize;%criteria for slow moving
%decide_slow=11; %distancce that ovelap happen on the differ_holo_stack

%-----------combine xyz

for combine_compare_loop=1:iter-2    
    behind_compare=combine_compare_loop+1;
    maxsize_cmb_cp=size(xyz{1,combine_compare_loop}); %read size
    for inft_loop=1:maxsize_cmb_cp(1,1)       
        [row_cbx,col_cbx,check_cbx]=find( xyz{1,behind_compare}(:,1) < xyz{1,combine_compare_loop}(inft_loop,1) + combine_limitation & xyz{1,behind_compare}(:,1) > xyz{1,combine_compare_loop}(inft_loop,1) - combine_limitation);       
        if(check_cbx == 1)            
            maxsize_cmb_cp2=size(row_cbx); %read size        
            for bhid_loop=1:maxsize_cmb_cp2(1,1)
                xx_shrt = xyz{1,combine_compare_loop}(inft_loop,1) - xyz{1,behind_compare}(row_cbx(bhid_loop),1); % comparing x coordinates
                yy_shrt = xyz{1,combine_compare_loop}(inft_loop,2) - xyz{1,behind_compare}(row_cbx(bhid_loop),2); % comparing y coordinates

                distance_shrt_limit= sqrt((xx_shrt^2)+(yy_shrt^2));
                if (distance_shrt_limit <= combine_limitation)         
                    xyz{1,behind_compare}(row_cbx(bhid_loop),1) = xyz{1,combine_compare_loop}(inft_loop,1);
                    xyz{1,behind_compare}(row_cbx(bhid_loop),2) = xyz{1,combine_compare_loop}(inft_loop,2);
                end
            end      
        end
    end
end

%------change xyz

fast_obj_start=1; %make a loop

%----------make a combined matrix

sum3_xyz=xyz{1,fast_obj_start};

for sum_fst_loop=fast_obj_start+1:iter-1

    comp_xyz_new=0;
    comp_xyz_new2=0;
    
    comp_xyz_new=xyz{1,sum_fst_loop-1};
    comp_xyz_new2=xyz{1,sum_fst_loop};

    size_cmp_n=size(comp_xyz_new);
    
    for loop_cmp_del=1:size_cmp_n(1,1)
    
        [row_del_x]=find(comp_xyz_new2(:,1) == comp_xyz_new(loop_cmp_del,1));
        [row_del_y]=find(comp_xyz_new2(:,2) == comp_xyz_new(loop_cmp_del,2));
    
        rr_ss=size(row_del_x);
        for check_loop_dlt=1:rr_ss
            [row_del_f,col_del_f,check_del_f]=find(row_del_y == row_del_x(check_loop_dlt));
            
            if(check_del_f ~= 0)
                comp_xyz_new2(row_del_x(check_loop_dlt),:)=nan;
            end
        end
        
    end
    
    sum3_xyz=[sum3_xyz; comp_xyz_new2];
end

%--------------------------------make starting point
xyz_1=0;
xyz_1=xyz;
xyz_change_loop=1;
    
comp_xyz_new=0;
comp_xyz_new2=0;
rsh_number=0;
researching_special=0;

comp_xyz_new=xyz{1,xyz_change_loop};
comp_xyz_new2=xyz{1,xyz_change_loop+1};

size_cmp_n=size(comp_xyz_new2);

for loop_cmp_del=1:size_cmp_n(1,1)

    [row_del_x]=find(comp_xyz_new(:,1) == comp_xyz_new2(loop_cmp_del,1));
    [row_del_y]=find(comp_xyz_new(:,2) == comp_xyz_new2(loop_cmp_del,2));

    rr_ss=size(row_del_x);
    for check_loop_dlt=1:rr_ss
        [row_del_f,col_del_f,check_del_f]=find(row_del_y == row_del_x(check_loop_dlt));

        if(check_del_f ~= 0)
            rsh_number=rsh_number+1;
            researching_special(rsh_number,1)=xyz{1,xyz_change_loop}(row_del_x(check_loop_dlt),1);
            researching_special(rsh_number,2)=xyz{1,xyz_change_loop}(row_del_x(check_loop_dlt),2);
            researching_special(rsh_number,3)=row_del_x(check_loop_dlt);
        end
    end
end

maxsize_rsch_xyz1=size(researching_special); %read size

for find_start_loop=1:maxsize_rsch_xyz1(1,1)-1
    [row_search_strt,col_search_strt,check_search_strt]=find( researching_special(:,1) < researching_special(find_start_loop,1) + distance_threshold_fast & researching_special(:,1) > researching_special(find_start_loop,1) - distance_threshold_fast); 
    maxsize_check_search_strt=size(check_search_strt);
    if(maxsize_check_search_strt(1,1) > 1)
        counting_strt=0;
        distance_strt=0;
        for row_search_strt_loop=1:maxsize_check_search_strt(1,1)
            srt_xx=researching_special((row_search_strt(row_search_strt_loop)),1)-researching_special(find_start_loop,1);
            srt_yy=researching_special((row_search_strt(row_search_strt_loop)),2)-researching_special(find_start_loop,2);
            if(srt_xx == 0)
                srt_xx=nan;
            end
            if(srt_yy == 0)
                srt_yy=nan;
            end
            counting_strt=counting_strt+1;
            distance_strt(counting_strt,1)=sqrt(srt_xx^2+srt_yy^2);
            distance_strt(counting_strt,2)=row_search_strt_loop;
        end
        [row_distrt,col_distrt,check_distrt]= find(distance_strt(:,1) <= distance_threshold_fast);
        if(check_distrt == 1)           
            left_only=nan;
            left_one_moving=0;
            left_one_moving(1,1)=min_intensities_xyz{1,xyz_change_loop}(researching_special(find_start_loop,3));
            left_one_moving(1,2)=researching_special(find_start_loop,3);
            maxsize_row_distrt=size(row_distrt);
            for row_distrt_loop=1:maxsize_row_distrt(1,1)
                left_one_moving(row_distrt_loop+1,1)=min_intensities_xyz{1,xyz_change_loop}(researching_special(row_search_strt(distance_strt(row_distrt(row_distrt_loop),2)),3));
                left_one_moving(row_distrt_loop+1,2)=researching_special(row_search_strt(distance_strt(row_distrt(row_distrt_loop),2)),3);
            end
            left_only=max(left_one_moving(:,1));
            maxsize_left_one_moving=size(left_one_moving);
            for left_one_moving_loop=1:maxsize_left_one_moving(1,1)
                if(left_one_moving(left_one_moving_loop,1) ~= left_only)
                    xyz_1{1,xyz_change_loop}(left_one_moving(left_one_moving_loop,2),:) = nan;
                end
            end
        else
            [row_search_strt2,col_search_strt2,check_search_strt2]=find( xyz{1,xyz_change_loop}(:,1) < xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1) + distance_threshold_fast & xyz{1,xyz_change_loop}(:,1) > xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1) - distance_threshold_fast);
            maxsize_check_search_strt2=size(check_search_strt2);
            if(maxsize_check_search_strt2(1,1) > 1)
                counting_strt=0;
                distance_strt=0;
                check_distrt=0;
                row_distrt=0;
                col_distrt=0;
                for row_search_strt_loop=1:maxsize_check_search_strt2(1,1)
                    srt_xx=xyz{1,xyz_change_loop}(row_search_strt2(row_search_strt_loop),1)-xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1);
                    srt_yy=xyz{1,xyz_change_loop}(row_search_strt2(row_search_strt_loop),2)-xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),2);
                    if(srt_xx == 0)
                        srt_xx=nan;
                    end
                    if(srt_yy == 0)
                        srt_yy=nan;
                    end
                    counting_strt=counting_strt+1;
                    distance_strt(counting_strt,1)=sqrt(srt_xx^2+srt_yy^2);
                    distance_strt(counting_strt,2)=row_search_strt_loop;
                end
                [row_distrt,col_distrt,check_distrt]= find(distance_strt(:,1) <= distance_threshold_fast);
                if(check_distrt == 1)
                    xyz_1{1,xyz_change_loop}(researching_special(find_start_loop,3),:)=nan;
                end
            end
        end
    else
        [row_search_strt3,col_search_strt3,check_search_strt3]=find( xyz{1,xyz_change_loop}(:,1) < xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1) + distance_threshold_fast & xyz{1,xyz_change_loop}(:,1) > xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1) - distance_threshold_fast);
        maxsize_check_search_strt3=size(check_search_strt3);
        if(maxsize_check_search_strt3(1,1) > 1)
            counting_strt=0;
            distance_strt=0;
            check_distrt=0;
            row_distrt=0;
            col_distrt=0;
            for row_search_strt_loop=1:maxsize_check_search_strt3(1,1)
                srt_xx=xyz{1,xyz_change_loop}(row_search_strt3(row_search_strt_loop),1)-xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1);
                srt_yy=xyz{1,xyz_change_loop}(row_search_strt3(row_search_strt_loop),2)-xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),2);
                if(srt_xx == 0)
                    srt_xx=nan;
                end
                if(srt_yy == 0)
                    srt_yy=nan;
                end
                counting_strt=counting_strt+1;
                distance_strt(counting_strt,1)=sqrt(srt_xx^2+srt_yy^2);
                distance_strt(counting_strt,2)=row_search_strt_loop;
            end
            [row_distrt,col_distrt,check_distrt]= find(distance_strt(:,1) <= distance_threshold_fast);
            if(check_distrt == 1)
                xyz_1{1,xyz_change_loop}(researching_special(find_start_loop,3),:)=nan;
            end
        end
    end
end

%****************************************main trajectory***************************
max_frame_size=fix((iter-2)/4);

frame_number=4*max_frame_size+1

object_number_count=0;
fast_total_coor={0};
angle_information=0;

for change_start_point= 0:2

    for rest_start_searching= 1:max_frame_size

        prevent_rest_check={0};
        fast_obj_start=(4*rest_start_searching-3)+change_start_point;
        if(fast_obj_start == 1+change_start_point)
            maxsize_fast_xyz=size(xyz{1,fast_obj_start}); %read size
            change_inside_loop = maxsize_fast_xyz(1,1);
        else
            change_inside_loop = object_number_count;
        end
        duple_num_count=0;
        for inside_loop=1:change_inside_loop
            stop_search=0;
            distance_fast_mat=0;        
            if(fast_obj_start == 1+change_start_point)
                [row_f1,col_f1,check_f1]=find( xyz{1,fast_obj_start+1}(:,1) < xyz_1{1,fast_obj_start}(inside_loop,1) + distance_threshold_fast & xyz{1,fast_obj_start+1}(:,1) > xyz_1{1,fast_obj_start}(inside_loop,1) - distance_threshold_fast);   
                xc_1=xyz_1{1,fast_obj_start}(inside_loop,1);
                yc_1=xyz_1{1,fast_obj_start}(inside_loop,2);
            else
                estm_ft_size=size(fast_total_coor{1,inside_loop});
                if( estm_ft_size(1,1) >= fast_obj_start-change_start_point )
                    [row_f1,col_f1,check_f1]=find( xyz{1,fast_obj_start+1}(:,1) < fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point ,1) + distance_threshold_fast & xyz{1,fast_obj_start+1}(:,1) > fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,1) - distance_threshold_fast);
                    xc_1=fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,1);
                    yc_1=fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,2);
                    xc_back_1=fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point-1,1);
                    yc_back_1=fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point-1,2);
                    xc_back_2=fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point-2,1);
                    yc_back_2=fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point-2,2);
                end
            end
            if(check_f1 == 1)  
                maxsize_fast2_xyz=size(row_f1); %read size
                distance_count1=0;
                for fast_second_mt=1:maxsize_fast2_xyz(1,1) %for fast sperm
                    if(fast_obj_start == 1+change_start_point)
                        xx_1 = xyz{1,fast_obj_start+1}(row_f1(fast_second_mt),1) - xc_1; % comparing x coordinates   
                        yy_1 = xyz{1,fast_obj_start+1}(row_f1(fast_second_mt),2) - yc_1; % comparing y coordinates
                        distance_count1=distance_count1+1;
                        distance_fast_mat(distance_count1,1) = sqrt((xx_1^2)+(yy_1^2));
                        distance_fast_mat(distance_count1,2) = row_f1(fast_second_mt);
                    else
                        estm_ft_size=size(fast_total_coor{1,inside_loop});
                        if( estm_ft_size(1,1) >= fast_obj_start-change_start_point )
                            xx_1 = xyz{1,fast_obj_start+1}(row_f1(fast_second_mt),1) - xc_1; % comparing x coordinates   
                            yy_1 = xyz{1,fast_obj_start+1}(row_f1(fast_second_mt),2) - yc_1; % comparing y coordinates
                        else     
                            xx_1=nan;
                            yy_1=nan;
                        end
                        distance_count1=distance_count1+1;
                        distance_fast_mat(distance_count1,1) = sqrt((xx_1^2)+(yy_1^2));
                        distance_fast_mat(distance_count1,2) = row_f1(fast_second_mt);
                    end
                end
                if(fast_obj_start == 1+change_start_point) 
                    [row_dist,col_dist,check_dist]=find(distance_fast_mat(:,1) <= distance_threshold_fast);     
                else
                    before_distance=sqrt((xc_1-xc_back_1)^2+(yc_1-yc_back_1)^2);
                    before_angle=atan((yc_1-yc_back_1)/(xc_1-xc_back_1))*180/pi;
                    if( (yc_1-yc_back_1) == 0 & (xc_1-xc_back_1) == 0)
                        before_angle = 0;
                    end 
                    if( before_distance <= decide_slow)
                        [row_dist,col_dist,check_dist]=find(distance_fast_mat(:,1) <= distance_threshold_fast);            
                    else
                        [row_dist,col_dist,check_dist]=find(distance_fast_mat(:,1) <= before_distance*threshold_diff_dislocation & distance_fast_mat(:,1) <= distance_threshold_fast);
                    end
                end
                if (check_dist == 1) %first distance threshold check
                    row_dist_size=size(row_dist); %read size  
                    angle_fast=0;
                    dist_thresh_sec=0;
                    priority_angle_order=0;
                    if(fast_obj_start == 1+change_start_point)  
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
                            angle_fast_test = atan((xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(try_loop),2),2)-yc_back_1)/(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(try_loop),2),1)-xc_back_1))*180/pi; % unit is degree   
                            priority_angle_order(try_loop,1)=abs(abs(angle_fast_test)-abs(before_angle))/priority_degree_rate + distance_fast_mat(row_dist(try_loop));                  
                            if( (xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(try_loop),2),2)-yc_back_1) == 0 & (xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(try_loop),2),1)-xc_back_1) == 0 )
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
                            xc_2=xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1);
                            yc_2=xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2);
                            if(fast_obj_start == 1+change_start_point)
                                %decide vector andgle
                                %each magnitude
                                angle_fast = atan((yc_2-yc_1)/(xc_2-xc_1))*180/pi; % unit is degree
                                if( (yc_2-yc_1) == 0 & (xc_2-xc_1) == 0 )
                                    angle_fast = 0;
                                end
                                dist_thresh_sec = distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)));
                                if(angle_fast <= 45 & angle_fast >= -45) %vetorc
                                    vector=1; % x vetorc
                                    diff_vector = xc_2 - xc_1;
                                    posit_nega=sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                else
                                    vector=0; % y vetorc
                                    diff_vector = yc_2 - yc_1; 
                                    posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                end
                                posit_nega2 = posit_nega;
                                pass_numfil=1;
                            else
                                if(before_angle <= 45 & before_angle >= -45) %vetorc
                                    vector=1; % x vetorc
                                    diff_vector = xc_1 - xc_back_1 ;
                                    posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                else 
                                    vector=0; % y vetorc
                                    diff_vector = yc_1 - yc_back_1 ;
                                    posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                end
                                diff_vector2=0;
                                posit_nega2=nan;
                                if(vector == 1)         
                                    diff_vector2= xc_2-xc_1;    
                                    posit_nega2=sign(diff_vector2);
                                elseif (vector == 0)            
                                    diff_vector2= yc_2-yc_1;
                                    posit_nega2=sign(diff_vector2);          
                                end                      
                                if( posit_nega == 0 )                     
                                    posit_nega = posit_nega2;
                                    xc_back_1=xc_back_2;
                                    yc_back_1=yc_back_2;
                                end
                                if(posit_nega == posit_nega2)                            
                                    before_angle=atan((yc_1-yc_back_1)/(xc_1-xc_back_1))*180/pi;
                                    if( (yc_1-yc_back_1) == 0 & (xc_1-xc_back_1) == 0)
                                        before_angle = 0;
                                    end        
                                    before_distance=sqrt((xc_1-xc_back_1)^2+(yc_1-yc_back_1)^2);
                                    angle_fast_test=atan((yc_2-yc_1)/(xc_2-xc_1))*180/pi;
                                    if( (yc_2-yc_1) == 0 & (xc_2-xc_1) == 0)
                                        angle_fast_test = 0;
                                    end
                                    if( before_distance > decide_slow_angle)           
                                        if(abs(before_angle) > tangent_control | abs(angle_fast_test) > tangent_control)              
                                            if(abs(before_angle)+control_vecter_angle > abs(atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi) & abs(before_angle)-control_vecter_angle < abs(atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi)) % unit is degree
                                                pass_numfil=1;
                                            end
                                        else
                                            if(before_angle+control_vecter_angle > atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi & before_angle-control_vecter_angle < atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi) % unit is degree
                                            pass_numfil=1;
                                            end
                                        end    
                                    else
                                        if(abs(before_angle) > tangent_control | abs(angle_fast_test) > tangent_control)              
                                            if(abs(before_angle)+control_vecter_angle*multilple_short_angle > abs(atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi) & abs(before_angle)-control_vecter_angle*multilple_short_angle < abs(atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi)) % unit is degree
                                                pass_numfil=1;
                                            end
                                        else
                                            if(before_angle+control_vecter_angle*multilple_short_angle > atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi & before_angle-control_vecter_angle*multilple_short_angle < atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi) % unit is degree
                                            pass_numfil=1;
                                            end
                                        end
                                    end                                   
                                end
                            end
                            increase_check=1;
                            if(fast_obj_start == 1+change_start_point)                         
                                increase_check=0; 
                            else
                                if((distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)))+before_distance)*increase_distance_threshold < sqrt( (xc_2-xc_back_1)^2+(yc_2-yc_back_1)^2 ))
                                    increase_check=0; 
                                else   
                                    increase_check=1;
                                end        
                            end
                            if(posit_nega2 == 0)
                                pass_numfil=1;
                                increase_check=0;
                            end
                            if(increase_check == 0)       
                                if(pass_numfil == 1) 
                                    if(fast_obj_start > 1+change_start_point)
                                        angle_fast = atan((yc_2-yc_1)/(xc_2-xc_1))*180/pi;
                                        if( (yc_2-yc_1) == 0 & (xc_2-xc_1) == 0)
                                            angle_fast = 0;
                                        end
                                        dist_thresh_sec = distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)));                                
                                    end   
                                    distance_repeat_rest=0;
                                    [row_f2,col_f2,check_f2]=find( xyz{1,fast_obj_start+2}(:,1) < xc_2 + distance_threshold_fast & xyz{1,fast_obj_start+2}(:,1) > xc_2 - distance_threshold_fast);       

                                    if(check_f2 == 1)
                                        maxsize_fast_xyz2=size(row_f2); %read size
                                        distance_count2=0;
                                        for dist_cmp_loop = 1: maxsize_fast_xyz2(1,1)
                                            xx_2 = xyz{1,fast_obj_start+2}(row_f2(dist_cmp_loop),1) - xc_2;                              
                                            yy_2 = xyz{1,fast_obj_start+2}(row_f2(dist_cmp_loop),2) - yc_2;                 
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
                                                angle_secd_fast_test = atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),2)-yc_1)/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),1)-xc_1))*180/pi;
                                                priority_angle_order2(find_secd_loop,1)=abs(abs(angle_secd_fast_test)-abs(angle_fast))/priority_degree_rate + distance_repeat_rest(row_secd(find_secd_loop));
                                                if( (xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),2)-yc_1) == 0 & (xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),1)-xc_1) == 0 )
                                                    priority_angle_order2(find_secd_loop,1) = priority_zero; 
                                                end
                                                priority_angle_order2(find_secd_loop,2)=find_secd_loop;%give priority for smaller angle!!
                                                [temp_pr1,ord_pr1] = sort(priority_angle_order2(:,1));
                                                priority_angle_order2 = priority_angle_order2(ord_pr1,:);
                                            end
                                            for find_secd_loop2 = 1: row_secd_size(1,1)
                                                if(stop_search == 0)
                                                    xc_3=xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1);
                                                    yc_3=xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2);
                                                    pass_numfil=0;
                                                    diff_vector2=0;
                                                    posit_nega2=nan;
                                                    if( fast_obj_start > 1+change_start_point)
                                                         if(angle_fast <= 45 & angle_fast >= -45) %vetorc
                                                            vector=1; % x vetorc
                                                            diff_vector = xc_2-xc_1;
                                                            posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                        else 
                                                            vector=0; % y vetorc
                                                            diff_vector = yc_2-yc_1;
                                                            posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                         end
                                                    end
                                                    if(vector == 1)
                                                        diff_vector2= xc_3-xc_2;
                                                        posit_nega2=sign(diff_vector2);
                                                    elseif (vector == 0)
                                                        diff_vector2= yc_3-yc_2;                            
                                                        posit_nega2=sign(diff_vector2);
                                                    end
                                                    if( posit_nega == 0 )
                                                        if( fast_obj_start > 1+change_start_point)
                                                            xc_1=xc_back_1;
                                                            yc_1=yc_back_1;
                                                        end
                                                        posit_nega = posit_nega2;        
                                                    end
                                                    if(posit_nega == posit_nega2)
                                                        angle_fast=atan((yc_2-yc_1)/(xc_2-xc_1))*180/pi;
                                                        if( (yc_2-yc_1) == 0 & (xc_2-xc_1) == 0 )
                                                            angle_fast=0;
                                                        end
                                                        dist_thresh_sec=sqrt((xc_2-xc_1)^2+(yc_2-yc_1)^2);
                                                        increase_check=1;
                                                        if(fast_obj_start == 1+change_start_point)
                                                            if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xc_3-xc_1)^2+(yc_3-yc_1)^2 ))
                                                                increase_check=0;
                                                            else
                                                                increase_check=1;
                                                            end
                                                        else
                                                            if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xc_3-xc_1)^2+(yc_3-yc_1)^2 ))
                                                                increase_check=0;
                                                            else
                                                                increase_check=1;
                                                            end
                                                        end
                                                        if(posit_nega2 == 0)
                                                            pass_numfil=1;
                                                            increase_check=0;
                                                        end
                                                        if(increase_check == 0)
                                                            angle_secd_fast_cmp=atan((yc_3-yc_2)/(xc_3-xc_2))*180/pi; 
                                                            if( (yc_3-yc_2) == 0 & (xc_3-xc_2) == 0 )
                                                                angle_secd_fast_cmp=0;
                                                            end
                                                            if(fast_obj_start==1+change_start_point)
                                                                if(dist_thresh_sec > decide_slow_angle)
                                                                    if(abs(angle_fast)>tangent_control | abs(angle_secd_fast_cmp) > tangent_control)                                                     
                                                                        if(abs(angle_fast)+control_vecter_angle > abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) & abs(angle_fast)-control_vecter_angle < abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi)) % unit is degree                       
                                                                            pass_numfil=1;
                                                                        end
                                                                    else
                                                                        if(angle_fast+control_vecter_angle > atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi & angle_fast-control_vecter_angle < atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) % unit is degree  
                                                                            pass_numfil=1;
                                                                        end
                                                                    end
                                                                else
                                                                    if(abs(angle_fast)>tangent_control | abs(angle_secd_fast_cmp) > tangent_control)                                                     
                                                                        if(abs(angle_fast)+control_vecter_angle*multilple_short_angle > abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) & abs(angle_fast)-control_vecter_angle*multilple_short_angle < abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi)) % unit is degree                       
                                                                            pass_numfil=1;
                                                                        end
                                                                    else
                                                                        if(angle_fast+control_vecter_angle*multilple_short_angle > atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi & angle_fast-control_vecter_angle*multilple_short_angle < atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) % unit is degree  
                                                                            pass_numfil=1;
                                                                        end
                                                                    end
                                                                end
                                                            else
                                                                if(dist_thresh_sec > decide_slow_angle)
                                                                    if(abs(angle_fast)>tangent_control | abs(angle_secd_fast_cmp) > tangent_control)                                                     
                                                                        if(abs(angle_fast)+control_vecter_angle > abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) & abs(angle_fast)-control_vecter_angle < abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi)) % unit is degree                       
                                                                            pass_numfil=1;
                                                                        end
                                                                    else
                                                                        if(angle_fast+control_vecter_angle > atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi & angle_fast-control_vecter_angle < atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) % unit is degree  
                                                                            pass_numfil=1;
                                                                        end
                                                                    end
                                                                else
                                                                    if(abs(angle_fast)>tangent_control | abs(angle_secd_fast_cmp) > tangent_control)                                                     
                                                                        if(abs(angle_fast)+control_vecter_angle*multilple_short_angle > abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) & abs(angle_fast)-control_vecter_angle*multilple_short_angle < abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi)) % unit is degree                       
                                                                            pass_numfil=1;
                                                                        end
                                                                    else
                                                                        if(angle_fast+control_vecter_angle*multilple_short_angle > atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi & angle_fast-control_vecter_angle*multilple_short_angle < atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) % unit is degree  
                                                                            pass_numfil=1;
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                            if(pass_numfil == 1) 
                                                                angle_secd_fast = atan((yc_3-yc_2)/(xc_3-xc_2))*180/pi;   
                                                                if( (yc_3-yc_2) == 0 & (xc_3-xc_2) == 0 )
                                                                    angle_secd_fast=0;
                                                                end
                                                                dist_thresh_third = distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)));                           
                                                                % finish third
                                                                distance_thrd_mat=0;
                                                                [row_f3,col_f3,check_f3]=find( xyz{1,fast_obj_start+3}(:,1) < xc_3 + distance_threshold_fast & xyz{1,fast_obj_start+3}(:,1) > xc_3 - distance_threshold_fast);       

                                                                if(check_f3 == 1)
                                                                    maxsize_fast_xyz3=size(row_f3);
                                                                    distance_count3=0;
                                                                    for dist_thrd_loop = 1: maxsize_fast_xyz3(1,1)                                   
                                                                        xx_3 = xyz{1,fast_obj_start+3}(row_f3(dist_thrd_loop),1) - xc_3;
                                                                        yy_3 = xyz{1,fast_obj_start+3}(row_f3(dist_thrd_loop),2) - yc_3;
                                                                        distance_count3=distance_count3+1;
                                                                        distance_thrd_mat(distance_count3,1) = sqrt((xx_3^2)+(yy_3^2));
                                                                        distance_thrd_mat(distance_count3,2) = row_f3(dist_thrd_loop);
                                                                    end
                                                                    if(dist_thresh_third <= decide_slow)
                                                                        [row_third,col_third,check_third]=find(distance_thrd_mat(:,1) < distance_threshold_fast);
                                                                    else
                                                                        [row_third,col_third,check_third]=find(distance_thrd_mat(:,1) < dist_thresh_third*threshold_diff_dislocation & distance_thrd_mat(:,1) <= distance_threshold_fast );
                                                                    end
                                                                    if (check_third == 1)
                                                                        row_third_size=size(row_third);
                                                                        angle_rest_th=0;
                                                                        dist_rest_th=0;                         
                                                                        priority_angle_order3=0;
                                                                        for find_rest_loop = 1: row_third_size(1,1)
                                                                            angle_rest_th_test=atan((xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),2)-yc_2)/(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),1)-xc_2))*180/pi;                                                           
                                                                            priority_angle_order3(find_rest_loop,1)=abs(abs(angle_rest_th_test)-abs(angle_secd_fast))/priority_degree_rate + distance_thrd_mat(row_third(find_rest_loop));
                                                                            if( (xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),2)-yc_2) == 0 & (xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),1)-xc_2) ==0 )
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
                                                                                xc_4=xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1);
                                                                                yc_4=xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2);
                                                                                if(angle_secd_fast <= 45 & angle_secd_fast >= -45) %vetorc
                                                                                    vector=1; % x vetorc
                                                                                    diff_vector = xc_3-xc_2;
                                                                                    posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                else 
                                                                                    vector=0; % y vetorc
                                                                                    diff_vector = yc_3-yc_2;
                                                                                    posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                end
                                                                                if(vector == 1)       
                                                                                    diff_vector2= xc_4 - xc_3;
                                                                                    posit_nega2=sign(diff_vector2);                                           
                                                                                elseif(vector == 0)                                            
                                                                                    diff_vector2= yc_4 - yc_3;                                            
                                                                                    posit_nega2=sign(diff_vector2);
                                                                                end
                                                                                if( posit_nega == 0 )
                                                                                    xc_2=xc_1;
                                                                                    yc_2=yc_1;
                                                                                    posit_nega = posit_nega2;        
                                                                                end
                                                                                if(posit_nega == posit_nega2)
                                                                                    angle_secd_fast=atan((yc_3-yc_2)/(xc_3-xc_2))*180/pi;
                                                                                    if( (yc_3-yc_2) == 0 & (xc_3-xc_2) == 0 )
                                                                                        angle_secd_fast = 0;
                                                                                    end
                                                                                    dist_thresh_third=sqrt((xc_3-xc_2)^2+(yc_3-yc_2)^2);
                                                                                    increase_check=1;
                                                                                    if((distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)))+dist_thresh_third)*increase_distance_threshold < sqrt( (xc_4-xc_2)^2+(yc_4-yc_2)^2 ))
                                                                                        increase_check=0;                                                     
                                                                                    else            
                                                                                        increase_check=1;
                                                                                    end 
                                                                                    if(posit_nega2 == 0)
                                                                                        pass_numfil=1;
                                                                                        increase_check=0;
                                                                                    end
                                                                                    if(increase_check == 0)  
                                                                                        angle_rest_th_cmp=atan((yc_4-yc_3)/(xc_4-xc_3))*180/pi;
                                                                                        if( (yc_4-yc_3) == 0 & (xc_4-xc_3)==0)
                                                                                            angle_rest_th_cmp=0;
                                                                                        end
                                                                                        if( dist_thresh_third > decide_slow_angle)
                                                                                            if(abs(angle_secd_fast)>tangent_control | abs(angle_rest_th_cmp)>tangent_control)    
                                                                                                if(abs(angle_secd_fast)+control_vecter_angle > abs(atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi) & abs(angle_secd_fast)-control_vecter_angle < abs(atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi))  
                                                                                                    pass_numfil=1;
                                                                                                end
                                                                                            else
                                                                                                if(angle_secd_fast+control_vecter_angle > atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi & angle_secd_fast-control_vecter_angle < atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi)  
                                                                                                    pass_numfil=1;
                                                                                                end
                                                                                            end
                                                                                        else
                                                                                            if(abs(angle_secd_fast)>tangent_control | abs(angle_rest_th_cmp)>tangent_control)    
                                                                                                if(abs(angle_secd_fast)+control_vecter_angle*multilple_short_angle > abs(atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi) & abs(angle_secd_fast)-control_vecter_angle*multilple_short_angle < abs(atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi))  
                                                                                                    pass_numfil=1;
                                                                                                end
                                                                                            else
                                                                                                if(angle_secd_fast+control_vecter_angle*multilple_short_angle > atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi & angle_secd_fast-control_vecter_angle*multilple_short_angle < atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi)  
                                                                                                    pass_numfil=1;
                                                                                                end
                                                                                            end 
                                                                                        end
                                                                                        if( pass_numfil == 1) 
                                                                                            angle_rest_th=atan((yc_4-yc_3)/(xc_4-xc_3))*180/pi;                 
                                                                                            if ( (yc_4-yc_3) == 0 & (xc_4-xc_3) == 0 )
                                                                                                angle_rest_th = 0;
                                                                                            end
                                                                                            dist_rest_th=distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)));
                                                                                            [row_f4,col_f4,check_f4]=find( xyz{1,fast_obj_start+4}(:,1) < xc_4 + distance_threshold_fast & xyz{1,fast_obj_start+4}(:,1) > xc_4 - distance_threshold_fast);       

                                                                                            if(check_f4 == 1)
                                                                                                distance_four_mat=0;
                                                                                                maxsize_fast_xyz4=size(row_f4);
                                                                                                distance_count4=0;
                                                                                                for dist_four_loop = 1: maxsize_fast_xyz4(1,1)
                                                                                                    xx_4 = xyz{1,fast_obj_start+4}(row_f4(dist_four_loop),1) - xc_4;    
                                                                                                    yy_4 = xyz{1,fast_obj_start+4}(row_f4(dist_four_loop),2) - yc_4;                                
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
                                                                                                        angle_final_th_test=atan((xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),2)-yc_3)/(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),1)-xc_3))*180/pi;
                                                                                                        priority_angle_order4(find_four_loop,1)=abs(abs(angle_rest_th)-abs(angle_final_th_test))/priority_degree_rate + distance_four_mat(row_four(find_four_loop));
                                                                                                        if( (xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),2)-yc_3) == 0 & (xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),1)-xc_3) == 0)
                                                                                                            priority_angle_order4(find_four_loop,1)=priority_zero;
                                                                                                        end
                                                                                                        priority_angle_order4(find_four_loop,2)=find_four_loop;%give priority for smaller angle!!
                                                                                                        [temp_pr1,ord_pr1] = sort(priority_angle_order4(:,1));
                                                                                                        priority_angle_order4 = priority_angle_order4(ord_pr1,:);
                                                                                                    end
                                                                                                    for find_four_loop2 = 1: row_four_size(1,1)                                                          
                                                                                                        if(stop_search == 0)
                                                                                                            xc_5=xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1);
                                                                                                            yc_5=xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2);
                                                                                                            pass_numfil=0;
                                                                                                            diff_vector2=0;
                                                                                                            posit_nega2=nan;
                                                                                                            if(angle_rest_th <= 45 & angle_rest_th >= -45) %vetorc
                                                                                                                vector=1; % x vetorc
                                                                                                                diff_vector = xc_4-xc_3;
                                                                                                                posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                            else 
                                                                                                                vector=0; % y vetorc
                                                                                                                diff_vector = yc_4-yc_3;
                                                                                                                posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                            end
                                                                                                            if(vector == 1)       
                                                                                                                diff_vector2= xc_5 - xc_4;      
                                                                                                                posit_nega2=sign(diff_vector2);        
                                                                                                            elseif(vector == 0)                                            
                                                                                                                diff_vector2= yc_5 - yc_4;                                            
                                                                                                                posit_nega2=sign(diff_vector2);
                                                                                                            end
                                                                                                            if( posit_nega == 0 )
                                                                                                                xc_3=xc_2;
                                                                                                                yc_3=yc_2;
                                                                                                                posit_nega = posit_nega2;        
                                                                                                            end
                                                                                                            if(posit_nega == posit_nega2)
                                                                                                                angle_rest_th=atan((yc_4-yc_3)/(xc_4-xc_3))*180/pi;
                                                                                                                if( (yc_4-yc_3) == 0 & (xc_4-xc_3) == 0 )
                                                                                                                    angle_rest_th = 0;
                                                                                                                end
                                                                                                                dist_rest_th= sqrt((xc_4-xc_3)^2+(yc_4-yc_3)^2);
                                                                                                                increase_check=1;
                                                                                                                if((distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)))+dist_rest_th)*increase_distance_threshold < sqrt((xc_5-xc_3)^2+(yc_5-yc_3)^2))       
                                                                                                                    increase_check=0;    
                                                                                                                else
                                                                                                                    increase_check=1;      
                                                                                                                end
                                                                                                                if(posit_nega2 == 0)
                                                                                                                    pass_numfil=1;
                                                                                                                    increase_check=0;
                                                                                                                end
                                                                                                                if(increase_check == 0)                                                                               
                                                                                                                    angle_final_th_cmp=atan((yc_5-yc_4)/(xc_5-xc_4))*180/pi;
                                                                                                                    if( (yc_5-yc_4) == 0 & (xc_5-xc_4) ==0 )
                                                                                                                        angle_final_th_cmp=0;
                                                                                                                    end
                                                                                                                    if( dist_rest_th > decide_slow_angle)
                                                                                                                        if(abs(angle_rest_th)>tangent_control | abs(angle_final_th_cmp) >tangent_control)
                                                                                                                            if(abs(angle_rest_th)+control_vecter_angle > abs(atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi) & abs(angle_rest_th)-control_vecter_angle < abs(atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi))
                                                                                                                                pass_numfil=1;
                                                                                                                            end
                                                                                                                        else
                                                                                                                            if(angle_rest_th+control_vecter_angle > atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi & angle_rest_th-control_vecter_angle < atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi)
                                                                                                                                pass_numfil=1;
                                                                                                                            end
                                                                                                                        end
                                                                                                                    else  
                                                                                                                         if(abs(angle_rest_th)>tangent_control | abs(angle_final_th_cmp) >tangent_control)
                                                                                                                            if(abs(angle_rest_th)+control_vecter_angle*multilple_short_angle > abs(atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi) & abs(angle_rest_th)-control_vecter_angle*multilple_short_angle  < abs(atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi))
                                                                                                                                pass_numfil=1;
                                                                                                                            end
                                                                                                                        else
                                                                                                                            if(angle_rest_th+control_vecter_angle*multilple_short_angle  > atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi & angle_rest_th-control_vecter_angle*multilple_short_angle  < atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi)
                                                                                                                                pass_numfil=1;
                                                                                                                            end
                                                                                                                        end
                                                                                                                    end
                                                                                                                    if( pass_numfil == 1)
                                                                                                                        angle_final_th=atan((yc_5-yc_4)/(xc_5-xc_4))*180/pi;                                    
                                                                                                                        if( (yc_5-yc_4) == 0 & (xc_5-xc_4)==0)
                                                                                                                            angle_final_th=0;
                                                                                                                        end
                                                                                                                        dist_final_th=distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)));
                                                                                                                        if(fast_obj_start == 1+change_start_point) 
                                                                                                                            object_number_count = object_number_count+1;
                                                                                                                            fst_mv_oj_coor=0;

                                                                                                                            fst_mv_oj_coor(1,1) = xyz{1,fast_obj_start}(inside_loop,1);                              
                                                                                                                            fst_mv_oj_coor(1,2) = xyz{1,fast_obj_start}(inside_loop,2);          
                                                                                                                            fst_mv_oj_coor(1,3) = object_number_count;
                                                                                                                            fst_mv_oj_coor(1,4) = inside_loop;

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
                                                                                                                                angle_information(object_number_count,1)=(abs(abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(1,1)))*180/pi)-abs(atan((fst_mv_oj_coor(2,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(2,1)-fst_mv_oj_coor(1,1)))*180/pi))+abs(abs(atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(2,1)))*180/pi)-abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(2,1)))*180/pi))+abs(abs(atan((fst_mv_oj_coor(5,2)-fst_mv_oj_coor(3,2))/(fst_mv_oj_coor(5,1)-fst_mv_oj_coor(3,1)))*180/pi)-abs(atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(3,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(3,1)))*180/pi)))/priority_degree_rate + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) );
                                                                                                                            else
                                                                                                                                angle_information(object_number_count,1)=(abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(1,1)))*180/pi-atan((fst_mv_oj_coor(2,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(2,1)-fst_mv_oj_coor(1,1)))*180/pi)+abs(atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(2,1)))*180/pi-atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(2,1)))*180/pi)+abs(atan((fst_mv_oj_coor(5,2)-fst_mv_oj_coor(3,2))/(fst_mv_oj_coor(5,1)-fst_mv_oj_coor(3,1)))*180/pi-atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(3,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(3,1)))*180/pi))/priority_degree_rate + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) );
                                                                                                                            end
                                                                                                                            distance_six_mat=0;
                                                                                                                            [row_f5,col_f5,check_f5]=find( xyz{1,fast_obj_start+4}(:,1) < xc_5 + distance_threshold_fast & xyz{1,fast_obj_start+4}(:,1) > xc_5 - distance_threshold_fast);       

                                                                                                                            if(check_f5 == 1)
                                                                                                                                stop_search2=0;
                                                                                                                                maxsize_fast_xyz6=size(row_f5);
                                                                                                                                distance_count6=0;
                                                                                                                                for dist_six_loop = 1: maxsize_fast_xyz6(1,1)                                   
                                                                                                                                    xx_6 = xyz{1,fast_obj_start+4}(row_f5(dist_six_loop),1) - xc_5;
                                                                                                                                    yy_6 = xyz{1,fast_obj_start+4}(row_f5(dist_six_loop),2) - yc_5;                                                                                                                                       
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
                                                                                                                                        angle_six_th_test=atan((xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),2)-yc_4)/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),1)-xc_4))*180/pi;                                                           
                                                                                                                                        priority_angle_order6(find_six_loop,1)=abs(abs(angle_final_th)-abs(angle_six_th_test))/priority_degree_rate + distance_six_mat(row_six(find_six_loop));
                                                                                                                                        if( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),2)-yc_4) == 0 & (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),1)- xc_4) )
                                                                                                                                            priority_angle_order6(find_six_loop,1)= priority_zero; 
                                                                                                                                        end
                                                                                                                                        priority_angle_order6(find_six_loop,2)=find_six_loop;%give priority for smaller angle!!
                                                                                                                                        [temp_pr1,ord_pr1] = sort(priority_angle_order6(:,1));
                                                                                                                                        priority_angle_order6 = priority_angle_order6(ord_pr1,:);
                                                                                                                                    end
                                                                                                                                    for find_six_loop2 = 1: row_six_size(1,1)  
                                                                                                                                        if(stop_search2 == 0)
                                                                                                                                            pass_numfil=0; 
                                                                                                                                            diff_vector2=0;
                                                                                                                                            posit_nega2=nan;
                                                                                                                                            xc_6=xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1);
                                                                                                                                            yc_6=xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2);
                                                                                                                                            if(angle_final_th <= 45 & angle_final_th >= -45) %vetorc
                                                                                                                                                vector=1; % x vetorc
                                                                                                                                                diff_vector = xc_5-xc_4;
                                                                                                                                                posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                                            else 
                                                                                                                                                vector=0; % y vetorc
                                                                                                                                                diff_vector = yc_5-yc_4;
                                                                                                                                                posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                                            end
                                                                                                                                            if(vector == 1)       
                                                                                                                                                diff_vector2= xc_6 - xc_5;
                                                                                                                                                posit_nega2=sign(diff_vector2);                                           
                                                                                                                                            elseif(vector == 0)                                            
                                                                                                                                                diff_vector2= yc_6 - yc_5;
                                                                                                                                                posit_nega2=sign(diff_vector2);
                                                                                                                                            end
                                                                                                                                            if( posit_nega == 0 )
                                                                                                                                                xc_4=xc_3;
                                                                                                                                                yc_4=yc_3;       
                                                                                                                                                posit_nega = posit_nega2;        
                                                                                                                                            end
                                                                                                                                            if(posit_nega == posit_nega2)
                                                                                                                                                angle_final_th=atan((yc_5-yc_4)/(xc_5-xc_4))*180/pi;
                                                                                                                                                if( (yc_5-yc_4) == 0 & (xc_5-xc_4) == 0 )
                                                                                                                                                    angle_final_th = 0;
                                                                                                                                                end
                                                                                                                                                dist_final_th=sqrt((xc_5-xc_4)^2+(yc_5-yc_4)^2);
                                                                                                                                                increase_check=1; 
                                                                                                                                                if((distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)))+dist_final_th)*increase_distance_threshold < sqrt( (xc_6-xc_4)^2+(yc_6-yc_4)^2 ))
                                                                                                                                                    increase_check=0;                                                     
                                                                                                                                                else            
                                                                                                                                                    increase_check=1;
                                                                                                                                                end  
                                                                                                                                                if(posit_nega2 == 0)
                                                                                                                                                    pass_numfil=1;
                                                                                                                                                    increase_check=0;
                                                                                                                                                end
                                                                                                                                                if(increase_check == 0)
                                                                                                                                                    angle_six_th_cmp= atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi;
                                                                                                                                                    if( (yc_6-yc_5) == 0 & (xc_6-xc_5) == 0 )
                                                                                                                                                        angle_six_th_cmp=0;
                                                                                                                                                    end
                                                                                                                                                    if( dist_final_th > decide_slow_angle)
                                                                                                                                                        if(abs(angle_final_th)>tangent_control | abs(angle_six_th_cmp)>tangent_control)    
                                                                                                                                                            if(abs(angle_final_th)+control_vecter_angle > abs(atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi) & abs(angle_final_th)-control_vecter_angle < abs(atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi))  
                                                                                                                                                                pass_numfil=1;
                                                                                                                                                            end
                                                                                                                                                        else
                                                                                                                                                            if(angle_final_th+control_vecter_angle > atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi & angle_final_th-control_vecter_angle < atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi)  
                                                                                                                                                                pass_numfil=1;
                                                                                                                                                            end
                                                                                                                                                        end
                                                                                                                                                    else
                                                                                                                                                        if(abs(angle_final_th)>tangent_control | abs(angle_six_th_cmp)>tangent_control)
                                                                                                                                                            if(abs(angle_final_th)+control_vecter_angle*multilple_short_angle > abs(atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi) & abs(angle_final_th)-control_vecter_angle*multilple_short_angle < abs(atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi))  
                                                                                                                                                                pass_numfil=1;
                                                                                                                                                            end
                                                                                                                                                        else
                                                                                                                                                            if(angle_final_th+control_vecter_angle*multilple_short_angle > atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi & angle_final_th-control_vecter_angle*multilple_short_angle < atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi)  
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
                                                                                                                                        duplication_number=prevent_duplication_loop;
                                                                                                                                    end
                                                                                                                                end
                                                                                                                                if(duplication_number ~= 0)                                                                                  
                                                                                                                                    if(angle_information(duplication_number,1) < angle_information(object_number_count,1))
                                                                                                                                        duplehappen=1;
                                                                                                                                    end
                                                                                                                                    if( duplehappen == 1)
                                                                                                                                        object_number_count= object_number_count-1; %duplication happen
                                                                                                                                    else       
                                                                                                                                        fast_total_coor(1,duplication_number) = {0};
                                                                                                                                        object_number_count = object_number_count-1;
                                                                                                                                        maxsize_fst_check=size(fst_mv_oj_coor);
                                                                                                                                        for thirdnumberchange=1:maxsize_fst_check(1,1)
                                                                                                                                            fst_mv_oj_coor(thirdnumberchange,3) = duplication_number;
                                                                                                                                        end
                                                                                                                                        fast_total_coor(1,duplication_number) = {fst_mv_oj_coor}; %smaller angle is right trajectory
                                                                                                                                        new_search1=fast_total_coor{1,duplication_number}(1,4);
                                                                                                                                        %****************************************main restore***************************

                                                                                                                                        distance_fast_mat=0;        

                                                                                                                                        [row_f1,col_f1,check_f1]=find( xyz{1,fast_obj_start+1}(:,1) < xyz_1{1,fast_obj_start}(new_search1,1) + distance_threshold_fast & xyz{1,fast_obj_start+1}(:,1) > xyz_1{1,fast_obj_start}(new_search1,1) - distance_threshold_fast);   
                                                                                                                                        xc_1=xyz_1{1,fast_obj_start}(new_search1,1);
                                                                                                                                        yc_1=xyz_1{1,fast_obj_start}(new_search1,2);

                                                                                                                                        if(check_f1 == 1)  
                                                                                                                                            maxsize_fast2_xyz=size(row_f1); %read size
                                                                                                                                            distance_count1=0;
                                                                                                                                            for fast_second_mt=1:maxsize_fast2_xyz(1,1) %for fast sperm
                                                                                                                                                xx_1 = xyz{1,fast_obj_start+1}(row_f1(fast_second_mt),1) - xc_1; % comparing x coordinates   
                                                                                                                                                yy_1 = xyz{1,fast_obj_start+1}(row_f1(fast_second_mt),2) - yc_1; % comparing y coordinates
                                                                                                                                                distance_count1=distance_count1+1;
                                                                                                                                                distance_fast_mat(distance_count1,1) = sqrt((xx_1^2)+(yy_1^2));
                                                                                                                                                distance_fast_mat(distance_count1,2) = row_f1(fast_second_mt);
                                                                                                                                            end
                                                                                                                                            [row_dist,col_dist,check_dist]=find(distance_fast_mat(:,1) <= distance_threshold_fast);     
                                                                                                                                            if (check_dist == 1) %first distance threshold check
                                                                                                                                                row_dist_size=size(row_dist); %read size  
                                                                                                                                                angle_fast=0;
                                                                                                                                                dist_thresh_sec=0;
                                                                                                                                                priority_angle_order=0;
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
                                                                                                                                                for try_loop2 = 1: row_dist_size(1,1)
                                                                                                                                                    if(stop_search == 0) 
                                                                                                                                                        pass_numfil=0;
                                                                                                                                                        xc_2=xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1);
                                                                                                                                                        yc_2=xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2);
                                                                                                                                                        %decide vector andgle
                                                                                                                                                        %each magnitude
                                                                                                                                                        angle_fast = atan((yc_2-yc_1)/(xc_2-xc_1))*180/pi; % unit is degree
                                                                                                                                                        if( (yc_2-yc_1) == 0 & (xc_2-xc_1) == 0 )
                                                                                                                                                            angle_fast = 0;
                                                                                                                                                        end
                                                                                                                                                        dist_thresh_sec = distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)));
                                                                                                                                                        if(angle_fast <= 45 & angle_fast >= -45) %vetorc
                                                                                                                                                            vector=1; % x vetorc
                                                                                                                                                            diff_vector = xc_2 - xc_1;
                                                                                                                                                            posit_nega=sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                                                        else
                                                                                                                                                            vector=0; % y vetorc
                                                                                                                                                            diff_vector = yc_2 - yc_1; 
                                                                                                                                                            posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                                                        end
                                                                                                                                                        posit_nega2 = posit_nega;
                                                                                                                                                        pass_numfil=1;                      
                                                                                                                                                        increase_check=0;      
                                                                                                                                                        if(increase_check == 0)       
                                                                                                                                                            if(pass_numfil == 1)   
                                                                                                                                                                distance_repeat_rest=0;
                                                                                                                                                                [row_f2,col_f2,check_f2]=find( xyz{1,fast_obj_start+2}(:,1) < xc_2 + distance_threshold_fast & xyz{1,fast_obj_start+2}(:,1) > xc_2 - distance_threshold_fast);       

                                                                                                                                                                if(check_f2 == 1)
                                                                                                                                                                    maxsize_fast_xyz2=size(row_f2); %read size
                                                                                                                                                                    distance_count2=0;
                                                                                                                                                                    for dist_cmp_loop = 1: maxsize_fast_xyz2(1,1)
                                                                                                                                                                        xx_2 = xyz{1,fast_obj_start+2}(row_f2(dist_cmp_loop),1) - xc_2;                              
                                                                                                                                                                        yy_2 = xyz{1,fast_obj_start+2}(row_f2(dist_cmp_loop),2) - yc_2;                 
                                                                                                                                                                        distance_count2=distance_count2+1;
                                                                                                                                                                        distance_repeat_rest(distance_count2,1) = sqrt((xx_2^2)+(yy_2^2));
                                                                                                                                                                        distance_repeat_rest(distance_count2,2) = row_f2(dist_cmp_loop);
                                                                                                                                                                    end 
                                                                                                                                                                    if(dist_thresh_sec <= decide_slow)
                                                                                                                                                                        [row_secd,col_secd,check_secd]=find(distance_repeat_rest(:,1) < distance_threshold_fast/2);
                                                                                                                                                                    else
                                                                                                                                                                        [row_secd,col_secd,check_secd]=find(distance_repeat_rest(:,1) < dist_thresh_sec*threshold_diff_dislocation & distance_repeat_rest(:,1) < distance_threshold_fast);
                                                                                                                                                                    end
                                                                                                                                                                    if (check_secd == 1)
                                                                                                                                                                        row_secd_size=size(row_secd); %read size  
                                                                                                                                                                        angle_secd_fast=0;
                                                                                                                                                                        dist_thresh_third=0;
                                                                                                                                                                        priority_angle_order2=0;
                                                                                                                                                                        for find_secd_loop = 1: row_secd_size(1,1)
                                                                                                                                                                            angle_secd_fast_test = atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),2)-yc_1)/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),1)-xc_1))*180/pi;
                                                                                                                                                                            priority_angle_order2(find_secd_loop,1)=abs(abs(angle_secd_fast_test)-abs(angle_fast))/priority_degree_rate + distance_repeat_rest(row_secd(find_secd_loop));
                                                                                                                                                                            if( (xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),2)-yc_1) == 0 & (xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),1)-xc_1) == 0 )
                                                                                                                                                                                priority_angle_order2(find_secd_loop,1) = priority_zero; 
                                                                                                                                                                            end
                                                                                                                                                                            priority_angle_order2(find_secd_loop,2)=find_secd_loop;%give priority for smaller angle!!
                                                                                                                                                                            [temp_pr1,ord_pr1] = sort(priority_angle_order2(:,1));
                                                                                                                                                                            priority_angle_order2 = priority_angle_order2(ord_pr1,:);
                                                                                                                                                                        end
                                                                                                                                                                        for find_secd_loop2 = 1: row_secd_size(1,1)
                                                                                                                                                                            if(stop_search == 0)
                                                                                                                                                                                xc_3=xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1);
                                                                                                                                                                                yc_3=xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2);
                                                                                                                                                                                pass_numfil=0;
                                                                                                                                                                                diff_vector2=0;
                                                                                                                                                                                posit_nega2=nan;

                                                                                                                                                                                if(vector == 1)
                                                                                                                                                                                    diff_vector2= xc_3-xc_2;
                                                                                                                                                                                    posit_nega2=sign(diff_vector2);
                                                                                                                                                                                elseif (vector == 0)
                                                                                                                                                                                    diff_vector2= yc_3-yc_2;                            
                                                                                                                                                                                    posit_nega2=sign(diff_vector2);
                                                                                                                                                                                end
                                                                                                                                                                                if( posit_nega == 0 )
                                                                                                                                                                                    if( fast_obj_start > 1+change_start_point)
                                                                                                                                                                                        xc_1=xc_back_1;
                                                                                                                                                                                        yc_1=yc_back_1;
                                                                                                                                                                                    end
                                                                                                                                                                                    posit_nega = posit_nega2;        
                                                                                                                                                                                end
                                                                                                                                                                                if(posit_nega == posit_nega2)
                                                                                                                                                                                    angle_fast=atan((yc_2-yc_1)/(xc_2-xc_1))*180/pi;
                                                                                                                                                                                    if( (yc_2-yc_1) == 0 & (xc_2-xc_1) == 0 )
                                                                                                                                                                                        angle_fast=0;
                                                                                                                                                                                    end
                                                                                                                                                                                    dist_thresh_sec=sqrt((xc_2-xc_1)^2+(yc_2-yc_1)^2);
                                                                                                                                                                                    increase_check=1;

                                                                                                                                                                                    if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xc_3-xc_1)^2+(yc_3-yc_1)^2 ))
                                                                                                                                                                                        increase_check=0;
                                                                                                                                                                                    else
                                                                                                                                                                                        increase_check=1;
                                                                                                                                                                                    end

                                                                                                                                                                                    if(posit_nega2 == 0)
                                                                                                                                                                                        pass_numfil=1;
                                                                                                                                                                                        increase_check=0;
                                                                                                                                                                                    end
                                                                                                                                                                                    if(increase_check == 0)
                                                                                                                                                                                        angle_secd_fast_cmp=atan((yc_3-yc_2)/(xc_3-xc_2))*180/pi; 
                                                                                                                                                                                        if( (yc_3-yc_2) == 0 & (xc_3-xc_2) == 0 )
                                                                                                                                                                                            angle_secd_fast_cmp=0;
                                                                                                                                                                                        end
                                                                                                                                                                                        if(dist_thresh_sec > decide_slow_angle)
                                                                                                                                                                                            if(abs(angle_fast)>tangent_control | abs(angle_secd_fast_cmp) > tangent_control)                                                     
                                                                                                                                                                                                if(abs(angle_fast)+control_vecter_angle > abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) & abs(angle_fast)-control_vecter_angle < abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi)) % unit is degree                       
                                                                                                                                                                                                    pass_numfil=1;
                                                                                                                                                                                                end
                                                                                                                                                                                            else
                                                                                                                                                                                                if(angle_fast+control_vecter_angle > atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi & angle_fast-control_vecter_angle < atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) % unit is degree  
                                                                                                                                                                                                    pass_numfil=1;
                                                                                                                                                                                                end
                                                                                                                                                                                            end
                                                                                                                                                                                        else
                                                                                                                                                                                            if(abs(angle_fast)>tangent_control | abs(angle_secd_fast_cmp) > tangent_control)                                                     
                                                                                                                                                                                                if(abs(angle_fast)+control_vecter_angle*multilple_short_angle > abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) & abs(angle_fast)-control_vecter_angle*multilple_short_angle < abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi)) % unit is degree                       
                                                                                                                                                                                                    pass_numfil=1;
                                                                                                                                                                                                end
                                                                                                                                                                                            else
                                                                                                                                                                                                if(angle_fast+control_vecter_angle*multilple_short_angle > atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi & angle_fast-control_vecter_angle*multilple_short_angle < atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) % unit is degree  
                                                                                                                                                                                                    pass_numfil=1;
                                                                                                                                                                                                end
                                                                                                                                                                                            end
                                                                                                                                                                                        end
                                                                                                                                                                                        if(pass_numfil == 1) 
                                                                                                                                                                                            angle_secd_fast = atan((yc_3-yc_2)/(xc_3-xc_2))*180/pi;   
                                                                                                                                                                                            if( (yc_3-yc_2) == 0 & (xc_3-xc_2) == 0 )
                                                                                                                                                                                                angle_secd_fast=0;
                                                                                                                                                                                            end
                                                                                                                                                                                            dist_thresh_third = distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)));                           
                                                                                                                                                                                            % finish third
                                                                                                                                                                                            distance_thrd_mat=0;
                                                                                                                                                                                            [row_f3,col_f3,check_f3]=find( xyz{1,fast_obj_start+3}(:,1) < xc_3 + distance_threshold_fast & xyz{1,fast_obj_start+3}(:,1) > xc_3 - distance_threshold_fast);       

                                                                                                                                                                                            if(check_f3 == 1)
                                                                                                                                                                                                maxsize_fast_xyz3=size(row_f3);
                                                                                                                                                                                                distance_count3=0;
                                                                                                                                                                                                for dist_thrd_loop = 1: maxsize_fast_xyz3(1,1)                                   
                                                                                                                                                                                                    xx_3 = xyz{1,fast_obj_start+3}(row_f3(dist_thrd_loop),1) - xc_3;
                                                                                                                                                                                                    yy_3 = xyz{1,fast_obj_start+3}(row_f3(dist_thrd_loop),2) - yc_3;
                                                                                                                                                                                                    distance_count3=distance_count3+1;
                                                                                                                                                                                                    distance_thrd_mat(distance_count3,1) = sqrt((xx_3^2)+(yy_3^2));
                                                                                                                                                                                                    distance_thrd_mat(distance_count3,2) = row_f3(dist_thrd_loop);
                                                                                                                                                                                                end
                                                                                                                                                                                                if(dist_thresh_third <= decide_slow)
                                                                                                                                                                                                    [row_third,col_third,check_third]=find(distance_thrd_mat(:,1) < distance_threshold_fast/2);
                                                                                                                                                                                                else
                                                                                                                                                                                                    [row_third,col_third,check_third]=find(distance_thrd_mat(:,1) < dist_thresh_third*threshold_diff_dislocation & distance_thrd_mat(:,1) <= distance_threshold_fast );
                                                                                                                                                                                                end
                                                                                                                                                                                                if (check_third == 1)
                                                                                                                                                                                                    row_third_size=size(row_third);
                                                                                                                                                                                                    angle_rest_th=0;
                                                                                                                                                                                                    dist_rest_th=0;                         
                                                                                                                                                                                                    priority_angle_order3=0;
                                                                                                                                                                                                    for find_rest_loop = 1: row_third_size(1,1)
                                                                                                                                                                                                        angle_rest_th_test=atan((xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),2)-yc_2)/(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),1)-xc_2))*180/pi;                                                           
                                                                                                                                                                                                        priority_angle_order3(find_rest_loop,1)=abs(abs(angle_rest_th_test)-abs(angle_secd_fast))/priority_degree_rate + distance_thrd_mat(row_third(find_rest_loop));
                                                                                                                                                                                                        if( (xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),2)-yc_2) == 0 & (xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),1)-xc_2) ==0 )
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
                                                                                                                                                                                                            xc_4=xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1);
                                                                                                                                                                                                            yc_4=xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2);
                                                                                                                                                                                                            if(angle_secd_fast <= 45 & angle_secd_fast >= -45) %vetorc
                                                                                                                                                                                                                vector=1; % x vetorc
                                                                                                                                                                                                                diff_vector = xc_3-xc_2;
                                                                                                                                                                                                                posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                                                                                                            else 
                                                                                                                                                                                                                vector=0; % y vetorc
                                                                                                                                                                                                                diff_vector = yc_3-yc_2;
                                                                                                                                                                                                                posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                                                                                                            end
                                                                                                                                                                                                            if(vector == 1)       
                                                                                                                                                                                                                diff_vector2= xc_4 - xc_3;
                                                                                                                                                                                                                posit_nega2=sign(diff_vector2);                                           
                                                                                                                                                                                                            elseif(vector == 0)                                            
                                                                                                                                                                                                                diff_vector2= yc_4 - yc_3;                                            
                                                                                                                                                                                                                posit_nega2=sign(diff_vector2);
                                                                                                                                                                                                            end
                                                                                                                                                                                                            if( posit_nega == 0 )
                                                                                                                                                                                                                xc_2=xc_1;
                                                                                                                                                                                                                yc_2=yc_1;
                                                                                                                                                                                                                posit_nega = posit_nega2;        
                                                                                                                                                                                                            end
                                                                                                                                                                                                            if(posit_nega == posit_nega2)
                                                                                                                                                                                                                angle_secd_fast=atan((yc_3-yc_2)/(xc_3-xc_2))*180/pi;
                                                                                                                                                                                                                if( (yc_3-yc_2) == 0 & (xc_3-xc_2) == 0 )
                                                                                                                                                                                                                    angle_secd_fast = 0;
                                                                                                                                                                                                                end
                                                                                                                                                                                                                dist_thresh_third=sqrt((xc_3-xc_2)^2+(yc_3-yc_2)^2);
                                                                                                                                                                                                                increase_check=1;
                                                                                                                                                                                                                if((distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)))+dist_thresh_third)*increase_distance_threshold < sqrt( (xc_4-xc_2)^2+(yc_4-yc_2)^2 ))
                                                                                                                                                                                                                    increase_check=0;                                                     
                                                                                                                                                                                                                else            
                                                                                                                                                                                                                    increase_check=1;
                                                                                                                                                                                                                end 
                                                                                                                                                                                                                if(posit_nega2 == 0)
                                                                                                                                                                                                                    pass_numfil=1;
                                                                                                                                                                                                                    increase_check=0;
                                                                                                                                                                                                                end
                                                                                                                                                                                                                if(increase_check == 0)  
                                                                                                                                                                                                                    angle_rest_th_cmp=atan((yc_4-yc_3)/(xc_4-xc_3))*180/pi;
                                                                                                                                                                                                                    if( (yc_4-yc_3) == 0 & (xc_4-xc_3)==0)
                                                                                                                                                                                                                        angle_rest_th_cmp=0;
                                                                                                                                                                                                                    end
                                                                                                                                                                                                                    if( dist_thresh_third > decide_slow_angle)
                                                                                                                                                                                                                        if(abs(angle_secd_fast)>tangent_control | abs(angle_rest_th_cmp)>tangent_control)    
                                                                                                                                                                                                                            if(abs(angle_secd_fast)+control_vecter_angle > abs(atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi) & abs(angle_secd_fast)-control_vecter_angle < abs(atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi))  
                                                                                                                                                                                                                                pass_numfil=1;
                                                                                                                                                                                                                            end
                                                                                                                                                                                                                        else
                                                                                                                                                                                                                            if(angle_secd_fast+control_vecter_angle > atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi & angle_secd_fast-control_vecter_angle < atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi)  
                                                                                                                                                                                                                                pass_numfil=1;
                                                                                                                                                                                                                            end
                                                                                                                                                                                                                        end
                                                                                                                                                                                                                    else
                                                                                                                                                                                                                        if(abs(angle_secd_fast)>tangent_control | abs(angle_rest_th_cmp)>tangent_control)    
                                                                                                                                                                                                                            if(abs(angle_secd_fast)+control_vecter_angle*multilple_short_angle > abs(atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi) & abs(angle_secd_fast)-control_vecter_angle*multilple_short_angle < abs(atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi))  
                                                                                                                                                                                                                                pass_numfil=1;
                                                                                                                                                                                                                            end
                                                                                                                                                                                                                        else
                                                                                                                                                                                                                            if(angle_secd_fast+control_vecter_angle*multilple_short_angle > atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi & angle_secd_fast-control_vecter_angle*multilple_short_angle < atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi)  
                                                                                                                                                                                                                                pass_numfil=1;
                                                                                                                                                                                                                            end
                                                                                                                                                                                                                        end 
                                                                                                                                                                                                                    end
                                                                                                                                                                                                                    if( pass_numfil == 1) 
                                                                                                                                                                                                                        angle_rest_th=atan((yc_4-yc_3)/(xc_4-xc_3))*180/pi;                 
                                                                                                                                                                                                                        if ( (yc_4-yc_3) == 0 & (xc_4-xc_3) == 0 )
                                                                                                                                                                                                                            angle_rest_th = 0;
                                                                                                                                                                                                                        end
                                                                                                                                                                                                                        dist_rest_th=distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)));
                                                                                                                                                                                                                        [row_f4,col_f4,check_f4]=find( xyz{1,fast_obj_start+4}(:,1) < xc_4 + distance_threshold_fast & xyz{1,fast_obj_start+4}(:,1) > xc_4 - distance_threshold_fast);       

                                                                                                                                                                                                                        if(check_f4 == 1)
                                                                                                                                                                                                                            distance_four_mat=0;
                                                                                                                                                                                                                            maxsize_fast_xyz4=size(row_f4);
                                                                                                                                                                                                                            distance_count4=0;
                                                                                                                                                                                                                            for dist_four_loop = 1: maxsize_fast_xyz4(1,1)
                                                                                                                                                                                                                                xx_4 = xyz{1,fast_obj_start+4}(row_f4(dist_four_loop),1) - xc_4;    
                                                                                                                                                                                                                                yy_4 = xyz{1,fast_obj_start+4}(row_f4(dist_four_loop),2) - yc_4;                                
                                                                                                                                                                                                                                distance_count4=distance_count4+1;
                                                                                                                                                                                                                                distance_four_mat(distance_count4,1) = sqrt((xx_4^2)+(yy_4^2));
                                                                                                                                                                                                                                distance_four_mat(distance_count4,2) = row_f4(dist_four_loop);
                                                                                                                                                                                                                            end
                                                                                                                                                                                                                            if(dist_rest_th <= decide_slow)     
                                                                                                                                                                                                                                [row_four,col_four,check_four]=find(distance_four_mat(:,1) < distance_threshold_fast/2);                              
                                                                                                                                                                                                                            else
                                                                                                                                                                                                                                [row_four,col_four,check_four]=find(distance_four_mat(:,1) < dist_rest_th*threshold_diff_dislocation & distance_four_mat(:,1) <= distance_threshold_fast);   
                                                                                                                                                                                                                            end
                                                                                                                                                                                                                            if (check_four == 1)
                                                                                                                                                                                                                                row_four_size=size(row_four);
                                                                                                                                                                                                                                angle_final_th=0;
                                                                                                                                                                                                                                dist_final_th=0; 
                                                                                                                                                                                                                                priority_angle_order4=0;
                                                                                                                                                                                                                                for find_four_loop = 1: row_four_size(1,1)
                                                                                                                                                                                                                                    angle_final_th_test=atan((xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),2)-yc_3)/(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),1)-xc_3))*180/pi;
                                                                                                                                                                                                                                    priority_angle_order4(find_four_loop,1)=abs(abs(angle_rest_th)-abs(angle_final_th_test))/priority_degree_rate + distance_four_mat(row_four(find_four_loop));
                                                                                                                                                                                                                                    if( (xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),2)-yc_3) == 0 & (xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),1)-xc_3) == 0)
                                                                                                                                                                                                                                        priority_angle_order4(find_four_loop,1)=priority_zero;
                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                    priority_angle_order4(find_four_loop,2)=find_four_loop;%give priority for smaller angle!!
                                                                                                                                                                                                                                    [temp_pr1,ord_pr1] = sort(priority_angle_order4(:,1));
                                                                                                                                                                                                                                    priority_angle_order4 = priority_angle_order4(ord_pr1,:);
                                                                                                                                                                                                                                end
                                                                                                                                                                                                                                for find_four_loop2 = 1: row_four_size(1,1)                                                          
                                                                                                                                                                                                                                    if(stop_search == 0)
                                                                                                                                                                                                                                        xc_5=xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1);
                                                                                                                                                                                                                                        yc_5=xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2);
                                                                                                                                                                                                                                        pass_numfil=0;
                                                                                                                                                                                                                                        diff_vector2=0;
                                                                                                                                                                                                                                        posit_nega2=nan;
                                                                                                                                                                                                                                        if(angle_rest_th <= 45 & angle_rest_th >= -45) %vetorc
                                                                                                                                                                                                                                            vector=1; % x vetorc
                                                                                                                                                                                                                                            diff_vector = xc_4-xc_3;
                                                                                                                                                                                                                                            posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                                                                                                                                        else 
                                                                                                                                                                                                                                            vector=0; % y vetorc
                                                                                                                                                                                                                                            diff_vector = yc_4-yc_3;
                                                                                                                                                                                                                                            posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                        if(vector == 1)       
                                                                                                                                                                                                                                            diff_vector2= xc_5 - xc_4;      
                                                                                                                                                                                                                                            posit_nega2=sign(diff_vector2);        
                                                                                                                                                                                                                                        elseif(vector == 0)                                            
                                                                                                                                                                                                                                            diff_vector2= yc_5 - yc_4;                                            
                                                                                                                                                                                                                                            posit_nega2=sign(diff_vector2);
                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                        if( posit_nega == 0 )
                                                                                                                                                                                                                                            xc_3=xc_2;
                                                                                                                                                                                                                                            yc_3=yc_2;
                                                                                                                                                                                                                                            posit_nega = posit_nega2;        
                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                        if(posit_nega == posit_nega2)
                                                                                                                                                                                                                                            angle_rest_th=atan((yc_4-yc_3)/(xc_4-xc_3))*180/pi;
                                                                                                                                                                                                                                            if( (yc_4-yc_3) == 0 & (xc_4-xc_3) == 0 )
                                                                                                                                                                                                                                                angle_rest_th = 0;
                                                                                                                                                                                                                                            end
                                                                                                                                                                                                                                            dist_rest_th= sqrt((xc_4-xc_3)^2+(yc_4-yc_3)^2);
                                                                                                                                                                                                                                            increase_check=1;
                                                                                                                                                                                                                                            if((distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)))+dist_rest_th)*increase_distance_threshold < sqrt((xc_5-xc_3)^2+(yc_5-yc_3)^2))       
                                                                                                                                                                                                                                                increase_check=0;    
                                                                                                                                                                                                                                            else
                                                                                                                                                                                                                                                increase_check=1;      
                                                                                                                                                                                                                                            end
                                                                                                                                                                                                                                            if(posit_nega2 == 0)
                                                                                                                                                                                                                                                pass_numfil=1;
                                                                                                                                                                                                                                                increase_check=0;
                                                                                                                                                                                                                                            end
                                                                                                                                                                                                                                            if(increase_check == 0)                                                                               
                                                                                                                                                                                                                                                angle_final_th_cmp=atan((yc_5-yc_4)/(xc_5-xc_4))*180/pi;
                                                                                                                                                                                                                                                if( (yc_5-yc_4) == 0 & (xc_5-xc_4) ==0 )
                                                                                                                                                                                                                                                    angle_final_th_cmp=0;
                                                                                                                                                                                                                                                end
                                                                                                                                                                                                                                                if( dist_rest_th > decide_slow_angle)
                                                                                                                                                                                                                                                    if(abs(angle_rest_th)>tangent_control | abs(angle_final_th_cmp) >tangent_control)
                                                                                                                                                                                                                                                        if(abs(angle_rest_th)+control_vecter_angle > abs(atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi) & abs(angle_rest_th)-control_vecter_angle < abs(atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi))
                                                                                                                                                                                                                                                            pass_numfil=1;
                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                    else
                                                                                                                                                                                                                                                        if(angle_rest_th+control_vecter_angle > atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi & angle_rest_th-control_vecter_angle < atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi)
                                                                                                                                                                                                                                                            pass_numfil=1;
                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                else  
                                                                                                                                                                                                                                                     if(abs(angle_rest_th)>tangent_control | abs(angle_final_th_cmp) >tangent_control)
                                                                                                                                                                                                                                                        if(abs(angle_rest_th)+control_vecter_angle*multilple_short_angle > abs(atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi) & abs(angle_rest_th)-control_vecter_angle*multilple_short_angle  < abs(atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi))
                                                                                                                                                                                                                                                            pass_numfil=1;
                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                    else
                                                                                                                                                                                                                                                        if(angle_rest_th+control_vecter_angle*multilple_short_angle  > atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi & angle_rest_th-control_vecter_angle*multilple_short_angle  < atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi)
                                                                                                                                                                                                                                                            pass_numfil=1;
                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                end
                                                                                                                                                                                                                                                if( pass_numfil == 1)
                                                                                                                                                                                                                                                    angle_final_th=atan((yc_5-yc_4)/(xc_5-xc_4))*180/pi;                                    
                                                                                                                                                                                                                                                    if( (yc_5-yc_4) == 0 & (xc_5-xc_4)==0)
                                                                                                                                                                                                                                                        angle_final_th=0;
                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                    dist_final_th=distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)));

                                                                                                                                                                                                                                                    object_number_count = object_number_count+1;
                                                                                                                                                                                                                                                    fst_mv_oj_coor=0;

                                                                                                                                                                                                                                                    fst_mv_oj_coor(1,1) = xyz{1,fast_obj_start}(new_search1,1);                              
                                                                                                                                                                                                                                                    fst_mv_oj_coor(1,2) = xyz{1,fast_obj_start}(new_search1,2);          
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
                                                                                                                                                                                                                                                        angle_information(object_number_count,1)=(abs(abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(1,1)))*180/pi)-abs(atan((fst_mv_oj_coor(2,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(2,1)-fst_mv_oj_coor(1,1)))*180/pi))+abs(abs(atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(2,1)))*180/pi)-abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(2,1)))*180/pi))+abs(abs(atan((fst_mv_oj_coor(5,2)-fst_mv_oj_coor(3,2))/(fst_mv_oj_coor(5,1)-fst_mv_oj_coor(3,1)))*180/pi)-abs(atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(3,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(3,1)))*180/pi)))/priority_degree_rate + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) );
                                                                                                                                                                                                                                                    else
                                                                                                                                                                                                                                                        angle_information(object_number_count,1)=(abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(1,1)))*180/pi-atan((fst_mv_oj_coor(2,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(2,1)-fst_mv_oj_coor(1,1)))*180/pi)+abs(atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(2,1)))*180/pi-atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(2,1)))*180/pi)+abs(atan((fst_mv_oj_coor(5,2)-fst_mv_oj_coor(3,2))/(fst_mv_oj_coor(5,1)-fst_mv_oj_coor(3,1)))*180/pi-atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(3,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(3,1)))*180/pi))/priority_degree_rate + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) );
                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                    distance_six_mat=0;
                                                                                                                                                                                                                                                    [row_f5,col_f5,check_f5]=find( xyz{1,fast_obj_start+4}(:,1) < xc_5 + distance_threshold_fast & xyz{1,fast_obj_start+4}(:,1) > xc_5 - distance_threshold_fast);       

                                                                                                                                                                                                                                                    if(check_f5 == 1)
                                                                                                                                                                                                                                                        stop_search2=0;
                                                                                                                                                                                                                                                        maxsize_fast_xyz6=size(row_f5);
                                                                                                                                                                                                                                                        distance_count6=0;
                                                                                                                                                                                                                                                        for dist_six_loop = 1: maxsize_fast_xyz6(1,1)                                   
                                                                                                                                                                                                                                                            xx_6 = xyz{1,fast_obj_start+4}(row_f5(dist_six_loop),1) - xc_5;
                                                                                                                                                                                                                                                            yy_6 = xyz{1,fast_obj_start+4}(row_f5(dist_six_loop),2) - yc_5;                                                                                                                                       
                                                                                                                                                                                                                                                            distance_count6=distance_count6+1;
                                                                                                                                                                                                                                                            distance_six_mat(distance_count6,1) = sqrt((xx_6^2)+(yy_6^2));
                                                                                                                                                                                                                                                            distance_six_mat(distance_count6,2) = row_f5(dist_six_loop);
                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                        if(dist_final_th <= decide_slow)
                                                                                                                                                                                                                                                            [row_six,col_six,check_six]=find(distance_six_mat(:,1) <= distance_threshold_fast/2 & distance_six_mat(:,1) ~= 0);
                                                                                                                                                                                                                                                        else
                                                                                                                                                                                                                                                            [row_six,col_six,check_six]=find(distance_six_mat(:,1) <= dist_final_th*threshold_diff_dislocation & distance_six_mat(:,1) <= distance_threshold_fast & distance_six_mat(:,1) ~= 0);
                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                        if (check_six == 1)
                                                                                                                                                                                                                                                            row_six_size=size(row_six);
                                                                                                                                                                                                                                                            angle_six_th=0;
                                                                                                                                                                                                                                                            dist_six_th=0;                         
                                                                                                                                                                                                                                                            priority_angle_order6=0;
                                                                                                                                                                                                                                                            for find_six_loop = 1: row_six_size(1,1)
                                                                                                                                                                                                                                                                angle_six_th_test=atan((xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),2)-yc_4)/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),1)-xc_4))*180/pi;                                                           
                                                                                                                                                                                                                                                                priority_angle_order6(find_six_loop,1)=abs(abs(angle_final_th)-abs(angle_six_th_test))/priority_degree_rate + distance_six_mat(row_six(find_six_loop));
                                                                                                                                                                                                                                                                if( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),2)-yc_4) == 0 & (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),1)- xc_4) )
                                                                                                                                                                                                                                                                    priority_angle_order6(find_six_loop,1)= priority_zero; 
                                                                                                                                                                                                                                                                end
                                                                                                                                                                                                                                                                priority_angle_order6(find_six_loop,2)=find_six_loop;%give priority for smaller angle!!
                                                                                                                                                                                                                                                                [temp_pr1,ord_pr1] = sort(priority_angle_order6(:,1));
                                                                                                                                                                                                                                                                priority_angle_order6 = priority_angle_order6(ord_pr1,:);
                                                                                                                                                                                                                                                            end
                                                                                                                                                                                                                                                            for find_six_loop2 = 1: row_six_size(1,1)  
                                                                                                                                                                                                                                                                if(stop_search2 == 0)
                                                                                                                                                                                                                                                                    pass_numfil=0; 
                                                                                                                                                                                                                                                                    diff_vector2=0;
                                                                                                                                                                                                                                                                    posit_nega2=nan;
                                                                                                                                                                                                                                                                    xc_6=xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1);
                                                                                                                                                                                                                                                                    yc_6=xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2);
                                                                                                                                                                                                                                                                    if(angle_final_th <= 45 & angle_final_th >= -45) %vetorc
                                                                                                                                                                                                                                                                        vector=1; % x vetorc
                                                                                                                                                                                                                                                                        diff_vector = xc_5-xc_4;
                                                                                                                                                                                                                                                                        posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                                                                                                                                                                    else 
                                                                                                                                                                                                                                                                        vector=0; % y vetorc
                                                                                                                                                                                                                                                                        diff_vector = yc_5-yc_4;
                                                                                                                                                                                                                                                                        posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                                    if(vector == 1)       
                                                                                                                                                                                                                                                                        diff_vector2= xc_6 - xc_5;
                                                                                                                                                                                                                                                                        posit_nega2=sign(diff_vector2);                                           
                                                                                                                                                                                                                                                                    elseif(vector == 0)                                            
                                                                                                                                                                                                                                                                        diff_vector2= yc_6 - yc_5;
                                                                                                                                                                                                                                                                        posit_nega2=sign(diff_vector2);
                                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                                    if( posit_nega == 0 )
                                                                                                                                                                                                                                                                        xc_4=xc_3;
                                                                                                                                                                                                                                                                        yc_4=yc_3;       
                                                                                                                                                                                                                                                                        posit_nega = posit_nega2;        
                                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                                    if(posit_nega == posit_nega2)
                                                                                                                                                                                                                                                                        angle_final_th=atan((yc_5-yc_4)/(xc_5-xc_4))*180/pi;
                                                                                                                                                                                                                                                                        if( (yc_5-yc_4) == 0 & (xc_5-xc_4) == 0 )
                                                                                                                                                                                                                                                                            angle_final_th = 0;
                                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                                        dist_final_th=sqrt((xc_5-xc_4)^2+(yc_5-yc_4)^2);
                                                                                                                                                                                                                                                                        increase_check=1; 
                                                                                                                                                                                                                                                                        if((distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)))+dist_final_th)*increase_distance_threshold < sqrt( (xc_6-xc_4)^2+(yc_6-yc_4)^2 ))
                                                                                                                                                                                                                                                                            increase_check=0;                                                     
                                                                                                                                                                                                                                                                        else            
                                                                                                                                                                                                                                                                            increase_check=1;
                                                                                                                                                                                                                                                                        end  
                                                                                                                                                                                                                                                                        if(posit_nega2 == 0)
                                                                                                                                                                                                                                                                            pass_numfil=1;
                                                                                                                                                                                                                                                                            increase_check=0;
                                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                                        if(increase_check == 0)
                                                                                                                                                                                                                                                                            angle_six_th_cmp= atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi;
                                                                                                                                                                                                                                                                            if( (yc_6-yc_5) == 0 & (xc_6-xc_5) == 0 )
                                                                                                                                                                                                                                                                                angle_six_th_cmp=0;
                                                                                                                                                                                                                                                                            end
                                                                                                                                                                                                                                                                            if( dist_final_th > decide_slow_angle)
                                                                                                                                                                                                                                                                                if(abs(angle_final_th)>tangent_control | abs(angle_six_th_cmp)>tangent_control)    
                                                                                                                                                                                                                                                                                    if(abs(angle_final_th)+control_vecter_angle > abs(atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi) & abs(angle_final_th)-control_vecter_angle < abs(atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi))  
                                                                                                                                                                                                                                                                                        pass_numfil=1;
                                                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                                                else
                                                                                                                                                                                                                                                                                    if(angle_final_th+control_vecter_angle > atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi & angle_final_th-control_vecter_angle < atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi)  
                                                                                                                                                                                                                                                                                        pass_numfil=1;
                                                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                                                end
                                                                                                                                                                                                                                                                            else
                                                                                                                                                                                                                                                                                if(abs(angle_final_th)>tangent_control | abs(angle_six_th_cmp)>tangent_control)
                                                                                                                                                                                                                                                                                    if(abs(angle_final_th)+control_vecter_angle*multilple_short_angle > abs(atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi) & abs(angle_final_th)-control_vecter_angle*multilple_short_angle < abs(atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi))  
                                                                                                                                                                                                                                                                                        pass_numfil=1;
                                                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                                                else
                                                                                                                                                                                                                                                                                    if(angle_final_th+control_vecter_angle*multilple_short_angle > atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi & angle_final_th-control_vecter_angle*multilple_short_angle < atan((yc_6-yc_5)/(xc_6-xc_5))*180/pi)  
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
                                                                                                                                                                                                                                                            if(duplication_happen > 1) 
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

                                                                                                                                        %-----finish restore---------------------
                                                                                                                                        stop_search=1;
                                                                                                                                    end
                                                                                                                                else
                                                                                                                                    fast_total_coor(1,object_number_count)={fst_mv_oj_coor};
                                                                                                                                    stop_search=1; %prevent duplicate trajectories
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
                                                                                                                                angle_information(object_number_count,1)=(abs(abs(atan((fst_mv_oj_coor(2,2)-fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,2))/(fst_mv_oj_coor(2,1)-fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,1)))*180/pi)-abs(atan((fst_mv_oj_coor(1,2)-fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,2))/(fst_mv_oj_coor(1,1)-fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,1)))*180/pi))+abs(abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(1,1)))*180/pi)-abs(atan((fst_mv_oj_coor(2,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(2,1)-fst_mv_oj_coor(1,1)))*180/pi))+abs(abs(atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(2,1)))*180/pi)-abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(2,1)))*180/pi)))/priority_degree_rate + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) );
                                                                                                                            else
                                                                                                                                angle_information(object_number_count,1)=(abs(atan((fst_mv_oj_coor(2,2)-fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,2))/(fst_mv_oj_coor(2,1)-fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,1)))*180/pi-atan((fst_mv_oj_coor(1,2)-fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,2))/(fst_mv_oj_coor(1,1)-fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,1)))*180/pi)+abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(1,1)))*180/pi-atan((fst_mv_oj_coor(2,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(2,1)-fst_mv_oj_coor(1,1)))*180/pi)+abs(atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(2,1)))*180/pi-atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(2,1)))*180/pi))/priority_degree_rate + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) );
                                                                                                                            end
                                                                                                                            distance_six_mat=0;
                                                                                                                            [row_f5,col_f5,check_f5]=find( xyz{1,fast_obj_start+4}(:,1) < xc_5 + distance_threshold_fast & xyz{1,fast_obj_start+4}(:,1) > xc_5 - distance_threshold_fast);       

                                                                                                                            if(check_f5 == 1)
                                                                                                                                stop_search2=0;
                                                                                                                                maxsize_fast_xyz6=size(row_f5);
                                                                                                                                distance_count6=0;
                                                                                                                                for dist_six_loop = 1: maxsize_fast_xyz6(1,1)                                   
                                                                                                                                    xx_6 = xyz{1,fast_obj_start+4}(row_f5(dist_six_loop),1) - xc_5;
                                                                                                                                    yy_6 = xyz{1,fast_obj_start+4}(row_f5(dist_six_loop),2) - yc_5;                                                                                                                                       
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
                                                                                                                                        angle_six_th_test=atan((xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),2)-yc_4)/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),1)-xc_4))*180/pi;                                                           
                                                                                                                                        priority_angle_order6(find_six_loop,1)=abs(abs(angle_final_th)-abs(angle_six_th_test))/priority_degree_rate + distance_six_mat(row_six(find_six_loop));
                                                                                                                                        if( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),2)-yc_4) == 0 & (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),1)- xc_4) )
                                                                                                                                            priority_angle_order6(find_six_loop,1)= priority_zero; 
                                                                                                                                        end
                                                                                                                                        priority_angle_order6(find_six_loop,2)=find_six_loop;%give priority for smaller angle!!
                                                                                                                                        [temp_pr1,ord_pr1] = sort(priority_angle_order6(:,1));
                                                                                                                                        priority_angle_order6 = priority_angle_order6(ord_pr1,:);
                                                                                                                                    end
                                                                                                                                    for find_six_loop2 = 1: row_six_size(1,1)  
                                                                                                                                        if(stop_search2 == 0)
                                                                                                                                            pass_numfil=0; 
                                                                                                                                            diff_vector2=0;
                                                                                                                                            posit_nega2=nan;
                                                                                                                                            xc_6=xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1);
                                                                                                                                            yc_6=xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2);
                                                                                                                                            if(angle_final_th <= 45 & angle_final_th >= -45) %vetorc
                                                                                                                                                vector=1; % x vetorc
                                                                                                                                                diff_vector = xc_5-xc_4;
                                                                                                                                                posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                                            else 
                                                                                                                                                vector=0; % y vetorc
                                                                                                                                                diff_vector = yc_5-yc_4;
                                                                                                                                                posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                                            end
                                                                                                                                            if(vector == 1)       
                                                                                                                                                diff_vector2= xc_6 - xc_5;
                                                                                                                                                posit_nega2=sign(diff_vector2);                                           
                                                                                                                                            elseif(vector == 0)                                            
                                                                                                                                                diff_vector2= yc_6 - yc_5;
                                                                                                                                                posit_nega2=sign(diff_vector2);
                                                                                                                                            end
                                                                                                                                            if( posit_nega == 0 )
                                                                                                                                                xc_4=xc_3;
                                                                                                                                                yc_4=yc_3;       
                                                                                                                                                posit_nega = posit_nega2;        
                                                                                                                                            end
                                                                                                                                            if(posit_nega == posit_nega2)
                                                                                                                                                angle_final_th=atan((yc_5-yc_4)/(xc_5-xc_4))*180/pi;
                                                                                                                                                if( (yc_5-yc_4) == 0 & (xc_5-xc_4) == 0 )
                                                                                                                                                    angle_final_th = 0;
                                                                                                                                                end
                                                                                                                                                dist_final_th=sqrt((xc_5-xc_4)^2+(yc_5-yc_4)^2);
                                                                                                                                                increase_check=1; 
                                                                                                                                                if((distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)))+dist_final_th)*increase_distance_threshold < sqrt( (xc_6-xc_4)^2+(yc_6-yc_4)^2 ))
                                                                                                                                                    increase_check=0;                                                     
                                                                                                                                                else            
                                                                                                                                                    increase_check=1;
                                                                                                                                                end  
                                                                                                                                                if(posit_nega2 == 0)
                                                                                                                                                    pass_numfil=1;
                                                                                                                                                    increase_check=0;
                                                                                                                                                end
                                                                                                                                                if(increase_check == 0)
                                                                                                                                                    angle_six_th_cmp= atan((xc_6-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2))/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)))*180/pi;
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
                                                                                                                                                        fst_mv_oj_coor(5,3) = inside_loop;
                                                                                                                                                        stop_search2=1;
                                                                                                                                                    end
                                                                                                                                                end
                                                                                                                                            end
                                                                                                                                        end
                                                                                                                                    end
                                                                                                                                end                                                                                                                                
                                                                                                                            end
                                                                                                                            %-----prevent duplication----
                                                                                                                            if( prevent_rest_check{1,1} == 0)

                                                                                                                                prevent_rest_check{1,inside_loop}(1,1) = fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,1);
                                                                                                                                prevent_rest_check{1,inside_loop}(1,2) = fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,2); 
                                                                                                                                prevent_rest_check{1,inside_loop}(2,1) = fst_mv_oj_coor(1,1);
                                                                                                                                prevent_rest_check{1,inside_loop}(2,2) = fst_mv_oj_coor(1,2); 
                                                                                                                                prevent_rest_check{1,inside_loop}(3,1) = fst_mv_oj_coor(2,1);
                                                                                                                                prevent_rest_check{1,inside_loop}(3,2) = fst_mv_oj_coor(2,2); 
                                                                                                                                prevent_rest_check{1,inside_loop}(4,1) = fst_mv_oj_coor(3,1);
                                                                                                                                prevent_rest_check{1,inside_loop}(4,2) = fst_mv_oj_coor(3,2);
                                                                                                                                prevent_rest_check{1,inside_loop}(5,1) = fst_mv_oj_coor(4,1);
                                                                                                                                prevent_rest_check{1,inside_loop}(5,2) = fst_mv_oj_coor(4,2);

                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+1,1) = fst_mv_oj_coor(1,1);
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+1,2) = fst_mv_oj_coor(1,2);       
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+1,3) = inside_loop;

                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+2,1) = fst_mv_oj_coor(2,1);
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+2,2) = fst_mv_oj_coor(2,2);         
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+2,3) = inside_loop;

                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+3,1) = fst_mv_oj_coor(3,1);
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+3,2) = fst_mv_oj_coor(3,2);       
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+3,3) = inside_loop;

                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+4,1) = fst_mv_oj_coor(4,1);
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+4,2) = fst_mv_oj_coor(4,2);
                                                                                                                                fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+4,3) = inside_loop;  

                                                                                                                                if(stop_search2==1)
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+5,1) = fst_mv_oj_coor(5,1);
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+5,2) = fst_mv_oj_coor(5,2);
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+5,3) = inside_loop;

                                                                                                                                    prevent_rest_check{1,inside_loop}(6,1)= fst_mv_oj_coor(5,1);           
                                                                                                                                    prevent_rest_check{1,inside_loop}(6,2)= fst_mv_oj_coor(5,2);
                                                                                                                                end
                                                                                                                                stop_search=1; %prevent duplicate trajectories
                                                                                                                            else
                                                                                                                                %------prevent duplication--------
                                                                                                                                maxsize_fstc=size(prevent_rest_check);     
                                                                                                                                maxsize_cmp_fst=size(fst_mv_oj_coor);
                                                                                                                                duplication_number=0;
                                                                                                                                duplecount=0;
                                                                                                                                duplehappen=0;
                                                                                                                                prevent_duplication_loop=0;
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
                                                                                                                                        duplication_number=prevent_duplication_loop;
                                                                                                                                    end
                                                                                                                                end
                                                                                                                                if(duplication_number ~= 0)
                                                                                                                                    if(angle_information(duplication_number,1) < angle_information(inside_loop,1))
                                                                                                                                        duplehappen=1;
                                                                                                                                    end
                                                                                                                                    if( duplehappen == 1)
                                                                                                                                        %duplication happen
                                                                                                                                    else %%%------------------------------check-------------------------------------------
                                                                                                                                        prevent_rest_check{1,inside_loop}(1,1) = fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,1);
                                                                                                                                        prevent_rest_check{1,inside_loop}(1,2) = fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,2); 
                                                                                                                                        prevent_rest_check{1,inside_loop}(2,1) = fst_mv_oj_coor(1,1);
                                                                                                                                        prevent_rest_check{1,inside_loop}(2,2) = fst_mv_oj_coor(1,2); 
                                                                                                                                        prevent_rest_check{1,inside_loop}(3,1) = fst_mv_oj_coor(2,1);
                                                                                                                                        prevent_rest_check{1,inside_loop}(3,2) = fst_mv_oj_coor(2,2); 
                                                                                                                                        prevent_rest_check{1,inside_loop}(4,1) = fst_mv_oj_coor(3,1);
                                                                                                                                        prevent_rest_check{1,inside_loop}(4,2) = fst_mv_oj_coor(3,2);
                                                                                                                                        prevent_rest_check{1,inside_loop}(5,1) = fst_mv_oj_coor(4,1);
                                                                                                                                        prevent_rest_check{1,inside_loop}(5,2) = fst_mv_oj_coor(4,2);

                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+1,1) = fst_mv_oj_coor(1,1);
                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+1,2) = fst_mv_oj_coor(1,2);       
                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+1,3) = inside_loop;

                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+2,1) = fst_mv_oj_coor(2,1);
                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+2,2) = fst_mv_oj_coor(2,2);         
                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+2,3) = inside_loop;

                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+3,1) = fst_mv_oj_coor(3,1);
                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+3,2) = fst_mv_oj_coor(3,2);       
                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+3,3) = inside_loop;

                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+4,1) = fst_mv_oj_coor(4,1);
                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+4,2) = fst_mv_oj_coor(4,2);
                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+4,3) = inside_loop;  

                                                                                                                                        if(stop_search2==1)
                                                                                                                                            prevent_rest_check{1,inside_loop}(6,1)= fst_mv_oj_coor(5,1);           
                                                                                                                                            prevent_rest_check{1,inside_loop}(6,2)= fst_mv_oj_coor(5,2);

                                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+5,1) = fst_mv_oj_coor(5,1);
                                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+5,2) = fst_mv_oj_coor(5,2);
                                                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+5,3) = inside_loop;  
                                                                                                                                        end

                                                                                                                                        %delete duplication before                                                                                                                                    
                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+1,1) = nan;
                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+1,2) = nan;      
                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+1,3) = duplication_number;

                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+2,1) = nan;
                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+2,2) = nan;        
                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+2,3) = duplication_number;

                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+3,1) = nan;
                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+3,2) = nan;       
                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+3,3) = duplication_number;

                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+4,1) = nan;
                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+4,2) = nan;
                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+4,3) = duplication_number; 

                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+5,1) = nan;
                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+5,2) = nan;
                                                                                                                                        fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+5,3) = duplication_number;

                                                                                                                                        duple_num_count=duple_num_count+1;
                                                                                                                                        duple_mat=duplication_number;
                                                                                                                                        %****************************************restore pprevented trajectory***************************
                                                                                                                                        prevent_rest_check{1,duple_mat}=0;
                                                                                                                                        prevent_rest_check{1,duple_mat}(1,1)=nan;
                                                                                                                                        prevent_rest_check{1,duple_mat}(1,2)=nan;                           
                                                                                                                                        distance_fast_mat=0;        
                                                                                                                                        estm_ft_size=size(fast_total_coor{1,duple_mat});
                                                                                                                                        if( estm_ft_size(1,1) >= fast_obj_start )
                                                                                                                                            [row_f1,col_f1,check_f1]=find( xyz{1,fast_obj_start+1}(:,1) < fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,1) + distance_threshold_fast & xyz{1,fast_obj_start+1}(:,1) > fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,1) - distance_threshold_fast);
                                                                                                                                            xc_1=fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,1);
                                                                                                                                            yc_1=fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,2);
                                                                                                                                            xc_back_1=fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point-1,1);
                                                                                                                                            yc_back_1=fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point-1,2);
                                                                                                                                            xc_back_2=fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point-2,1);
                                                                                                                                            yc_back_2=fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point-2,2);
                                                                                                                                        end
                                                                                                                                        if(check_f1 == 1)  
                                                                                                                                            maxsize_fast2_xyz=size(row_f1); %read size
                                                                                                                                            distance_count1=0;
                                                                                                                                            for fast_second_mt=1:maxsize_fast2_xyz(1,1) %for fast sperm
                                                                                                                                                estm_ft_size=size(fast_total_coor{1,duple_mat});
                                                                                                                                                if( estm_ft_size(1,1) >= fast_obj_start )
                                                                                                                                                    xx_1 = xyz{1,fast_obj_start+1}(row_f1(fast_second_mt),1) - xc_1; % comparing x coordinates   
                                                                                                                                                    yy_1 = xyz{1,fast_obj_start+1}(row_f1(fast_second_mt),2) - yc_1; % comparing y coordinates
                                                                                                                                                else     
                                                                                                                                                    xx_1=nan;
                                                                                                                                                    yy_1=nan;
                                                                                                                                                end
                                                                                                                                                distance_count1=distance_count1+1;
                                                                                                                                                distance_fast_mat(distance_count1,1) = sqrt((xx_1^2)+(yy_1^2));
                                                                                                                                                distance_fast_mat(distance_count1,2) = row_f1(fast_second_mt);
                                                                                                                                            end
                                                                                                                                            before_distance=sqrt((xc_1-xc_back_1)^2+(yc_1-yc_back_1)^2);
                                                                                                                                            before_angle=atan((yc_1-yc_back_1)/(xc_1-xc_back_1))*180/pi;
                                                                                                                                            if( (yc_1-yc_back_1) == 0 & (xc_1-xc_back_1) == 0)
                                                                                                                                                before_angle = 0;
                                                                                                                                            end 
                                                                                                                                            if( before_distance <= decide_slow)
                                                                                                                                                [row_dist,col_dist,check_dist]=find(distance_fast_mat(:,1) <= distance_threshold_fast);            
                                                                                                                                            else
                                                                                                                                                [row_dist,col_dist,check_dist]=find(distance_fast_mat(:,1) <= before_distance*threshold_diff_dislocation & distance_fast_mat(:,1) <= distance_threshold_fast);
                                                                                                                                            end
                                                                                                                                            if (check_dist == 1) %first distance threshold check
                                                                                                                                                row_dist_size=size(row_dist); %read size  
                                                                                                                                                angle_fast=0;
                                                                                                                                                dist_thresh_sec=0;
                                                                                                                                                priority_angle_order=0;
                                                                                                                                                for try_loop = 1: row_dist_size(1,1) 
                                                                                                                                                    angle_fast_test = atan((xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(try_loop),2),2)-yc_back_1)/(xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(try_loop),2),1)-xc_back_1))*180/pi; % unit is degree   
                                                                                                                                                    priority_angle_order(try_loop,1)=abs(abs(angle_fast_test)-abs(before_angle))/priority_degree_rate + distance_fast_mat(row_dist(try_loop));                  
                                                                                                                                                    if( (xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(try_loop),2),2)-yc_back_1) == 0 & (xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(try_loop),2),1)-xc_back_1) == 0 )
                                                                                                                                                        priority_angle_order(try_loop,1) = priority_zero; 
                                                                                                                                                    end
                                                                                                                                                    priority_angle_order(try_loop,2)=try_loop;%give priority for smaller angle!!
                                                                                                                                                    [temp_pr1,ord_pr1] = sort(priority_angle_order(:,1));
                                                                                                                                                    priority_angle_order = priority_angle_order(ord_pr1,:);
                                                                                                                                                end
                                                                                                                                                for try_loop2 = 1: row_dist_size(1,1)
                                                                                                                                                    if(stop_search == 0) 
                                                                                                                                                        pass_numfil=0;
                                                                                                                                                        xc_2=xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1);
                                                                                                                                                        yc_2=xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2);
                                                                                                                                                        if(before_angle <= 45 & before_angle >= -45) %vetorc
                                                                                                                                                            vector=1; % x vetorc
                                                                                                                                                            diff_vector = xc_1 - xc_back_1 ;
                                                                                                                                                            posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                                                        else 
                                                                                                                                                            vector=0; % y vetorc
                                                                                                                                                            diff_vector = yc_1 - yc_back_1 ;
                                                                                                                                                            posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                                                        end
                                                                                                                                                        diff_vector2=0;
                                                                                                                                                        posit_nega2=nan;
                                                                                                                                                        if(vector == 1)         
                                                                                                                                                            diff_vector2= xc_2-xc_1;    
                                                                                                                                                            posit_nega2=sign(diff_vector2);
                                                                                                                                                        elseif (vector == 0)            
                                                                                                                                                            diff_vector2= yc_2-yc_1;
                                                                                                                                                            posit_nega2=sign(diff_vector2);          
                                                                                                                                                        end                      
                                                                                                                                                        if( posit_nega == 0 )                     
                                                                                                                                                            posit_nega = posit_nega2;
                                                                                                                                                            xc_back_1=xc_back_2;
                                                                                                                                                            yc_back_1=yc_back_2;
                                                                                                                                                        end
                                                                                                                                                        if(posit_nega == posit_nega2)                            
                                                                                                                                                            before_angle=atan((yc_1-yc_back_1)/(xc_1-xc_back_1))*180/pi;
                                                                                                                                                            if( (yc_1-yc_back_1) == 0 & (xc_1-xc_back_1) == 0)
                                                                                                                                                                before_angle = 0;
                                                                                                                                                            end        
                                                                                                                                                            before_distance=sqrt((xc_1-xc_back_1)^2+(yc_1-yc_back_1)^2);
                                                                                                                                                            angle_fast_test=atan((yc_2-yc_1)/(xc_2-xc_1))*180/pi;
                                                                                                                                                            if( (yc_2-yc_1) == 0 & (xc_2-xc_1) == 0)
                                                                                                                                                                angle_fast_test = 0;
                                                                                                                                                            end
                                                                                                                                                            if( before_distance > decide_slow_angle)           
                                                                                                                                                                if(abs(before_angle) > tangent_control | abs(angle_fast_test) > tangent_control)              
                                                                                                                                                                    if(abs(before_angle)+control_vecter_angle > abs(atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi) & abs(before_angle)-control_vecter_angle < abs(atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi)) % unit is degree
                                                                                                                                                                        pass_numfil=1;
                                                                                                                                                                    end
                                                                                                                                                                else
                                                                                                                                                                    if(before_angle+control_vecter_angle > atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi & before_angle-control_vecter_angle < atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi) % unit is degree
                                                                                                                                                                    pass_numfil=1;
                                                                                                                                                                    end
                                                                                                                                                                end    
                                                                                                                                                            else
                                                                                                                                                                if(abs(before_angle) > tangent_control | abs(angle_fast_test) > tangent_control)              
                                                                                                                                                                    if(abs(before_angle)+control_vecter_angle*multilple_short_angle > abs(atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi) & abs(before_angle)-control_vecter_angle*multilple_short_angle < abs(atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi)) % unit is degree
                                                                                                                                                                        pass_numfil=1;
                                                                                                                                                                    end
                                                                                                                                                                else
                                                                                                                                                                    if(before_angle+control_vecter_angle*multilple_short_angle > atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi & before_angle-control_vecter_angle*multilple_short_angle < atan((yc_2-yc_back_1)/(xc_2-xc_back_1))*180/pi) % unit is degree
                                                                                                                                                                    pass_numfil=1;
                                                                                                                                                                    end
                                                                                                                                                                end
                                                                                                                                                            end                                   
                                                                                                                                                        end
                                                                                                                                                        increase_check=1;
                                                                                                                                                        if((distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)))+before_distance)*increase_distance_threshold < sqrt( (xc_2-xc_back_1)^2+(yc_2-yc_back_1)^2 ))
                                                                                                                                                            increase_check=0; 
                                                                                                                                                        else   
                                                                                                                                                            increase_check=1;
                                                                                                                                                        end        
                                                                                                                                                        if(posit_nega2 == 0)
                                                                                                                                                            pass_numfil=1;
                                                                                                                                                            increase_check=0;
                                                                                                                                                        end
                                                                                                                                                        if(increase_check == 0)       
                                                                                                                                                            if(pass_numfil == 1) 
                                                                                                                                                                if(fast_obj_start > 1+change_start_point)
                                                                                                                                                                    angle_fast = atan((yc_2-yc_1)/(xc_2-xc_1))*180/pi;
                                                                                                                                                                    if( (yc_2-yc_1) == 0 & (xc_2-xc_1) == 0)
                                                                                                                                                                        angle_fast = 0;
                                                                                                                                                                    end
                                                                                                                                                                    dist_thresh_sec = distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)));                                
                                                                                                                                                                end   
                                                                                                                                                                distance_repeat_rest=0;
                                                                                                                                                                [row_f2,col_f2,check_f2]=find( xyz{1,fast_obj_start+2}(:,1) < xc_2 + distance_threshold_fast & xyz{1,fast_obj_start+2}(:,1) > xc_2 - distance_threshold_fast);       

                                                                                                                                                                if(check_f2 == 1)
                                                                                                                                                                    maxsize_fast_xyz2=size(row_f2); %read size
                                                                                                                                                                    distance_count2=0;
                                                                                                                                                                    for dist_cmp_loop = 1: maxsize_fast_xyz2(1,1)
                                                                                                                                                                        xx_2 = xyz{1,fast_obj_start+2}(row_f2(dist_cmp_loop),1) - xc_2;                              
                                                                                                                                                                        yy_2 = xyz{1,fast_obj_start+2}(row_f2(dist_cmp_loop),2) - yc_2;                 
                                                                                                                                                                        distance_count2=distance_count2+1;
                                                                                                                                                                        distance_repeat_rest(distance_count2,1) = sqrt((xx_2^2)+(yy_2^2));
                                                                                                                                                                        distance_repeat_rest(distance_count2,2) = row_f2(dist_cmp_loop);
                                                                                                                                                                    end 
                                                                                                                                                                    if(dist_thresh_sec <= decide_slow)
                                                                                                                                                                        [row_secd,col_secd,check_secd]=find(distance_repeat_rest(:,1) < distance_threshold_fast/2);
                                                                                                                                                                    else
                                                                                                                                                                        [row_secd,col_secd,check_secd]=find(distance_repeat_rest(:,1) < dist_thresh_sec*threshold_diff_dislocation & distance_repeat_rest(:,1) < distance_threshold_fast);
                                                                                                                                                                    end
                                                                                                                                                                    if (check_secd == 1)
                                                                                                                                                                        row_secd_size=size(row_secd); %read size  
                                                                                                                                                                        angle_secd_fast=0;
                                                                                                                                                                        dist_thresh_third=0;
                                                                                                                                                                        priority_angle_order2=0;
                                                                                                                                                                        for find_secd_loop = 1: row_secd_size(1,1)
                                                                                                                                                                            angle_secd_fast_test = atan((xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),2)-yc_1)/(xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),1)-xc_1))*180/pi;
                                                                                                                                                                            priority_angle_order2(find_secd_loop,1)=abs(abs(angle_secd_fast_test)-abs(angle_fast))/priority_degree_rate + distance_repeat_rest(row_secd(find_secd_loop));
                                                                                                                                                                            if( (xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),2)-yc_1) == 0 & (xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(find_secd_loop),2),1)-xc_1) == 0 )
                                                                                                                                                                                priority_angle_order2(find_secd_loop,1) = priority_zero; 
                                                                                                                                                                            end
                                                                                                                                                                            priority_angle_order2(find_secd_loop,2)=find_secd_loop;%give priority for smaller angle!!
                                                                                                                                                                            [temp_pr1,ord_pr1] = sort(priority_angle_order2(:,1));
                                                                                                                                                                            priority_angle_order2 = priority_angle_order2(ord_pr1,:);
                                                                                                                                                                        end
                                                                                                                                                                        for find_secd_loop2 = 1: row_secd_size(1,1)
                                                                                                                                                                            if(stop_search == 0)
                                                                                                                                                                                xc_3=xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1);
                                                                                                                                                                                yc_3=xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2);
                                                                                                                                                                                pass_numfil=0;
                                                                                                                                                                                diff_vector2=0;
                                                                                                                                                                                posit_nega2=nan;
                                                                                                                                                                                if( fast_obj_start > 1+change_start_point)
                                                                                                                                                                                     if(angle_fast <= 45 & angle_fast >= -45) %vetorc
                                                                                                                                                                                        vector=1; % x vetorc
                                                                                                                                                                                        diff_vector = xc_2-xc_1;
                                                                                                                                                                                        posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                                                                                    else 
                                                                                                                                                                                        vector=0; % y vetorc
                                                                                                                                                                                        diff_vector = yc_2-yc_1;
                                                                                                                                                                                        posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                                                                                     end
                                                                                                                                                                                end
                                                                                                                                                                                if(vector == 1)
                                                                                                                                                                                    diff_vector2= xc_3-xc_2;
                                                                                                                                                                                    posit_nega2=sign(diff_vector2);
                                                                                                                                                                                elseif (vector == 0)
                                                                                                                                                                                    diff_vector2= yc_3-yc_2;                            
                                                                                                                                                                                    posit_nega2=sign(diff_vector2);
                                                                                                                                                                                end
                                                                                                                                                                                if( posit_nega == 0 )
                                                                                                                                                                                    if( fast_obj_start > 1+change_start_point)
                                                                                                                                                                                        xc_1=xc_back_1;
                                                                                                                                                                                        yc_1=yc_back_1;
                                                                                                                                                                                    end
                                                                                                                                                                                    posit_nega = posit_nega2;        
                                                                                                                                                                                end
                                                                                                                                                                                if(posit_nega == posit_nega2)
                                                                                                                                                                                    angle_fast=atan((yc_2-yc_1)/(xc_2-xc_1))*180/pi;
                                                                                                                                                                                    if( (yc_2-yc_1) == 0 & (xc_2-xc_1) == 0 )
                                                                                                                                                                                        angle_fast=0;
                                                                                                                                                                                    end
                                                                                                                                                                                    dist_thresh_sec=sqrt((xc_2-xc_1)^2+(yc_2-yc_1)^2);
                                                                                                                                                                                    increase_check=1;
                                                                                                                                                                                    if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xc_3-xc_1)^2+(yc_3-yc_1)^2 ))
                                                                                                                                                                                        increase_check=0;
                                                                                                                                                                                    else
                                                                                                                                                                                        increase_check=1;
                                                                                                                                                                                    end
                                                                                                                                                                                    if(posit_nega2 == 0)
                                                                                                                                                                                        pass_numfil=1;
                                                                                                                                                                                        increase_check=0;
                                                                                                                                                                                    end
                                                                                                                                                                                    if(increase_check == 0)
                                                                                                                                                                                        angle_secd_fast_cmp=atan((yc_3-yc_2)/(xc_3-xc_2))*180/pi; 
                                                                                                                                                                                        if( (yc_3-yc_2) == 0 & (xc_3-xc_2) == 0 )
                                                                                                                                                                                            angle_secd_fast_cmp=0;
                                                                                                                                                                                        end                                                
                                                                                                                                                                                        if(dist_thresh_sec > decide_slow_angle)
                                                                                                                                                                                            if(abs(angle_fast)>tangent_control | abs(angle_secd_fast_cmp) > tangent_control)                                                     
                                                                                                                                                                                                if(abs(angle_fast)+control_vecter_angle > abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) & abs(angle_fast)-control_vecter_angle < abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi)) % unit is degree                       
                                                                                                                                                                                                    pass_numfil=1;
                                                                                                                                                                                                end
                                                                                                                                                                                            else
                                                                                                                                                                                                if(angle_fast+control_vecter_angle > atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi & angle_fast-control_vecter_angle < atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) % unit is degree  
                                                                                                                                                                                                    pass_numfil=1;
                                                                                                                                                                                                end
                                                                                                                                                                                            end
                                                                                                                                                                                        else
                                                                                                                                                                                            if(abs(angle_fast)>tangent_control | abs(angle_secd_fast_cmp) > tangent_control)                                                     
                                                                                                                                                                                                if(abs(angle_fast)+control_vecter_angle*multilple_short_angle > abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) & abs(angle_fast)-control_vecter_angle*multilple_short_angle < abs(atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi)) % unit is degree                       
                                                                                                                                                                                                    pass_numfil=1;
                                                                                                                                                                                                end
                                                                                                                                                                                            else
                                                                                                                                                                                                if(angle_fast+control_vecter_angle*multilple_short_angle > atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi & angle_fast-control_vecter_angle*multilple_short_angle < atan((yc_3-yc_1 )/(xc_3-xc_1))*180/pi) % unit is degree  
                                                                                                                                                                                                    pass_numfil=1;
                                                                                                                                                                                                end
                                                                                                                                                                                            end
                                                                                                                                                                                        end
                                                                                                                                                                                        if(pass_numfil == 1) 
                                                                                                                                                                                            angle_secd_fast = atan((yc_3-yc_2)/(xc_3-xc_2))*180/pi;   
                                                                                                                                                                                            if( (yc_3-yc_2) == 0 & (xc_3-xc_2) == 0 )
                                                                                                                                                                                                angle_secd_fast=0;
                                                                                                                                                                                            end
                                                                                                                                                                                            dist_thresh_third = distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)));                           
                                                                                                                                                                                            % finish third
                                                                                                                                                                                            distance_thrd_mat=0;
                                                                                                                                                                                            [row_f3,col_f3,check_f3]=find( xyz{1,fast_obj_start+3}(:,1) < xc_3 + distance_threshold_fast & xyz{1,fast_obj_start+3}(:,1) > xc_3 - distance_threshold_fast);       

                                                                                                                                                                                            if(check_f3 == 1)
                                                                                                                                                                                                maxsize_fast_xyz3=size(row_f3);
                                                                                                                                                                                                distance_count3=0;
                                                                                                                                                                                                for dist_thrd_loop = 1: maxsize_fast_xyz3(1,1)                                   
                                                                                                                                                                                                    xx_3 = xyz{1,fast_obj_start+3}(row_f3(dist_thrd_loop),1) - xc_3;
                                                                                                                                                                                                    yy_3 = xyz{1,fast_obj_start+3}(row_f3(dist_thrd_loop),2) - yc_3;
                                                                                                                                                                                                    distance_count3=distance_count3+1;
                                                                                                                                                                                                    distance_thrd_mat(distance_count3,1) = sqrt((xx_3^2)+(yy_3^2));
                                                                                                                                                                                                    distance_thrd_mat(distance_count3,2) = row_f3(dist_thrd_loop);
                                                                                                                                                                                                end
                                                                                                                                                                                                if(dist_thresh_third <= decide_slow)
                                                                                                                                                                                                    [row_third,col_third,check_third]=find(distance_thrd_mat(:,1) < distance_threshold_fast);
                                                                                                                                                                                                else
                                                                                                                                                                                                    [row_third,col_third,check_third]=find(distance_thrd_mat(:,1) < dist_thresh_third*threshold_diff_dislocation & distance_thrd_mat(:,1) <= distance_threshold_fast );
                                                                                                                                                                                                end
                                                                                                                                                                                                if (check_third == 1)
                                                                                                                                                                                                    row_third_size2=size(row_third);
                                                                                                                                                                                                    angle_rest_th=0;
                                                                                                                                                                                                    dist_rest_th=0;                         
                                                                                                                                                                                                    priority_angle_order3=0;
                                                                                                                                                                                                    for find_rest_loop = 1: row_third_size2(1,1)
                                                                                                                                                                                                        angle_rest_th_test=atan((xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),2)-yc_2)/(xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),1)-xc_2))*180/pi;                                                           
                                                                                                                                                                                                        priority_angle_order3(find_rest_loop,1)=abs(abs(angle_rest_th_test)-abs(angle_secd_fast))/priority_degree_rate + distance_thrd_mat(row_third(find_rest_loop));
                                                                                                                                                                                                        if( (xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),2)-yc_2) == 0 & (xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(find_rest_loop),2),1)-xc_2) ==0 )
                                                                                                                                                                                                            priority_angle_order3(find_rest_loop,1) = priority_zero;
                                                                                                                                                                                                        end
                                                                                                                                                                                                        priority_angle_order3(find_rest_loop,2)=find_rest_loop;%give priority for smaller angle!!
                                                                                                                                                                                                        [temp_pr1,ord_pr1] = sort(priority_angle_order3(:,1));
                                                                                                                                                                                                        priority_angle_order3 = priority_angle_order3(ord_pr1,:);
                                                                                                                                                                                                    end
                                                                                                                                                                                                    for find_rest_loop2 = 1: row_third_size2(1,1) 
                                                                                                                                                                                                        if(stop_search == 0)
                                                                                                                                                                                                            pass_numfil=0; 
                                                                                                                                                                                                            diff_vector2=0;
                                                                                                                                                                                                            posit_nega2=nan;
                                                                                                                                                                                                            xc_4=xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1);
                                                                                                                                                                                                            yc_4=xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2);
                                                                                                                                                                                                            if(angle_secd_fast <= 45 & angle_secd_fast >= -45) %vetorc
                                                                                                                                                                                                                vector=1; % x vetorc
                                                                                                                                                                                                                diff_vector = xc_3-xc_2;
                                                                                                                                                                                                                posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                                                                                                            else 
                                                                                                                                                                                                                vector=0; % y vetorc
                                                                                                                                                                                                                diff_vector = yc_3-yc_2;
                                                                                                                                                                                                                posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                                                                                                            end
                                                                                                                                                                                                            if(vector == 1)       
                                                                                                                                                                                                                diff_vector2= xc_4 - xc_3;
                                                                                                                                                                                                                posit_nega2=sign(diff_vector2);                                           
                                                                                                                                                                                                            elseif(vector == 0)                                            
                                                                                                                                                                                                                diff_vector2= yc_4 - yc_3;                                            
                                                                                                                                                                                                                posit_nega2=sign(diff_vector2);
                                                                                                                                                                                                            end
                                                                                                                                                                                                            if( posit_nega == 0 )
                                                                                                                                                                                                                xc_2=xc_1;
                                                                                                                                                                                                                yc_2=yc_1;
                                                                                                                                                                                                                posit_nega = posit_nega2;        
                                                                                                                                                                                                            end
                                                                                                                                                                                                            if(posit_nega == posit_nega2)
                                                                                                                                                                                                                angle_secd_fast=atan((yc_3-yc_2)/(xc_3-xc_2))*180/pi;
                                                                                                                                                                                                                if( (yc_3-yc_2) == 0 & (xc_3-xc_2) == 0 )
                                                                                                                                                                                                                    angle_secd_fast = 0;
                                                                                                                                                                                                                end
                                                                                                                                                                                                                dist_thresh_third=sqrt((xc_3-xc_2)^2+(yc_3-yc_2)^2);
                                                                                                                                                                                                                increase_check=1;
                                                                                                                                                                                                                if((distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)))+dist_thresh_third)*increase_distance_threshold < sqrt( (xc_4-xc_2)^2+(yc_4-yc_2)^2 ))
                                                                                                                                                                                                                    increase_check=0;                                                     
                                                                                                                                                                                                                else            
                                                                                                                                                                                                                    increase_check=1;
                                                                                                                                                                                                                end 
                                                                                                                                                                                                                if(posit_nega2 == 0)
                                                                                                                                                                                                                    pass_numfil=1;
                                                                                                                                                                                                                    increase_check=0;
                                                                                                                                                                                                                end
                                                                                                                                                                                                                if(increase_check == 0)  
                                                                                                                                                                                                                    angle_rest_th_cmp=atan((yc_4-yc_3)/(xc_4-xc_3))*180/pi;
                                                                                                                                                                                                                    if( (yc_4-yc_3) == 0 & (xc_4-xc_3)==0)
                                                                                                                                                                                                                        angle_rest_th_cmp=0;
                                                                                                                                                                                                                    end
                                                                                                                                                                                                                    if( dist_thresh_third > decide_slow_angle)
                                                                                                                                                                                                                        if(abs(angle_secd_fast)>tangent_control | abs(angle_rest_th_cmp)>tangent_control)    
                                                                                                                                                                                                                            if(abs(angle_secd_fast)+control_vecter_angle > abs(atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi) & abs(angle_secd_fast)-control_vecter_angle < abs(atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi))  
                                                                                                                                                                                                                                pass_numfil=1;
                                                                                                                                                                                                                            end
                                                                                                                                                                                                                        else
                                                                                                                                                                                                                            if(angle_secd_fast+control_vecter_angle > atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi & angle_secd_fast-control_vecter_angle < atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi)  
                                                                                                                                                                                                                                pass_numfil=1;
                                                                                                                                                                                                                            end
                                                                                                                                                                                                                        end
                                                                                                                                                                                                                    else
                                                                                                                                                                                                                        if(abs(angle_secd_fast)>tangent_control | abs(angle_rest_th_cmp)>tangent_control)    
                                                                                                                                                                                                                            if(abs(angle_secd_fast)+control_vecter_angle*multilple_short_angle > abs(atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi) & abs(angle_secd_fast)-control_vecter_angle*multilple_short_angle < abs(atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi))  
                                                                                                                                                                                                                                pass_numfil=1;
                                                                                                                                                                                                                            end
                                                                                                                                                                                                                        else
                                                                                                                                                                                                                            if(angle_secd_fast+control_vecter_angle*multilple_short_angle > atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi & angle_secd_fast-control_vecter_angle*multilple_short_angle < atan((yc_4-yc_2)/(xc_4-xc_2))*180/pi)  
                                                                                                                                                                                                                                pass_numfil=1;
                                                                                                                                                                                                                            end
                                                                                                                                                                                                                        end 
                                                                                                                                                                                                                    end
                                                                                                                                                                                                                    if( pass_numfil == 1) 
                                                                                                                                                                                                                        angle_rest_th=atan((yc_4-yc_3)/(xc_4-xc_3))*180/pi;                 
                                                                                                                                                                                                                        if ( (yc_4-yc_3) == 0 & (xc_4-xc_3) == 0 )
                                                                                                                                                                                                                            angle_rest_th = 0;
                                                                                                                                                                                                                        end
                                                                                                                                                                                                                        dist_rest_th=distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)));
                                                                                                                                                                                                                        [row_f4,col_f4,check_f4]=find( xyz{1,fast_obj_start+4}(:,1) < xc_4 + distance_threshold_fast & xyz{1,fast_obj_start+4}(:,1) > xc_4 - distance_threshold_fast);       

                                                                                                                                                                                                                        if(check_f4 == 1)
                                                                                                                                                                                                                            distance_four_mat=0;
                                                                                                                                                                                                                            maxsize_fast_xyz4=size(row_f4);
                                                                                                                                                                                                                            distance_count4=0;
                                                                                                                                                                                                                            for dist_four_loop = 1: maxsize_fast_xyz4(1,1)
                                                                                                                                                                                                                                xx_4 = xyz{1,fast_obj_start+4}(row_f4(dist_four_loop),1) - xc_4;    
                                                                                                                                                                                                                                yy_4 = xyz{1,fast_obj_start+4}(row_f4(dist_four_loop),2) - yc_4;                                
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
                                                                                                                                                                                                                                    angle_final_th_test=atan((xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),2)-yc_3)/(xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),1)-xc_3))*180/pi;
                                                                                                                                                                                                                                    priority_angle_order4(find_four_loop,1)=abs(abs(angle_rest_th)-abs(angle_final_th_test))/priority_degree_rate + distance_four_mat(row_four(find_four_loop));
                                                                                                                                                                                                                                    if( (xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),2)-yc_3) == 0 & (xyz{1,fast_obj_start+4}(distance_four_mat(row_four(find_four_loop),2),1)-xc_3) == 0)
                                                                                                                                                                                                                                        priority_angle_order4(find_four_loop,1)=priority_zero;
                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                    priority_angle_order4(find_four_loop,2)=find_four_loop;%give priority for smaller angle!!
                                                                                                                                                                                                                                    [temp_pr1,ord_pr1] = sort(priority_angle_order4(:,1));
                                                                                                                                                                                                                                    priority_angle_order4 = priority_angle_order4(ord_pr1,:);
                                                                                                                                                                                                                                end
                                                                                                                                                                                                                                for find_four_loop2 = 1: row_four_size(1,1)                                                          
                                                                                                                                                                                                                                    if(stop_search == 0)
                                                                                                                                                                                                                                        xc_5=xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1);
                                                                                                                                                                                                                                        yc_5=xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2);
                                                                                                                                                                                                                                        pass_numfil=0;
                                                                                                                                                                                                                                        diff_vector2=0;
                                                                                                                                                                                                                                        posit_nega2=nan;
                                                                                                                                                                                                                                        if(angle_rest_th <= 45 & angle_rest_th >= -45) %vetorc
                                                                                                                                                                                                                                            vector=1; % x vetorc
                                                                                                                                                                                                                                            diff_vector = xc_4-xc_3;
                                                                                                                                                                                                                                            posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                                                                                                                                        else 
                                                                                                                                                                                                                                            vector=0; % y vetorc
                                                                                                                                                                                                                                            diff_vector = yc_4-yc_3;
                                                                                                                                                                                                                                            posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                        if(vector == 1)       
                                                                                                                                                                                                                                            diff_vector2= xc_5 - xc_4;      
                                                                                                                                                                                                                                            posit_nega2=sign(diff_vector2);        
                                                                                                                                                                                                                                        elseif(vector == 0)                                            
                                                                                                                                                                                                                                            diff_vector2= yc_5 - yc_4;                                            
                                                                                                                                                                                                                                            posit_nega2=sign(diff_vector2);
                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                        if( posit_nega == 0 )
                                                                                                                                                                                                                                            xc_3=xc_2;
                                                                                                                                                                                                                                            yc_3=yc_2;
                                                                                                                                                                                                                                            posit_nega = posit_nega2;        
                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                        if(posit_nega == posit_nega2)
                                                                                                                                                                                                                                            angle_rest_th=atan((yc_4-yc_3)/(xc_4-xc_3))*180/pi;
                                                                                                                                                                                                                                            if( (yc_4-yc_3) == 0 & (xc_4-xc_3) == 0 )
                                                                                                                                                                                                                                                angle_rest_th = 0;
                                                                                                                                                                                                                                            end
                                                                                                                                                                                                                                            dist_rest_th= sqrt((xc_4-xc_3)^2+(yc_4-yc_3)^2);
                                                                                                                                                                                                                                            increase_check=1;
                                                                                                                                                                                                                                            if((distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)))+dist_rest_th)*increase_distance_threshold < sqrt((xc_5-xc_3)^2+(yc_5-yc_3)^2))       
                                                                                                                                                                                                                                                increase_check=0;    
                                                                                                                                                                                                                                            else
                                                                                                                                                                                                                                                increase_check=1;      
                                                                                                                                                                                                                                            end
                                                                                                                                                                                                                                            if(posit_nega2 == 0)
                                                                                                                                                                                                                                                pass_numfil=1;
                                                                                                                                                                                                                                                increase_check=0;
                                                                                                                                                                                                                                            end
                                                                                                                                                                                                                                            if(increase_check == 0)                                                                               
                                                                                                                                                                                                                                                angle_final_th_cmp=atan((yc_5-yc_4)/(xc_5-xc_4))*180/pi;
                                                                                                                                                                                                                                                if( (yc_5-yc_4) == 0 & (xc_5-xc_4) ==0 )
                                                                                                                                                                                                                                                    angle_final_th_cmp=0;
                                                                                                                                                                                                                                                end
                                                                                                                                                                                                                                                if( dist_rest_th > decide_slow_angle)
                                                                                                                                                                                                                                                    if(abs(angle_rest_th)>tangent_control | abs(angle_final_th_cmp) >tangent_control)
                                                                                                                                                                                                                                                        if(abs(angle_rest_th)+control_vecter_angle > abs(atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi) & abs(angle_rest_th)-control_vecter_angle < abs(atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi))
                                                                                                                                                                                                                                                            pass_numfil=1;
                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                    else
                                                                                                                                                                                                                                                        if(angle_rest_th+control_vecter_angle > atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi & angle_rest_th-control_vecter_angle < atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi)
                                                                                                                                                                                                                                                            pass_numfil=1;
                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                else  
                                                                                                                                                                                                                                                     if(abs(angle_rest_th)>tangent_control | abs(angle_final_th_cmp) >tangent_control)
                                                                                                                                                                                                                                                        if(abs(angle_rest_th)+control_vecter_angle*multilple_short_angle > abs(atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi) & abs(angle_rest_th)-control_vecter_angle*multilple_short_angle  < abs(atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi))
                                                                                                                                                                                                                                                            pass_numfil=1;
                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                    else
                                                                                                                                                                                                                                                        if(angle_rest_th+control_vecter_angle*multilple_short_angle  > atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi & angle_rest_th-control_vecter_angle*multilple_short_angle  < atan((yc_5-yc_3)/(xc_5-xc_3))*180/pi)
                                                                                                                                                                                                                                                            pass_numfil=1;
                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                end
                                                                                                                                                                                                                                                if( pass_numfil == 1)
                                                                                                                                                                                                                                                    angle_final_th=atan((yc_5-yc_4)/(xc_5-xc_4))*180/pi;                                    
                                                                                                                                                                                                                                                    if( (yc_5-yc_4) == 0 & (xc_5-xc_4)==0)
                                                                                                                                                                                                                                                        angle_final_th=0;
                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                    dist_final_th=distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)));

                                                                                                                                                                                                                                                    fst_mv_oj_coor=0;

                                                                                                                                                                                                                                                    fst_mv_oj_coor(1,1) = xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1);
                                                                                                                                                                                                                                                    fst_mv_oj_coor(1,2) = xyz{1,fast_obj_start+1}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2);       
                                                                                                                                                                                                                                                    fst_mv_oj_coor(1,3) = duple_mat;

                                                                                                                                                                                                                                                    fst_mv_oj_coor(2,1) = xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1);
                                                                                                                                                                                                                                                    fst_mv_oj_coor(2,2) = xyz{1,fast_obj_start+2}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2);          
                                                                                                                                                                                                                                                    fst_mv_oj_coor(2,3) = duple_mat;

                                                                                                                                                                                                                                                    fst_mv_oj_coor(3,1) = xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1);
                                                                                                                                                                                                                                                    fst_mv_oj_coor(3,2) = xyz{1,fast_obj_start+3}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2); 
                                                                                                                                                                                                                                                    fst_mv_oj_coor(3,3) = duple_mat;

                                                                                                                                                                                                                                                    fst_mv_oj_coor(4,1) = xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1);
                                                                                                                                                                                                                                                    fst_mv_oj_coor(4,2) = xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2);
                                                                                                                                                                                                                                                    fst_mv_oj_coor(4,3) = duple_mat; 

                                                                                                                                                                                                                                                    if( abs(angle_final_th) > tangent_control)
                                                                                                                                                                                                                                                        angle_information(object_number_count,1)=(abs(abs(atan((fst_mv_oj_coor(2,2)-fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,2))/(fst_mv_oj_coor(2,1)-fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,1)))*180/pi)-abs(atan((fst_mv_oj_coor(1,2)-fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,2))/(fst_mv_oj_coor(1,1)-fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,1)))*180/pi))+abs(abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(1,1)))*180/pi)-abs(atan((fst_mv_oj_coor(2,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(2,1)-fst_mv_oj_coor(1,1)))*180/pi))+abs(abs(atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(2,1)))*180/pi)-abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(2,1)))*180/pi)))/priority_degree_rate + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) );
                                                                                                                                                                                                                                                    else
                                                                                                                                                                                                                                                        angle_information(object_number_count,1)=(abs(atan((fst_mv_oj_coor(2,2)-fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,2))/(fst_mv_oj_coor(2,1)-fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,1)))*180/pi-atan((fst_mv_oj_coor(1,2)-fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,2))/(fst_mv_oj_coor(1,1)-fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,1)))*180/pi)+abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(1,1)))*180/pi-atan((fst_mv_oj_coor(2,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(2,1)-fst_mv_oj_coor(1,1)))*180/pi)+abs(atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(2,1)))*180/pi-atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(2,1)))*180/pi))/priority_degree_rate + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) );
                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                    distance_six_mat=0;
                                                                                                                                                                                                                                                    [row_f5,col_f5,check_f5]=find( xyz{1,fast_obj_start+4}(:,1) < xc_5 + distance_threshold_fast & xyz{1,fast_obj_start+4}(:,1) > xc_5 - distance_threshold_fast);       

                                                                                                                                                                                                                                                    if(check_f5 == 1)
                                                                                                                                                                                                                                                        stop_search2=0;
                                                                                                                                                                                                                                                        maxsize_fast_xyz6=size(row_f5);
                                                                                                                                                                                                                                                        distance_count6=0;
                                                                                                                                                                                                                                                        for dist_six_loop = 1: maxsize_fast_xyz6(1,1)                                   
                                                                                                                                                                                                                                                            xx_6 = xyz{1,fast_obj_start+4}(row_f5(dist_six_loop),1) - xc_5;
                                                                                                                                                                                                                                                            yy_6 = xyz{1,fast_obj_start+4}(row_f5(dist_six_loop),2) - yc_5;                                                                                                                                       
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
                                                                                                                                                                                                                                                                angle_six_th_test=atan((xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),2)-yc_4)/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),1)-xc_4))*180/pi;                                                           
                                                                                                                                                                                                                                                                priority_angle_order6(find_six_loop,1)=abs(abs(angle_final_th)-abs(angle_six_th_test))/priority_degree_rate + distance_six_mat(row_six(find_six_loop));
                                                                                                                                                                                                                                                                if( (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),2)-yc_4) == 0 & (xyz{1,fast_obj_start+4}(distance_six_mat(row_six(find_six_loop),2),1)- xc_4) )
                                                                                                                                                                                                                                                                    priority_angle_order6(find_six_loop,1)= priority_zero; 
                                                                                                                                                                                                                                                                end
                                                                                                                                                                                                                                                                priority_angle_order6(find_six_loop,2)=find_six_loop;%give priority for smaller angle!!
                                                                                                                                                                                                                                                                [temp_pr1,ord_pr1] = sort(priority_angle_order6(:,1));
                                                                                                                                                                                                                                                                priority_angle_order6 = priority_angle_order6(ord_pr1,:);
                                                                                                                                                                                                                                                            end
                                                                                                                                                                                                                                                            for find_six_loop2 = 1: row_six_size(1,1)  
                                                                                                                                                                                                                                                                if(stop_search2 == 0)
                                                                                                                                                                                                                                                                    pass_numfil=0; 
                                                                                                                                                                                                                                                                    diff_vector2=0;
                                                                                                                                                                                                                                                                    posit_nega2=nan;
                                                                                                                                                                                                                                                                    xc_6=xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1);
                                                                                                                                                                                                                                                                    yc_6=xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2);
                                                                                                                                                                                                                                                                    if(angle_final_th <= 45 & angle_final_th >= -45) %vetorc
                                                                                                                                                                                                                                                                        vector=1; % x vetorc
                                                                                                                                                                                                                                                                        diff_vector = xc_5-xc_4;
                                                                                                                                                                                                                                                                        posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                                                                                                                                                                                                                                                                    else 
                                                                                                                                                                                                                                                                        vector=0; % y vetorc
                                                                                                                                                                                                                                                                        diff_vector = yc_5-yc_4;
                                                                                                                                                                                                                                                                        posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                                    if(vector == 1)       
                                                                                                                                                                                                                                                                        diff_vector2= xc_6 - xc_5;
                                                                                                                                                                                                                                                                        posit_nega2=sign(diff_vector2);                                           
                                                                                                                                                                                                                                                                    elseif(vector == 0)                                            
                                                                                                                                                                                                                                                                        diff_vector2= yc_6 - yc_5;
                                                                                                                                                                                                                                                                        posit_nega2=sign(diff_vector2);
                                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                                    if( posit_nega == 0 )
                                                                                                                                                                                                                                                                        xc_4=xc_3;
                                                                                                                                                                                                                                                                        yc_4=yc_3;       
                                                                                                                                                                                                                                                                        posit_nega = posit_nega2;        
                                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                                    if(posit_nega == posit_nega2)
                                                                                                                                                                                                                                                                        angle_final_th=atan((yc_5-yc_4)/(xc_5-xc_4))*180/pi;
                                                                                                                                                                                                                                                                        if( (yc_5-yc_4) == 0 & (xc_5-xc_4) == 0 )
                                                                                                                                                                                                                                                                            angle_final_th = 0;
                                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                                        dist_final_th=sqrt((xc_5-xc_4)^2+(yc_5-yc_4)^2);
                                                                                                                                                                                                                                                                        increase_check=1; 
                                                                                                                                                                                                                                                                        if((distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)))+dist_final_th)*increase_distance_threshold < sqrt( (xc_6-xc_4)^2+(yc_6-yc_4)^2 ))
                                                                                                                                                                                                                                                                            increase_check=0;                                                     
                                                                                                                                                                                                                                                                        else            
                                                                                                                                                                                                                                                                            increase_check=1;
                                                                                                                                                                                                                                                                        end  
                                                                                                                                                                                                                                                                        if(posit_nega2 == 0)
                                                                                                                                                                                                                                                                            pass_numfil=1;
                                                                                                                                                                                                                                                                            increase_check=0;
                                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                                        if(increase_check == 0)
                                                                                                                                                                                                                                                                            angle_six_th_cmp= atan((xc_6-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2))/(xyz{1,fast_obj_start+4}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1)-xyz{1,fast_obj_start+4}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1)))*180/pi;
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
                                                                                                                                                                                                                                                                                fst_mv_oj_coor(5,3) = duple_mat;
                                                                                                                                                                                                                                                                                stop_search2=1;
                                                                                                                                                                                                                                                                            end
                                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                                end
                                                                                                                                                                                                                                                            end
                                                                                                                                                                                                                                                        end                                                                                                                                
                                                                                                                                                                                                                                                    end
                                                                                                                                                                                                                                                    %-----prevent duplication----
                                                                                                                                                                                                                                                    if( prevent_rest_check{1,1} == 0)

                                                                                                                                                                                                                                                        prevent_rest_check{1,duple_mat}(1,1) = fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,1);
                                                                                                                                                                                                                                                        prevent_rest_check{1,duple_mat}(1,2) = fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,2); 
                                                                                                                                                                                                                                                        prevent_rest_check{1,duple_mat}(2,1) = fst_mv_oj_coor(1,1);
                                                                                                                                                                                                                                                        prevent_rest_check{1,duple_mat}(2,2) = fst_mv_oj_coor(1,2); 
                                                                                                                                                                                                                                                        prevent_rest_check{1,duple_mat}(3,1) = fst_mv_oj_coor(2,1);
                                                                                                                                                                                                                                                        prevent_rest_check{1,duple_mat}(3,2) = fst_mv_oj_coor(2,2); 
                                                                                                                                                                                                                                                        prevent_rest_check{1,duple_mat}(4,1) = fst_mv_oj_coor(3,1);
                                                                                                                                                                                                                                                        prevent_rest_check{1,duple_mat}(4,2) = fst_mv_oj_coor(3,2);
                                                                                                                                                                                                                                                        prevent_rest_check{1,duple_mat}(5,1) = fst_mv_oj_coor(4,1);
                                                                                                                                                                                                                                                        prevent_rest_check{1,duple_mat}(5,2) = fst_mv_oj_coor(4,2);

                                                                                                                                                                                                                                                        fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+1,1) = fst_mv_oj_coor(1,1);
                                                                                                                                                                                                                                                        fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+1,2) = fst_mv_oj_coor(1,2);       
                                                                                                                                                                                                                                                        fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+1,3) = duple_mat;

                                                                                                                                                                                                                                                        fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+2,1) = fst_mv_oj_coor(2,1);
                                                                                                                                                                                                                                                        fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+2,2) = fst_mv_oj_coor(2,2);         
                                                                                                                                                                                                                                                        fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+2,3) = duple_mat;

                                                                                                                                                                                                                                                        fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+3,1) = fst_mv_oj_coor(3,1);
                                                                                                                                                                                                                                                        fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+3,2) = fst_mv_oj_coor(3,2);       
                                                                                                                                                                                                                                                        fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+3,3) = duple_mat;

                                                                                                                                                                                                                                                        fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+4,1) = fst_mv_oj_coor(4,1);
                                                                                                                                                                                                                                                        fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+4,2) = fst_mv_oj_coor(4,2);
                                                                                                                                                                                                                                                        fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+4,3) = duple_mat;  

                                                                                                                                                                                                                                                        if(stop_search2==1)
                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+5,1) = fst_mv_oj_coor(5,1);
                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+5,2) = fst_mv_oj_coor(5,2);
                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+5,3) = duple_mat;

                                                                                                                                                                                                                                                            prevent_rest_check{1,duple_mat}(6,1)= fst_mv_oj_coor(5,1);           
                                                                                                                                                                                                                                                            prevent_rest_check{1,duple_mat}(6,2)= fst_mv_oj_coor(5,2);
                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                        stop_search=1; %prevent duplicate trajectories
                                                                                                                                                                                                                                                    else
                                                                                                                                                                                                                                                        %------prevent duplication--------
                                                                                                                                                                                                                                                        maxsize_fstc=size(prevent_rest_check);     
                                                                                                                                                                                                                                                        maxsize_cmp_fst=size(fst_mv_oj_coor);
                                                                                                                                                                                                                                                        duplication_number=0;
                                                                                                                                                                                                                                                        duplecount=0;
                                                                                                                                                                                                                                                        duplehappen=0;
                                                                                                                                                                                                                                                        prevent_duplication_loop=0;
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
                                                                                                                                                                                                                                                                duplication_number=prevent_duplication_loop;
                                                                                                                                                                                                                                                            end
                                                                                                                                                                                                                                                        end
                                                                                                                                                                                                                                                        if(duplication_number ~= 0)
                                                                                                                                                                                                                                                            if(angle_information(duplication_number,1) < angle_information(duple_mat,1))
                                                                                                                                                                                                                                                                duplehappen=1;
                                                                                                                                                                                                                                                            end
                                                                                                                                                                                                                                                            if( duplehappen == 1)
                                                                                                                                                                                                                                                                %duplication happen
                                                                                                                                                                                                                                                            else %%%------------------------------check-------------------------------------------
                                                                                                                                                                                                                                                                prevent_rest_check{1,duple_mat}(1,1) = fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,1);
                                                                                                                                                                                                                                                                prevent_rest_check{1,duple_mat}(1,2) = fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,2); 
                                                                                                                                                                                                                                                                prevent_rest_check{1,duple_mat}(2,1) = fst_mv_oj_coor(1,1);
                                                                                                                                                                                                                                                                prevent_rest_check{1,duple_mat}(2,2) = fst_mv_oj_coor(1,2); 
                                                                                                                                                                                                                                                                prevent_rest_check{1,duple_mat}(3,1) = fst_mv_oj_coor(2,1);
                                                                                                                                                                                                                                                                prevent_rest_check{1,duple_mat}(3,2) = fst_mv_oj_coor(2,2); 
                                                                                                                                                                                                                                                                prevent_rest_check{1,duple_mat}(4,1) = fst_mv_oj_coor(3,1);
                                                                                                                                                                                                                                                                prevent_rest_check{1,duple_mat}(4,2) = fst_mv_oj_coor(3,2);
                                                                                                                                                                                                                                                                prevent_rest_check{1,duple_mat}(5,1) = fst_mv_oj_coor(4,1);
                                                                                                                                                                                                                                                                prevent_rest_check{1,duple_mat}(5,2) = fst_mv_oj_coor(4,2);

                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+1,1) = fst_mv_oj_coor(1,1);
                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+1,2) = fst_mv_oj_coor(1,2);       
                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+1,3) = duple_mat;

                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+2,1) = fst_mv_oj_coor(2,1);
                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+2,2) = fst_mv_oj_coor(2,2);         
                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+2,3) = duple_mat;

                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+3,1) = fst_mv_oj_coor(3,1);
                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+3,2) = fst_mv_oj_coor(3,2);       
                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+3,3) = duple_mat;

                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+4,1) = fst_mv_oj_coor(4,1);
                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+4,2) = fst_mv_oj_coor(4,2);
                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+4,3) = duple_mat;  

                                                                                                                                                                                                                                                                if(stop_search2==1)
                                                                                                                                                                                                                                                                    prevent_rest_check{1,duple_mat}(6,1)= fst_mv_oj_coor(5,1);           
                                                                                                                                                                                                                                                                    prevent_rest_check{1,duple_mat}(6,2)= fst_mv_oj_coor(5,2);

                                                                                                                                                                                                                                                                    fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+5,1) = fst_mv_oj_coor(5,1);
                                                                                                                                                                                                                                                                    fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+5,2) = fst_mv_oj_coor(5,2);
                                                                                                                                                                                                                                                                    fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+5,3) = duple_mat;  
                                                                                                                                                                                                                                                                end

                                                                                                                                                                                                                                                                %delete duplication before                                                                                                                                    
                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+1,1) = nan;
                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+1,2) = nan;      
                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+1,3) = duplication_number;

                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+2,1) = nan;
                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+2,2) = nan;        
                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+2,3) = duplication_number;

                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+3,1) = nan;
                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+3,2) = nan;       
                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+3,3) = duplication_number;

                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+4,1) = nan;
                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+4,2) = nan;
                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+4,3) = duplication_number; 

                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+5,1) = nan;
                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+5,2) = nan;
                                                                                                                                                                                                                                                                fast_total_coor{1,duplication_number}(fast_obj_start-change_start_point+5,3) = duplication_number;

                                                                                                                                                                                                                                                                duple_mat=duplication_number;
                                                                                                                                                                                                                                                                stop_search=1; %prevent duplicate trajectories
                                                                                                                                                                                                                                                            end
                                                                                                                                                                                                                                                        else
                                                                                                                                                                                                                                                            prevent_rest_check{1,duple_mat}(1,1) = fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,1);
                                                                                                                                                                                                                                                            prevent_rest_check{1,duple_mat}(1,2) = fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point,2); 
                                                                                                                                                                                                                                                            prevent_rest_check{1,duple_mat}(2,1) = fst_mv_oj_coor(1,1);
                                                                                                                                                                                                                                                            prevent_rest_check{1,duple_mat}(2,2) = fst_mv_oj_coor(1,2); 
                                                                                                                                                                                                                                                            prevent_rest_check{1,duple_mat}(3,1) = fst_mv_oj_coor(2,1);
                                                                                                                                                                                                                                                            prevent_rest_check{1,duple_mat}(3,2) = fst_mv_oj_coor(2,2); 
                                                                                                                                                                                                                                                            prevent_rest_check{1,duple_mat}(4,1) = fst_mv_oj_coor(3,1);
                                                                                                                                                                                                                                                            prevent_rest_check{1,duple_mat}(4,2) = fst_mv_oj_coor(3,2);
                                                                                                                                                                                                                                                            prevent_rest_check{1,duple_mat}(5,1) = fst_mv_oj_coor(4,1);
                                                                                                                                                                                                                                                            prevent_rest_check{1,duple_mat}(5,2) = fst_mv_oj_coor(4,2);

                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+1,1) = fst_mv_oj_coor(1,1);
                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+1,2) = fst_mv_oj_coor(1,2);       
                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+1,3) = duple_mat;

                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+2,1) = fst_mv_oj_coor(2,1);
                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+2,2) = fst_mv_oj_coor(2,2);         
                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+2,3) = duple_mat;

                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+3,1) = fst_mv_oj_coor(3,1);
                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+3,2) = fst_mv_oj_coor(3,2);       
                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+3,3) = duple_mat;

                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+4,1) = fst_mv_oj_coor(4,1);
                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+4,2) = fst_mv_oj_coor(4,2);
                                                                                                                                                                                                                                                            fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+4,3) = duple_mat;  

                                                                                                                                                                                                                                                            if(stop_search2==1)
                                                                                                                                                                                                                                                                prevent_rest_check{1,duple_mat}(6,1)= fst_mv_oj_coor(5,1);           
                                                                                                                                                                                                                                                                prevent_rest_check{1,duple_mat}(6,2)= fst_mv_oj_coor(5,2);

                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+5,1) = fst_mv_oj_coor(5,1);
                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+5,2) = fst_mv_oj_coor(5,2);
                                                                                                                                                                                                                                                                fast_total_coor{1,duple_mat}(fast_obj_start-change_start_point+5,3) = duple_mat;
                                                                                                                                                                                                                                                            end
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
                                                                                                                                        if( stop_search == 0 )
                                                                                                                                            prevent_rest_check{1,duple_mat}(1,1) = nan;
                                                                                                                                            prevent_rest_check{1,duple_mat}(1,2) = nan;
                                                                                                                                        end
                                                                                                                                        stop_search=1; %prevent duplicate trajectories
                                                                                                                                        %---restore finish---
                                                                                                                                    end
                                                                                                                                else
                                                                                                                                    prevent_rest_check{1,inside_loop}(1,1) = fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,1);
                                                                                                                                    prevent_rest_check{1,inside_loop}(1,2) = fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point,2); 
                                                                                                                                    prevent_rest_check{1,inside_loop}(2,1) = fst_mv_oj_coor(1,1);
                                                                                                                                    prevent_rest_check{1,inside_loop}(2,2) = fst_mv_oj_coor(1,2); 
                                                                                                                                    prevent_rest_check{1,inside_loop}(3,1) = fst_mv_oj_coor(2,1);
                                                                                                                                    prevent_rest_check{1,inside_loop}(3,2) = fst_mv_oj_coor(2,2); 
                                                                                                                                    prevent_rest_check{1,inside_loop}(4,1) = fst_mv_oj_coor(3,1);
                                                                                                                                    prevent_rest_check{1,inside_loop}(4,2) = fst_mv_oj_coor(3,2);
                                                                                                                                    prevent_rest_check{1,inside_loop}(5,1) = fst_mv_oj_coor(4,1);
                                                                                                                                    prevent_rest_check{1,inside_loop}(5,2) = fst_mv_oj_coor(4,2);

                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+1,1) = fst_mv_oj_coor(1,1);
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+1,2) = fst_mv_oj_coor(1,2);       
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+1,3) = inside_loop;

                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+2,1) = fst_mv_oj_coor(2,1);
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+2,2) = fst_mv_oj_coor(2,2);         
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+2,3) = inside_loop;

                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+3,1) = fst_mv_oj_coor(3,1);
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+3,2) = fst_mv_oj_coor(3,2);       
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+3,3) = inside_loop;

                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+4,1) = fst_mv_oj_coor(4,1);
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+4,2) = fst_mv_oj_coor(4,2);
                                                                                                                                    fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+4,3) = inside_loop;  

                                                                                                                                    if(stop_search2==1)
                                                                                                                                        prevent_rest_check{1,inside_loop}(6,1)= fst_mv_oj_coor(5,1);           
                                                                                                                                        prevent_rest_check{1,inside_loop}(6,2)= fst_mv_oj_coor(5,2);

                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+5,1) = fst_mv_oj_coor(5,1);
                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+5,2) = fst_mv_oj_coor(5,2);
                                                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start-change_start_point+5,3) = inside_loop;
                                                                                                                                    end
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
            if( fast_obj_start > 1+change_start_point)
                if( stop_search == 0 )
                    prevent_rest_check{1,inside_loop}(1,1) = nan;
                    prevent_rest_check{1,inside_loop}(1,2) = nan;
                end
            end
        end
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


max_frame_size=fix((iter-2)/2);

skip_frame_nb=2*max_frame_size+1

new_number_count=object_number_count;

for change_slow_point= 0:2
    for rest_start_searching= 1:max_frame_size-1

        fast_obj_start=2*rest_start_searching-1+change_slow_point;

        if(fast_obj_start == 1+change_slow_point)
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
            maxsize_fast2_xyz=size(xyz{1,fast_obj_start+1}); %read size    

            for fast_second_mt=1:maxsize_fast2_xyz(1,1) %for fast sperm

                if(fast_obj_start == 1+change_slow_point)

                    xx_1 = xyz{1,fast_obj_start+1}(fast_second_mt,1) - xyz_11(inside_loop,1); % comparing x coordinates   
                    yy_1 = xyz{1,fast_obj_start+1}(fast_second_mt,2) - xyz_11(inside_loop,2); % comparing y coordinates

                    distance_fast_mat(fast_second_mt,1) = sqrt((xx_1^2)+(yy_1^2));
                else
                    estm_ft_size=size(fast_total_coor{1,inside_loop});
                    if( estm_ft_size(1,1) > rest_start_searching )             
                        xx_1 = xyz{1,fast_obj_start+1}(fast_second_mt,1) - fast_total_coor{1,inside_loop}(rest_start_searching,1); % comparing x coordinates   
                        yy_1 = xyz{1,fast_obj_start+1}(fast_second_mt,2) - fast_total_coor{1,inside_loop}(rest_start_searching,2); % comparing y coordinates                  
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

                if(fast_obj_start == 1+change_slow_point)  
                    diff_vector=0;
                    posit_nega=nan;
                    vector=nan;
                    for try_loop = 1: row_dist_size(1,1)
                        priority_angle_order(try_loop,1)=try_loop;
                        priority_angle_order(try_loop,2)=try_loop;%give priority for smaller angle!!
                    end
                else 
                    for try_loop = 1: row_dist_size(1,1)
                        angle_fast_test = atan((xyz{1,fast_obj_start+1}(row_dist(try_loop),2)-fast_total_coor{1,inside_loop}(rest_start_searching,2))/(xyz{1,fast_obj_start+1}(row_dist(try_loop),1)-fast_total_coor{1,inside_loop}(rest_start_searching,1)))*180/pi; % unit is degree   
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

                    if(fast_obj_start == 1+change_slow_point)
                        %decide vector andgle
                        %each magnitude
                        angle_fast = atan((xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2)-xyz_11(inside_loop,2))/(xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)-xyz_11(inside_loop,1)))*180/pi; % unit is degree
                        dist_thresh_sec = distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)));

                        if(angle_fast <= 45 & angle_fast >= -45) %vetorc
                            vector=1; % x vetorc
                            diff_vector = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1) - xyz_11(inside_loop,1);
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
                            diff_vector = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1) - fast_total_coor{1,inside_loop}(rest_start_searching,1);
                            posit_nega=sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                        else 
                            vector=0; % y vetorc
                            diff_vector = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2) - fast_total_coor{1,inside_loop}(rest_start_searching,2); 
                            posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                        end

                        diff_vector2=0;
                        posit_nega2=nan;  

                        if(vector == 1)         
                            diff_vector2= xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(rest_start_searching,1);    
                            posit_nega2=sign(diff_vector2);
                        elseif (vector == 0)            
                            diff_vector2= xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(rest_start_searching,2);
                            posit_nega2=sign(diff_vector2);          
                        end

                        angle_fast_test=atan((xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(rest_start_searching,2))/(xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(rest_start_searching,1)))*180/pi;

                        if(abs(cont_angle_th(inside_loop,1))+control_vecter_angle*multilple_short_angle > abs(angle_fast_test) & abs(cont_angle_th(inside_loop,1))-control_vecter_angle*multilple_short_angle < abs(angle_fast_test)) % unit is degree
                            pass_numfil=1;
                        end                   
                    end

                    if(posit_nega == posit_nega2) 
                        pass_value=1;   
                        if(fast_obj_start == 1+change_slow_point)                                
                            pass_value=0;   
                        else 
                            if((distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)))+cont_rest_th(inside_loop,1))*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(rest_start_searching-1,1))^2+(xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(rest_start_searching-1,2))^2 ))    
                                pass_value=0; 
                            else 
                                pass_value=1; 
                            end 
                        end
                        if(pass_value == 0) 
                            if(pass_numfil == 1) 
                                if(fast_obj_start > 1+change_slow_point)      
                                    angle_fast = atan((xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(rest_start_searching,2))/(xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(rest_start_searching,1)))*180/pi;   
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
                                [row_secd,col_secd,check_secd]=find(distance_repeat_rest <= distance_threshold_skip);

                                if (check_secd == 1)
                                    row_secd_size=size(row_secd); %read size   
                                    angle_secd_fast=0; 
                                    dist_thresh_third=0;  
                                    priority_angle_order2=0;

                                    for find_secd_loop = 1: row_secd_size(1,1) 
                                        angle_secd_fast_test = atan((xyz{1,fast_obj_start+2}(row_secd(find_secd_loop),2)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2))/(xyz{1,fast_obj_start+2}(row_secd(find_secd_loop),1)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)))*180/pi;

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
                                                diff_vector2= xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1);
                                                posit_nega2=sign(diff_vector2);
                                            elseif (vector == 0)
                                                diff_vector2= xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2);                            
                                                posit_nega2=sign(diff_vector2);
                                            end

                                            if(posit_nega == posit_nega2)
                                                pass_value=1;

                                                if(fast_obj_start == 1+change_slow_point)
                                                    if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz_11(inside_loop,1))^2+(xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz_11(inside_loop,2))^2 ))
                                                        pass_value=0;
                                                    else
                                                        pass_value=1;
                                                    end

                                                else
                                                    if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-fast_total_coor{1,inside_loop}(rest_start_searching,1))^2+(xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-fast_total_coor{1,inside_loop}(rest_start_searching,2))^2 ))
                                                        pass_value=0;
                                                    else
                                                        pass_value=1;
                                                    end
                                                end

                                                if(pass_value == 0)                                       
                                                    angle_secd_fast_cmp=atan((xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2))/(xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)))*180/pi; 

                                                    if(abs(angle_fast)+control_vecter_angle*multilple_short_angle > abs(angle_secd_fast_cmp) & abs(angle_fast)-control_vecter_angle*multilple_short_angle < abs(angle_secd_fast_cmp)) % unit is degree                       
                                                        pass_numfil=1;
                                                    end

                                                    if(pass_numfil == 1) 
                                                        angle_secd_fast = atan((xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2))/(xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)))*180/pi;   
                                                        dist_thresh_third = distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)));                           
                                                        % finish third

                                                        if(fast_obj_start == 1+change_slow_point)

                                                            new_number_count = new_number_count+1;
                                                            fst_mv_oj_coor=0;

                                                            fst_mv_oj_coor(1,1) = xyz_11(inside_loop,1);                              
                                                            fst_mv_oj_coor(1,2) = xyz_11(inside_loop,2);          
                                                            fst_mv_oj_coor(1,3) = new_number_count;

                                                            fst_mv_oj_coor(2,1) = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1);
                                                            fst_mv_oj_coor(2,2) = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2);           
                                                            fst_mv_oj_coor(2,3) = new_number_count;

                                                            fst_mv_oj_coor(3,1) = xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1); 
                                                            fst_mv_oj_coor(3,2) = xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2);                        
                                                            fst_mv_oj_coor(3,3) = new_number_count;

                                                            angle_information(new_number_count,1)=abs(abs(angle_fast-angle_secd_fast));

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
                                                                new_number_count = new_number_count-1; %duplication happen                         
                                                            else                           
                                                                fast_total_coor(1,new_number_count)={fst_mv_oj_coor};   
                                                                cont_rest_th(new_number_count,1)=dist_thresh_sec;                        
                                                                cont_angle_th(new_number_count,1)=angle_fast;                    
                                                                stop_search=1; %prevent duplicate trajectories                              
                                                            end

                                                        else

                                                            fast_total_coor{1,inside_loop}(rest_start_searching+1,1) = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1);
                                                            fast_total_coor{1,inside_loop}(rest_start_searching+1,2) = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2);                                                      
                                                            fast_total_coor{1,inside_loop}(rest_start_searching+1,3) = inside_loop;

                                                            fast_total_coor{1,inside_loop}(rest_start_searching+2,1) = xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1);
                                                            fast_total_coor{1,inside_loop}(rest_start_searching+2,2) = xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2);                                 
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


max_frame_size=fix((iter-2)/4);

skip_frame_nb=4*max_frame_size+1

new_number_count2=new_number_count;

for rest_start_searching= 1:max_frame_size

    fast_obj_start=4*rest_start_searching-3;

    if(fast_obj_start == 1)
        maxsize_fast_xyz=size(xyz_11); %read size
        change_inside_loop = maxsize_fast_xyz(1,1);
        start_number_skip=1;
    else
        change_inside_loop = new_number_count2; 
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
                        diff_vector = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2) - xyz_11(inside_loop,2); 
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
    
                                                        new_number_count2 = new_number_count2+1;
                                                        fst_mv_oj_coor=0;
                                         
                                                        fst_mv_oj_coor(1,1) = xyz_11(inside_loop,1);                              
                                                        fst_mv_oj_coor(1,2) = xyz_11(inside_loop,2);          
                                                        fst_mv_oj_coor(1,3) = new_number_count2;
                                       
                                                        fst_mv_oj_coor(2,1) = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1);
                                                        fst_mv_oj_coor(2,2) = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2);           
                                                        fst_mv_oj_coor(2,3) = new_number_count2;
     
                                                        fst_mv_oj_coor(3,1) = xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),1); 
                                                        fst_mv_oj_coor(3,2) = xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),2);                        
                                                        fst_mv_oj_coor(3,3) = new_number_count2;
                    
                                                        angle_information(new_number_count2,1)=abs(abs(angle_fast-angle_secd_fast));
  
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
                                                            new_number_count2 = new_number_count2-1; %duplication happen                         
                                                        else                           
                                                            fast_total_coor(1,new_number_count2)={fst_mv_oj_coor};   
                                                            cont_rest_th(new_number_count2,1)=dist_thresh_sec;                        
                                                            cont_angle_th(new_number_count2,1)=angle_fast;                    
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

% figure
% imshow(abs(recon_stack(:,:,1)),[])
% hold on
% plot(xyz{1,1}(:,1),xyz{1,1}(:,2),'g.')
% set(gca, 'XLim', [800, 900], 'YLim', [100,200])
% axis on
% 
% figure
% imshow(abs(recon_stack(:,:,1)),[])
% hold on
% plot(sum3_xyz(:,1),sum3_xyz(:,2),'g.')
% set(gca, 'XLim', [650, 680], 'YLim', [20,70])
% axis on
% 
% figure, imshow(recon_stack_diff(:,:,1),[])
% hold on
% plot(xyz{1,1}(:,1),xyz{1,1}(:,2),'g.')
% set(gca, 'XLim', [2420, 2520], 'YLim', [550,650])
% axis on
% 
% set(gca, 'XLim', [900,960], 'YLim', [0,60])
% axis on

% set(gca, 'XLim', [550,900], 'YLim', [350,600])
 
hold on
plot([1182; 2462], [170; 170], '--w',  [1182; 2462], [1130; 1130], '--w', 'LineWidth', 2)
plot([1182; 1182], [170; 1130], '--w',  [2462; 2462], [170; 1130], '--w', 'LineWidth', 2)

hold on
plot([2000; 2454.545455], [1800; 1800], '-w', 'LineWidth', 4)
text(2227.272727,1835, '1 mm', 'HorizontalAlignment','center','color','w','Fontsize',25)