%-------main trajectories-----------------------------
function [fast_total_coor,object_number_count]=change_start_point_frame(xyz,fast_total_coor,distance_threshold_fast,priority_degree_rate,threshold_diff_dislocation,priority_zero,decide_slow,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,increase_distance_threshold,rest_start_searching,max_frame_size,frame_number,iter,object_number_count,xyz_1)

xc_back_1=nan;
yc_back_1=nan;

for change_start_frame=1:2

    %------------------------delete finding--------------------------------------------------
    prevent_duplication_strt2=0;
    prevent_duplication_strt3=0;
    maxsize_fst_coor= size(fast_total_coor);

    for delete_strt_loop= 1:maxsize_fst_coor(1,2)
        prevent_duplication_strt2(delete_strt_loop,1)=fast_total_coor{1,delete_strt_loop}(change_start_frame,1);
        prevent_duplication_strt2(delete_strt_loop,2)=fast_total_coor{1,delete_strt_loop}(change_start_frame,2);
    end

    for delete_strt_loop= 1:maxsize_fst_coor(1,2)
        prevent_duplication_strt3(delete_strt_loop,1)=fast_total_coor{1,delete_strt_loop}(change_start_frame+1,1);
        prevent_duplication_strt3(delete_strt_loop,2)=fast_total_coor{1,delete_strt_loop}(change_start_frame+1,2);
    end

    for delete_strt_loop= 1:maxsize_fst_coor(1,2)
        prevent_duplication_strt4(delete_strt_loop,1)=fast_total_coor{1,delete_strt_loop}(change_start_frame+2,1);
        prevent_duplication_strt4(delete_strt_loop,2)=fast_total_coor{1,delete_strt_loop}(change_start_frame+2,2);
    end

    xyz_new_check=0;
    comp_prevent_strt=0;

    xyz_new_check=xyz{1,change_start_frame+1};
    comp_prevent_strt=[prevent_duplication_strt2; prevent_duplication_strt3];

    size_cmp_prevent_1=size(comp_prevent_strt);

    for loop_cmp_del=1:size_cmp_prevent_1(1,1)

        [row_del_x]=find(xyz_new_check(:,1) == comp_prevent_strt(loop_cmp_del,1));
        [row_del_y]=find(xyz_new_check(:,2) == comp_prevent_strt(loop_cmp_del,2));

        rr_ss=size(row_del_x);
        for check_loop_dlt=1:rr_ss
            [row_del_f,col_del_f,check_del_f]=find(row_del_y == row_del_x(check_loop_dlt));
            if(check_del_f ~= 0)
                xyz_new_check(row_del_x(check_loop_dlt),:)=nan;
            end
        end
    end

    %-----processed here-----thinking! you can do it
    maxsize_fast_xyz=size(xyz_new_check);

    for inside_loop=1:maxsize_fast_xyz(1,1)  
        fast_obj_start=(4*rest_start_searching-3)+change_start_frame;
        origin_size_fst=fast_obj_start;
        stop_search=0;
        increase_number=1;
        jump_number=0;

        if(stop_search==0)
            ans_xyz= change_first_start_point(xyz_new_check,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,priority_zero);
            if(isnan(ans_xyz))
                increase_number=increase_number+1;
                jump_number=jump_number+1;
                ans_xyz= change_first_start_point(xyz_new_check,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,priority_zero);
            end
            if(isnan(ans_xyz))
                increase_number=increase_number-1;
                jump_number=jump_number-1;
                ans_xyz= change_first_start_point(xyz_new_check,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,priority_zero);
            end

            if(~isnan(ans_xyz))
                try_ans=size(ans_xyz);
                keep_in_number_1=increase_number;
                jump_remember_1=jump_number;
                for try_ans1_loop=1:try_ans(1,1)
                    if(stop_search==0)
                        xc_1=ans_xyz(try_ans1_loop,1);
                        yc_1=ans_xyz(try_ans1_loop,2);
                        xc_2=ans_xyz(try_ans1_loop,3);
                        yc_2=ans_xyz(try_ans1_loop,4);
                        posit_nega=ans_xyz(try_ans1_loop,5);
                        vector_check=ans_xyz(try_ans1_loop,6);
                        angle_fast=ans_xyz(try_ans1_loop,7);
                        dist_thresh_sec=ans_xyz(try_ans1_loop,8);

                        jump_number=jump_remember_1;
                        %finish second
                        increase_number=keep_in_number_1+1;
                        ans_xyz2= second_point(xc_back_1,yc_back_1,xc_1,yc_1,xc_2,yc_2,posit_nega,angle_fast,vector_check,dist_thresh_sec,decide_slow,threshold_diff_dislocation,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,increase_distance_threshold,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,priority_zero);
                        if( (rest_start_searching ~= max_frame_size) & (fast_obj_start+increase_number+2 < frame_number) ) 
                            if(isnan(ans_xyz2))
                                increase_number=increase_number+1;
                                jump_number=jump_number+1;
                                ans_xyz2= second_point(xc_back_1,yc_back_1,xc_1,yc_1,xc_2,yc_2,posit_nega,angle_fast,vector_check,dist_thresh_sec,decide_slow,threshold_diff_dislocation,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,increase_distance_threshold,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,priority_zero);
                            end
                        end
                        if(isnan(ans_xyz2))
                            increase_number=increase_number-1;
                            jump_number=jump_number-1;
                            ans_xyz2= second_point(xc_back_1,yc_back_1,xc_1,yc_1,xc_2,yc_2,posit_nega,angle_fast,vector_check,dist_thresh_sec,decide_slow,threshold_diff_dislocation,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,increase_distance_threshold,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,priority_zero);
                        end
                        if(~isnan(ans_xyz2))
                            try_ans2=size(ans_xyz2);
                            keep_in_number_2=increase_number;
                            jump_remember_2=jump_number;
                            for try_ans2_loop=1:try_ans2(1,1)
                                if(stop_search==0)
                                    xc_3=ans_xyz2(try_ans2_loop,1);
                                    yc_3=ans_xyz2(try_ans2_loop,2);
                                    angle_secd_fast=ans_xyz2(try_ans2_loop,3);
                                    dist_thresh_third=ans_xyz2(try_ans2_loop,4);

                                    jump_number=jump_remember_2;
                                    %finish third
                                    increase_number=keep_in_number_2+1;
                                    ans_xyz3= third_point(xc_1,yc_1,xc_2,yc_2,xc_3,yc_3,angle_secd_fast,dist_thresh_third,decide_slow,threshold_diff_dislocation,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,increase_distance_threshold,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,priority_zero);
                                    if( (rest_start_searching ~= max_frame_size) & (fast_obj_start+increase_number+1 < frame_number) )
                                        if(isnan(ans_xyz3))
                                            increase_number=increase_number+1;
                                            jump_number=jump_number+1;
                                            ans_xyz3= third_point(xc_1,yc_1,xc_2,yc_2,xc_3,yc_3,angle_secd_fast,dist_thresh_third,decide_slow,threshold_diff_dislocation,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,increase_distance_threshold,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,priority_zero);
                                        end
                                    end
                                    if(isnan(ans_xyz3))
                                        increase_number=increase_number-1;
                                        jump_number=jump_number-1;
                                        ans_xyz3= third_point(xc_1,yc_1,xc_2,yc_2,xc_3,yc_3,angle_secd_fast,dist_thresh_third,decide_slow,threshold_diff_dislocation,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,increase_distance_threshold,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,priority_zero);
                                    end
                                    if(~isnan(ans_xyz3))
                                        try_ans3=size(ans_xyz3);
                                        keep_in_number_3=increase_number;
                                        jump_remember_3=jump_number;
                                        for try_ans3_loop=1:try_ans3(1,1)
                                            if(stop_search==0)
                                                xc_4=ans_xyz3(try_ans3_loop,1);
                                                yc_4=ans_xyz3(try_ans3_loop,2);
                                                angle_rest_th=ans_xyz3(try_ans3_loop,3);
                                                dist_rest_th=ans_xyz3(try_ans3_loop,4);

                                                jump_number=jump_remember_3;
                                                %finish fourth
                                                increase_number=keep_in_number_3+1;
                                                ans_xyz4= fourth_point(xc_2,yc_2,xc_3,yc_3,xc_4,yc_4,angle_rest_th,dist_rest_th,decide_slow,threshold_diff_dislocation,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,increase_distance_threshold,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,priority_zero);
                                                if( (rest_start_searching ~= max_frame_size) & (fast_obj_start+increase_number < frame_number) )
                                                    if(isnan(ans_xyz4))
                                                        increase_number=increase_number+1;
                                                        jump_number=jump_number+1;
                                                        ans_xyz4= fourth_point(xc_2,yc_2,xc_3,yc_3,xc_4,yc_4,angle_rest_th,dist_rest_th,decide_slow,threshold_diff_dislocation,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,increase_distance_threshold,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,priority_zero);
                                                    end
                                                end
                                                if(isnan(ans_xyz4))
                                                    increase_number=increase_number-1;
                                                    jump_number=jump_number-1;
                                                    ans_xyz4= fourth_point(xc_2,yc_2,xc_3,yc_3,xc_4,yc_4,angle_rest_th,dist_rest_th,decide_slow,threshold_diff_dislocation,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,increase_distance_threshold,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,priority_zero);
                                                end
                                                if(~isnan(ans_xyz4))
                                                    try_ans4=size(ans_xyz4);
                                                    keep_in_number_4=increase_number;
                                                    jump_remember_4=jump_number;
                                                    for try_ans4_loop=1:try_ans4(1,1)
                                                        if(stop_search==0)
                                                            xc_5=ans_xyz4(try_ans4_loop,1);
                                                            yc_5=ans_xyz4(try_ans4_loop,2);
                                                            angle_final_th=ans_xyz4(try_ans4_loop,3);
                                                            dist_final_th=ans_xyz4(try_ans4_loop,4);

                                                            jump_number=jump_remember_4+change_start_frame;
                                                            %finish fifth

                                                            object_number_count = object_number_count+1;
                                                            fst_mv_oj_coor=0;

                                                            fst_mv_oj_coor(1,1) = xc_1;                              
                                                            fst_mv_oj_coor(1,2) = yc_1;          
                                                            fst_mv_oj_coor(1,3) = object_number_count;
                                                            fst_mv_oj_coor(1,4) = jump_number;

                                                            fst_mv_oj_coor(2,1) = xc_2;
                                                            fst_mv_oj_coor(2,2) = yc_2;       
                                                            fst_mv_oj_coor(2,3) = object_number_count;
                                                            fst_mv_oj_coor(2,4) = jump_number;

                                                            fst_mv_oj_coor(3,1) = xc_3;
                                                            fst_mv_oj_coor(3,2) = yc_3;          
                                                            fst_mv_oj_coor(3,3) = object_number_count;
                                                            fst_mv_oj_coor(3,4) = jump_number;

                                                            fst_mv_oj_coor(4,1) = xc_4;
                                                            fst_mv_oj_coor(4,2) = yc_4;          
                                                            fst_mv_oj_coor(4,3) = object_number_count;
                                                            fst_mv_oj_coor(4,4) = jump_number;

                                                            fst_mv_oj_coor(5,1) = xc_5;
                                                            fst_mv_oj_coor(5,2) = yc_5;
                                                            fst_mv_oj_coor(5,3) = object_number_count;
                                                            fst_mv_oj_coor(5,4) = jump_number;

                                                            if( abs(angle_final_th) > tangent_control)
                                                                angle_information(object_number_count,1)=(abs(abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(1,1)))*180/pi)-abs(atan((fst_mv_oj_coor(2,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(2,1)-fst_mv_oj_coor(1,1)))*180/pi))+abs(abs(atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(2,1)))*180/pi)-abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(2,1)))*180/pi))+abs(abs(atan((fst_mv_oj_coor(5,2)-fst_mv_oj_coor(3,2))/(fst_mv_oj_coor(5,1)-fst_mv_oj_coor(3,1)))*180/pi)-abs(atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(3,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(3,1)))*180/pi))) + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) )*priority_degree_rate;
                                                            else
                                                                angle_information(object_number_count,1)=(abs(atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(1,1)))*180/pi-atan((fst_mv_oj_coor(2,2)-fst_mv_oj_coor(1,2))/(fst_mv_oj_coor(2,1)-fst_mv_oj_coor(1,1)))*180/pi)+abs(atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(2,1)))*180/pi-atan((fst_mv_oj_coor(3,2)-fst_mv_oj_coor(2,2))/(fst_mv_oj_coor(3,1)-fst_mv_oj_coor(2,1)))*180/pi)+abs(atan((fst_mv_oj_coor(5,2)-fst_mv_oj_coor(3,2))/(fst_mv_oj_coor(5,1)-fst_mv_oj_coor(3,1)))*180/pi-atan((fst_mv_oj_coor(4,2)-fst_mv_oj_coor(3,2))/(fst_mv_oj_coor(4,1)-fst_mv_oj_coor(3,1)))*180/pi)) + ( abs(dist_thresh_sec - dist_thresh_third) + abs(dist_thresh_third-dist_rest_th) +abs(dist_rest_th-dist_final_th) )*priority_degree_rate;
                                                            end

                                                            %finish fifth
                                                            increase_number=keep_in_number_4;
                                                            ans_xyz5= fifth_point(xc_3,yc_3,xc_4,yc_4,xc_5,yc_5,angle_final_th,dist_final_th,decide_slow,threshold_diff_dislocation,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,increase_distance_threshold,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,priority_zero);
                                                            if(~isnan(ans_xyz5))
                                                                keep_in_number_4=increase_number;
                                                                fst_mv_oj_coor(6,1)=ans_xyz5(1,1);
                                                                fst_mv_oj_coor(6,2)=ans_xyz5(1,2);
                                                                fst_mv_oj_coor(6,3)=object_number_count;
                                                                fst_mv_oj_coor(6,4) = jump_number;
                                                            end
                                                            %finish sixth

                                                            %-----prevent duplication----
                                                            [fast_total_coor,stop_search,object_number_count]=prevent_duple_changed(fst_mv_oj_coor,fast_total_coor,object_number_count,stop_search); 
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