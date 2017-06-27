function ans_xyz2= second_point(xc_back_1,yc_back_1,xc_1,yc_1,xc_2,yc_2,posit_nega,angle_fast,vector_check,dist_thresh_sec,decide_slow,threshold_diff_dislocation,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,increase_distance_threshold,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,priority_zero)
count_passe=0;
distance_repeat_rest=0;
[row_f2,col_f2,check_f2]=find( xyz{1,fast_obj_start+increase_number}(:,1) < xc_2 + distance_threshold_fast & xyz{1,fast_obj_start+increase_number}(:,1) > xc_2 - distance_threshold_fast);

if(check_f2 == 1)
    maxsize_fast_xyz2=size(row_f2); %read size
    distance_count2=0;
    for dist_cmp_loop = 1: maxsize_fast_xyz2(1,1)
        xx_2 = xyz{1,fast_obj_start+increase_number}(row_f2(dist_cmp_loop),1) - xc_2;                              
        yy_2 = xyz{1,fast_obj_start+increase_number}(row_f2(dist_cmp_loop),2) - yc_2;                 
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
            angle_secd_fast_test = atan((xyz{1,fast_obj_start+increase_number}(distance_repeat_rest(row_secd(find_secd_loop),2),2)-yc_1)/(xyz{1,fast_obj_start+increase_number}(distance_repeat_rest(row_secd(find_secd_loop),2),1)-xc_1))*180/pi;
            priority_angle_order2(find_secd_loop,1)=abs(abs(angle_secd_fast_test)-abs(angle_fast)) + distance_repeat_rest(row_secd(find_secd_loop))*priority_degree_rate;
            if( (xyz{1,fast_obj_start+increase_number}(distance_repeat_rest(row_secd(find_secd_loop),2),2)-yc_1) == 0 & (xyz{1,fast_obj_start+increase_number}(distance_repeat_rest(row_secd(find_secd_loop),2),1)-xc_1) == 0 )
                priority_angle_order2(find_secd_loop,1) = priority_zero; 
            end
            priority_angle_order2(find_secd_loop,2)=find_secd_loop;%give priority for smaller angle!!
            [temp_pr1,ord_pr1] = sort(priority_angle_order2(:,1));
            priority_angle_order2 = priority_angle_order2(ord_pr1,:);
        end
        for find_secd_loop2 = 1: row_secd_size(1,1)
            xc_3=xyz{1,fast_obj_start+increase_number}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),1);
            yc_3=xyz{1,fast_obj_start+increase_number}(distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)),2),2);
            pass_numfil=0;
            diff_vector2=0;
            posit_nega2=nan;
            if( fast_obj_start > 1)
                 if(angle_fast <= 45 & angle_fast >= -45) %vetorc
                    vector_check=1; % x vetorc
                    diff_vector = xc_2-xc_1;
                    posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                else 
                    vector_check=0; % y vetorc
                    diff_vector = yc_2-yc_1;
                    posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                 end
            end
            if(vector_check == 1)
                diff_vector2= xc_3-xc_2;
                posit_nega2=sign(diff_vector2);
            elseif (vector_check == 0)
                diff_vector2= yc_3-yc_2;                            
                posit_nega2=sign(diff_vector2);
            end
%             if( posit_nega == 0 )
%                 if( fast_obj_start > 1)
%                     xc_1=xc_back_1;
%                     yc_1=yc_back_1;
%                 end
%                 posit_nega = posit_nega2;        
%             end
            if(posit_nega == posit_nega2)
                angle_fast=atan((yc_2-yc_1)/(xc_2-xc_1))*180/pi;
                if( (yc_2-yc_1) == 0 & (xc_2-xc_1) == 0 )
                    angle_fast=0;
                end
                dist_thresh_sec=sqrt((xc_2-xc_1)^2+(yc_2-yc_1)^2);
                increase_check=1;
                if(fast_obj_start == 1)
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
%                 if(posit_nega2 == 0)
%                     pass_numfil=1;
%                     increase_check=0;
%                 end
                if(increase_check == 0)
                    angle_secd_fast_cmp=atan((yc_3-yc_2)/(xc_3-xc_2))*180/pi; 
                    if( (yc_3-yc_2) == 0 & (xc_3-xc_2) == 0 )
                        angle_secd_fast_cmp=0;
                    end
                    if(fast_obj_start==1)
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
                        count_passe=count_passe+1;
                        ans_xyz2(count_passe,1)=xc_3;
                        ans_xyz2(count_passe,2)=yc_3; 
                        ans_xyz2(count_passe,3)=angle_secd_fast;
                        ans_xyz2(count_passe,4)=dist_thresh_third;
                    end
                end
            end      
        end
    end
    if(count_passe==0)
        ans_xyz2=nan;
    end
else
    ans_xyz2=nan;
end
  
                     
end