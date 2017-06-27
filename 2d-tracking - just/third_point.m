function ans_xyz3= third_point(xc_1,yc_1,xc_2,yc_2,xc_3,yc_3,angle_secd_fast,dist_thresh_third,decide_slow,threshold_diff_dislocation,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,increase_distance_threshold,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,priority_zero)
count_passe=0;
distance_thrd_mat=0;
[row_f3,col_f3,check_f3]=find( xyz{1,fast_obj_start+increase_number}(:,1) < xc_3 + distance_threshold_fast & xyz{1,fast_obj_start+increase_number}(:,1) > xc_3 - distance_threshold_fast);

if(check_f3 == 1)
    maxsize_fast_xyz3=size(row_f3);
    distance_count3=0;
    for dist_thrd_loop = 1: maxsize_fast_xyz3(1,1)                                   
        xx_3 = xyz{1,fast_obj_start+increase_number}(row_f3(dist_thrd_loop),1) - xc_3;
        yy_3 = xyz{1,fast_obj_start+increase_number}(row_f3(dist_thrd_loop),2) - yc_3;
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
            angle_rest_th_test=atan((xyz{1,fast_obj_start+increase_number}(distance_thrd_mat(row_third(find_rest_loop),2),2)-yc_2)/(xyz{1,fast_obj_start+increase_number}(distance_thrd_mat(row_third(find_rest_loop),2),1)-xc_2))*180/pi;                                                           
            priority_angle_order3(find_rest_loop,1)=abs(abs(angle_rest_th_test)-abs(angle_secd_fast)) + distance_thrd_mat(row_third(find_rest_loop))*priority_degree_rate;
            if( (xyz{1,fast_obj_start+increase_number}(distance_thrd_mat(row_third(find_rest_loop),2),2)-yc_2) == 0 & (xyz{1,fast_obj_start+increase_number}(distance_thrd_mat(row_third(find_rest_loop),2),1)-xc_2) ==0 )
                priority_angle_order3(find_rest_loop,1) = priority_zero;
            end
            priority_angle_order3(find_rest_loop,2)=find_rest_loop;%give priority for smaller angle!!
            [temp_pr1,ord_pr1] = sort(priority_angle_order3(:,1));
            priority_angle_order3 = priority_angle_order3(ord_pr1,:);
        end
        for find_rest_loop2 = 1: row_third_size(1,1) 
            pass_numfil=0; 
            diff_vector2=0;
            posit_nega2=nan;
            xc_4=xyz{1,fast_obj_start+increase_number}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),1);
            yc_4=xyz{1,fast_obj_start+increase_number}(distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)),2),2);
            if(angle_secd_fast <= 45 & angle_secd_fast >= -45) %vetorc
                vector_check=1; % x vetorc
                diff_vector = xc_3-xc_2;
                posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
            else 
                vector_check=0; % y vetorc
                diff_vector = yc_3-yc_2;
                posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
            end
            if(vector_check == 1)       
                diff_vector2= xc_4 - xc_3;
                posit_nega2=sign(diff_vector2);                                           
            elseif(vector_check == 0)                                            
                diff_vector2= yc_4 - yc_3;                                            
                posit_nega2=sign(diff_vector2);
            end
%             if( posit_nega == 0 )
%                 xc_2=xc_1;
%                 yc_2=yc_1;
%                 posit_nega = posit_nega2;        
%             end
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
%                 if(posit_nega2 == 0)
%                     pass_numfil=1;
%                     increase_check=0;
%                 end
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
                        % finish fourth
                        count_passe=count_passe+1;
                        ans_xyz3(count_passe,1)=xc_4;
                        ans_xyz3(count_passe,2)=yc_4; 
                        ans_xyz3(count_passe,3)=angle_rest_th;
                        ans_xyz3(count_passe,4)=dist_rest_th;
                    end
                end
            end
        end
    end
    if(count_passe==0)
        ans_xyz3=nan;
    end
else
    ans_xyz3=nan;
end

end