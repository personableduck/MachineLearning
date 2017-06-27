function ans_xyz5= fifth_point(xc_3,yc_3,xc_4,yc_4,xc_5,yc_5,angle_final_th,dist_final_th,decide_slow,threshold_diff_dislocation,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,increase_distance_threshold,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,priority_zero)
count_passe=0;
distance_six_mat=0;
[row_f5,col_f5,check_f5]=find( xyz{1,fast_obj_start+increase_number}(:,1) < xc_5 + distance_threshold_fast & xyz{1,fast_obj_start+increase_number}(:,1) > xc_5 - distance_threshold_fast);       

if(check_f5 == 1)
    stop_search2=0;
    maxsize_fast_xyz6=size(row_f5);
    distance_count6=0;
    for dist_six_loop = 1: maxsize_fast_xyz6(1,1)                                   
        xx_6 = xyz{1,fast_obj_start+increase_number}(row_f5(dist_six_loop),1) - xc_5;
        yy_6 = xyz{1,fast_obj_start+increase_number}(row_f5(dist_six_loop),2) - yc_5;                                                                                                                                       
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
            angle_six_th_test=atan((xyz{1,fast_obj_start+increase_number}(distance_six_mat(row_six(find_six_loop),2),2)-yc_4)/(xyz{1,fast_obj_start+increase_number}(distance_six_mat(row_six(find_six_loop),2),1)-xc_4))*180/pi;                                                           
            priority_angle_order6(find_six_loop,1)=abs(abs(angle_final_th)-abs(angle_six_th_test)) + distance_six_mat(row_six(find_six_loop))*priority_degree_rate;
            if( (xyz{1,fast_obj_start+increase_number}(distance_six_mat(row_six(find_six_loop),2),2)-yc_4) == 0 & (xyz{1,fast_obj_start+increase_number}(distance_six_mat(row_six(find_six_loop),2),1)- xc_4) )
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
                xc_6=xyz{1,fast_obj_start+increase_number}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),1);
                yc_6=xyz{1,fast_obj_start+increase_number}(distance_six_mat(row_six(priority_angle_order6(find_six_loop2,2)),2),2);
                if(angle_final_th <= 45 & angle_final_th >= -45) %vetorc
                    vector_check=1; % x vetorc
                    diff_vector = xc_5-xc_4;
                    posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                else 
                    vector_check=0; % y vetorc
                    diff_vector = yc_5-yc_4;
                    posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                end
                if(vector_check == 1)       
                    diff_vector2= xc_6 - xc_5;
                    posit_nega2=sign(diff_vector2);                                           
                elseif(vector_check == 0)                                            
                    diff_vector2= yc_6 - yc_5;
                    posit_nega2=sign(diff_vector2);
                end
%                 if( posit_nega == 0 )
%                     xc_4=xc_3;
%                     yc_4=yc_3;       
%                     posit_nega = posit_nega2;        
%                 end
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
%                     if(posit_nega2 == 0)
%                         pass_numfil=1;
%                         increase_check=0;
%                     end
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
                            stop_search2=1;
                            count_passe=count_passe+1;
                            ans_xyz5(count_passe,1)=xc_6;
                            ans_xyz5(count_passe,2)=yc_6; 
                        end
                    end
                end
            end
        end
    end 
    if(count_passe==0)
        ans_xyz5=nan;
    end
else
    ans_xyz5=nan;
end


end
