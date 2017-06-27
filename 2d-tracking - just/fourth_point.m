function ans_xyz4= fourth_point(xc_2,yc_2,xc_3,yc_3,xc_4,yc_4,angle_rest_th,dist_rest_th,decide_slow,threshold_diff_dislocation,xyz,fast_obj_start,inside_loop,increase_number,distance_threshold_fast,priority_degree_rate,increase_distance_threshold,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,priority_zero)
count_passe=0;
distance_four_mat=0;
[row_f4,col_f4,check_f4]=find( xyz{1,fast_obj_start+increase_number}(:,1) < xc_4 + distance_threshold_fast & xyz{1,fast_obj_start+increase_number}(:,1) > xc_4 - distance_threshold_fast);       

if(check_f4 == 1)
    maxsize_fast_xyz4=size(row_f4);
    distance_count4=0;
    for dist_four_loop = 1: maxsize_fast_xyz4(1,1)
        xx_4 = xyz{1,fast_obj_start+increase_number}(row_f4(dist_four_loop),1) - xc_4;    
        yy_4 = xyz{1,fast_obj_start+increase_number}(row_f4(dist_four_loop),2) - yc_4;                                
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
            angle_final_th_test=atan((xyz{1,fast_obj_start+increase_number}(distance_four_mat(row_four(find_four_loop),2),2)-yc_3)/(xyz{1,fast_obj_start+increase_number}(distance_four_mat(row_four(find_four_loop),2),1)-xc_3))*180/pi;
            priority_angle_order4(find_four_loop,1)=abs(abs(angle_rest_th)-abs(angle_final_th_test)) + distance_four_mat(row_four(find_four_loop))*priority_degree_rate;
            if( (xyz{1,fast_obj_start+increase_number}(distance_four_mat(row_four(find_four_loop),2),2)-yc_3) == 0 & (xyz{1,fast_obj_start+increase_number}(distance_four_mat(row_four(find_four_loop),2),1)-xc_3) == 0)
                priority_angle_order4(find_four_loop,1)=priority_zero;
            end
            priority_angle_order4(find_four_loop,2)=find_four_loop;%give priority for smaller angle!!
            [temp_pr1,ord_pr1] = sort(priority_angle_order4(:,1));
            priority_angle_order4 = priority_angle_order4(ord_pr1,:);
        end
        for find_four_loop2 = 1: row_four_size(1,1)                                                          
            xc_5=xyz{1,fast_obj_start+increase_number}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),1);
            yc_5=xyz{1,fast_obj_start+increase_number}(distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)),2),2);
            pass_numfil=0;
            diff_vector2=0;
            posit_nega2=nan;
            if(angle_rest_th <= 45 & angle_rest_th >= -45) %vetorc
                vector_check=1; % x vetorc
                diff_vector = xc_4-xc_3;
                posit_nega= sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
            else 
                vector_check=0; % y vetorc
                diff_vector = yc_4-yc_3;
                posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
            end
            if(vector_check == 1)       
                diff_vector2= xc_5 - xc_4;      
                posit_nega2=sign(diff_vector2);        
            elseif(vector_check == 0)                                            
                diff_vector2= yc_5 - yc_4;                                            
                posit_nega2=sign(diff_vector2);
            end
%             if( posit_nega == 0 )
%                 xc_3=xc_2;
%                 yc_3=yc_2;
%                 posit_nega = posit_nega2;        
%             end
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
%                 if(posit_nega2 == 0)
%                     pass_numfil=1;
%                     increase_check=0;
%                 end
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

                        % finish fifth
                        count_passe=count_passe+1;
                        ans_xyz4(count_passe,1)=xc_5;
                        ans_xyz4(count_passe,2)=yc_5; 
                        ans_xyz4(count_passe,3)=angle_final_th;
                        ans_xyz4(count_passe,4)=dist_final_th;
                    end
                end
            end
        end
    end
    if(count_passe==0)
        ans_xyz4=nan;
    end
else
    ans_xyz4=nan;
end

                                                                                                  
end