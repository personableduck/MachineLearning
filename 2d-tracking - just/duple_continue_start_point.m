function ans_xyz= duple_continue_start_point(origin_size_fst,xyz,fast_total_coor,fast_obj_start,duple_mat,increase_number,distance_threshold_fast,priority_degree_rate,threshold_diff_dislocation,priority_zero,decide_slow,decide_slow_angle,tangent_control,control_vecter_angle,multilple_short_angle,increase_distance_threshold)

count_passe=0;

distance_fast_mat=0;
estm_ft_size=size(fast_total_coor{1,duple_mat});
if( estm_ft_size(1,1) >= origin_size_fst )
    [row_f1,col_f1,check_f1]=find( xyz{1,fast_obj_start+increase_number}(:,1) < fast_total_coor{1,duple_mat}(origin_size_fst ,1) + distance_threshold_fast & xyz{1,fast_obj_start+increase_number}(:,1) > fast_total_coor{1,duple_mat}(origin_size_fst,1) - distance_threshold_fast);
    xc_1=fast_total_coor{1,duple_mat}(origin_size_fst,1);
    yc_1=fast_total_coor{1,duple_mat}(origin_size_fst,2);
    xc_back_1=fast_total_coor{1,duple_mat}(origin_size_fst-1,1);
    yc_back_1=fast_total_coor{1,duple_mat}(origin_size_fst-1,2);
    xc_back_2=fast_total_coor{1,duple_mat}(origin_size_fst-2,1);
    yc_back_2=fast_total_coor{1,duple_mat}(origin_size_fst-2,2);
end

if(check_f1 == 1)  
    maxsize_fast2_xyz=size(row_f1); %read size
    distance_count1=0;
    for fast_second_mt=1:maxsize_fast2_xyz(1,1) %for fast sperm     
        estm_ft_size=size(fast_total_coor{1,duple_mat});
        if( estm_ft_size(1,1) >= origin_size_fst)
            xx_1 = xyz{1,fast_obj_start+increase_number}(row_f1(fast_second_mt),1) - xc_1; % comparing x coordinates   
            yy_1 = xyz{1,fast_obj_start+increase_number}(row_f1(fast_second_mt),2) - yc_1; % comparing y coordinates
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
            angle_fast_test = atan((xyz{1,fast_obj_start+increase_number}(distance_fast_mat(row_dist(try_loop),2),2)-yc_back_1)/(xyz{1,fast_obj_start+increase_number}(distance_fast_mat(row_dist(try_loop),2),1)-xc_back_1))*180/pi; % unit is degree   
            priority_angle_order(try_loop,1)=abs(abs(angle_fast_test)-abs(before_angle))/priority_degree_rate + distance_fast_mat(row_dist(try_loop));                  
            if( (xyz{1,fast_obj_start+increase_number}(distance_fast_mat(row_dist(try_loop),2),2)-yc_back_1) == 0 & (xyz{1,fast_obj_start+increase_number}(distance_fast_mat(row_dist(try_loop),2),1)-xc_back_1) == 0 )
                priority_angle_order(try_loop,1) = priority_zero; 
            end
            priority_angle_order(try_loop,2)=try_loop;%give priority for smaller angle!!
            [temp_pr1,ord_pr1] = sort(priority_angle_order(:,1));
            priority_angle_order = priority_angle_order(ord_pr1,:);
        end
        for try_loop2 = 1: row_dist_size(1,1) 
            pass_numfil=0;
            xc_2=xyz{1,fast_obj_start+increase_number}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1);
            yc_2=xyz{1,fast_obj_start+increase_number}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2);
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
%             if( posit_nega == 0 )                     
%                 posit_nega = posit_nega2;
%                 xc_back_1=xc_back_2;
%                 yc_back_1=yc_back_2;
%             end
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
%             if(posit_nega2 == 0)
%                 pass_numfil=1;
%                 increase_check=0;
%             end
            if(increase_check == 0)       
                if(pass_numfil == 1) 
                    angle_fast = atan((yc_2-yc_1)/(xc_2-xc_1))*180/pi;
                    if( (yc_2-yc_1) == 0 & (xc_2-xc_1) == 0)
                        angle_fast = 0;
                    end
                    dist_thresh_sec = distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)));
                    count_passe=count_passe+1;
                    ans_xyz(count_passe,1)=xc_back_1;
                    ans_xyz(count_passe,2)=yc_back_1;
                    ans_xyz(count_passe,3)=xc_1;
                    ans_xyz(count_passe,4)=yc_1;
                    ans_xyz(count_passe,5)=xc_2;
                    ans_xyz(count_passe,6)=yc_2;   
                    ans_xyz(count_passe,7)=angle_fast;  
                    ans_xyz(count_passe,8)=dist_thresh_sec;  
                end
            end    
        end
    end
    if(count_passe==0)
        ans_xyz=nan;
    end
else
    ans_xyz=nan;
end


end