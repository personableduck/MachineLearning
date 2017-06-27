function ans_xyz= new_search1_slow_start_point(xyz_11,fast_obj_start,new_search1,increase_number,distance_threshold_fast,priority_degree_rate,priority_zero)

count_passe=0;
distance_fast_mat=0;    
[row_f1,col_f1,check_f1]=find( xyz_11{1,fast_obj_start+increase_number}(:,1) < xyz_11{1,fast_obj_start}(new_search1,1) + distance_threshold_fast & xyz_11{1,fast_obj_start+increase_number}(:,1) > xyz_11{1,fast_obj_start}(new_search1,1) - distance_threshold_fast);   
xc_1=xyz_11{1,fast_obj_start}(new_search1,1);
yc_1=xyz_11{1,fast_obj_start}(new_search1,2);

if(check_f1 == 1)  
    maxsize_fast2_xyz=size(row_f1); %read size
    distance_count1=0;
    for fast_second_mt=1:maxsize_fast2_xyz(1,1) %for fast sperm
        xx_1 = xyz_11{1,fast_obj_start+increase_number}(row_f1(fast_second_mt),1) - xc_1; % comparing x coordinates   
        yy_1 = xyz_11{1,fast_obj_start+increase_number}(row_f1(fast_second_mt),2) - yc_1; % comparing y coordinates
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
        vector_check=nan;
        for try_loop = 1: row_dist_size(1,1)
            priority_angle_order(try_loop,1)=distance_fast_mat(row_dist(try_loop))*priority_degree_rate;
            if(priority_angle_order(try_loop,1) == 0)
                priority_angle_order(try_loop,1) = priority_zero; 
            end
            priority_angle_order(try_loop,2)=try_loop;%give priority for smaller distance!!
            [temp_pr1,ord_pr1] = sort(priority_angle_order(:,1));
            priority_angle_order = priority_angle_order(ord_pr1,:);
        end
        for try_loop2 = 1: row_dist_size(1,1)
            pass_numfil=0;
            xc_2=xyz_11{1,fast_obj_start+increase_number}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),1);
            yc_2=xyz_11{1,fast_obj_start+increase_number}(distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)),2),2);
            %decide vector andgle
            %each magnitude
            angle_fast = atan((yc_2-yc_1)/(xc_2-xc_1))*180/pi; % unit is degree
            if( (yc_2-yc_1) == 0 & (xc_2-xc_1) == 0 )
                angle_fast = 0;
            end
            dist_thresh_sec = distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)));
            if(angle_fast <= 45 & angle_fast >= -45) %vetorc
                vector_check=1; % x vetorc
                diff_vector = xc_2 - xc_1;
                posit_nega=sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
            else
                vector_check=0; % y vetorc
                diff_vector = yc_2 - yc_1; 
                posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
            end
            posit_nega2 = posit_nega;
            pass_numfil=1;                     
            increase_check=0; 
            if(increase_check == 0)       
                if(pass_numfil == 1) 
                end
            end
            count_passe=count_passe+1;
            ans_xyz(count_passe,1)=xc_1;
            ans_xyz(count_passe,2)=yc_1;
            ans_xyz(count_passe,3)=xc_2;
            ans_xyz(count_passe,4)=yc_2;
            ans_xyz(count_passe,5)=posit_nega;
            ans_xyz(count_passe,6)=vector_check;
            ans_xyz(count_passe,7)=angle_fast;
            ans_xyz(count_passe,8)=dist_thresh_sec;
        end
    end
    if(count_passe==0)
        ans_xyz=nan;
    end
else
    ans_xyz=nan;
end

end