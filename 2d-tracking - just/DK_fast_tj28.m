%New trajectories5 by DUCK-HA HWANG

%------solve diff holo dislocation, change xyz coordinate information.

constant_distance_rate=0.4; % endure 40% difference
threshold_diff_dislocation= 1+constant_distance_rate; %*****check it!! important value (1+1*40%) for combining close coordinates
threshold_diff_dislocation2= 1-constant_distance_rate;
threshold_diff_dislocation_shrt=2; % 200% for slow moving
threshold_diff_dislocation_shrt2=1/(threshold_diff_dislocation_shrt); % 200%
control_vecter_angle = 23; %angle control
multilple_short_angle = 2; %find wiggling moving for slow one
increase_distance_threshold=0.8; %make sure about linear moving for fast
combine_limitation=1; % combine xyz noise distance
tangent_control=60;%tangent absolute value 60 degree
increase_max_limit=2; % limitation for rapid incerease movement

fps=3; %fps means frame per second
distance_threshold_fast = (120/fps)/pixelsize; %******threshold object tracking. need to control. maximum ditance to distinguish as a same object
%unit is um. // I limited that the fasted sperm's speed is 110 um/sec.
%fps means frame per second

decide_slow=(80/fps)/pixelsize; %criteria for slow moving
%decide_slow=11; %distancce that ovelap happen on the differ_holo_stack

%-----rest
max_frame_size=fix((iter-4)/2);
frame_number=2*max_frame_size+1
%-----skip
%distance_threshold_skip=6;
distance_threshold_skip=(40/fps)/pixelsize;

%-----------combine xyz

for combine_compare_loop=1:iter-2
    
    behind_compare=combine_compare_loop+1;
    
    maxsize_cmb_cp=size(xyz{1,combine_compare_loop}); %read size
    
    maxsize_cmb_cp2=size(xyz{1,behind_compare}); %read size

    
    for inft_loop=1:maxsize_cmb_cp(1,1)
        for bhid_loop=1:maxsize_cmb_cp2(1,1)

            xx_shrt = xyz{1,combine_compare_loop}(inft_loop,1) - xyz{1,behind_compare}(bhid_loop,1); % comparing x coordinates
            yy_shrt = xyz{1,combine_compare_loop}(inft_loop,2) - xyz{1,behind_compare}(bhid_loop,2); % comparing y coordinates
        
            distance_shrt_limit= sqrt((xx_shrt^2)+(yy_shrt^2));

            if (distance_shrt_limit <= combine_limitation)         
                xyz{1,behind_compare}(bhid_loop,1) = xyz{1,combine_compare_loop}(inft_loop,1);
                xyz{1,behind_compare}(bhid_loop,2) = xyz{1,combine_compare_loop}(inft_loop,2);
                
            end
        end
    end
    
end

%------change xyz


fast_obj_start=1; %make a loop

%make a combined matrix

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

for change_xyz_loop=1:2

comp_xyz_new=0;
comp_xyz_new2=0;
rsh_number=0;
researching_xyz1=0;
    
comp_xyz_new=xyz{1,change_xyz_loop};
comp_xyz_new2=xyz{1,change_xyz_loop+1};

size_cmp_n=size(comp_xyz_new2);
    
for loop_cmp_del=1:size_cmp_n(1,1)
    
    [row_del_x]=find(comp_xyz_new(:,1) == comp_xyz_new2(loop_cmp_del,1));
    [row_del_y]=find(comp_xyz_new(:,2) == comp_xyz_new2(loop_cmp_del,2));
    
    rr_ss=size(row_del_x);
    for check_loop_dlt=1:rr_ss
        [row_del_f,col_del_f,check_del_f]=find(row_del_y == row_del_x(check_loop_dlt));
            
        if(check_del_f ~= 0)
            rsh_number=rsh_number+1;
            researching_xyz1(rsh_number,1)= row_del_x(check_loop_dlt);
        end
    end
end

maxsize_rsch_xyz1=size(researching_xyz1); %read size
maxsize_behinding=size(xyz{1,change_xyz_loop});

count_change_num=0;
change_nan=0;
distance_strt=0;
skip_strt_dlp=0;

for behinding_loop= 1: maxsize_behinding(1,1)
    
    distance_strt=0;
      
    for find_start_loop=1:maxsize_rsch_xyz1(1,1)
        
        srt_xx=xyz{1,change_xyz_loop}(behinding_loop,1)-xyz{1,change_xyz_loop}(researching_xyz1(find_start_loop,1),1);
        srt_yy=xyz{1,change_xyz_loop}(behinding_loop,2)-xyz{1,change_xyz_loop}(researching_xyz1(find_start_loop,1),2);
        
        [row_skip_dpl,col_skip_dpl,check_skip_dpl]=find(skip_strt_dlp == behinding_loop);
        
        if(check_skip_dpl == 1)
            srt_xx=nan;
            srt_yy=nan;
        end
        
        distance_strt(find_start_loop,1)=sqrt(srt_xx^2+srt_yy^2);
    end
    [row_strt,col_strt,check_strt]=find(distance_strt <= distance_threshold_fast & distance_strt ~= 0);
    if(check_strt == 1)
        row_strt_size=size(row_strt); %read size
        count_skip_strt=0;
        for row_srt_dpl_loop=1:row_strt_size(1,1)
            [row_srt_dpl,col_srt_dpl,check_dpl]=find(researching_xyz1(:,1) == researching_xyz1(row_strt(row_srt_dpl_loop)));
            if(check_dpl == 1)
                row_strt_dpl_size=size(row_srt_dpl); %read size
                for row_strt_dpl_loop=1:row_strt_dpl_size(1,1);
                    xyz_1{1,change_xyz_loop}(researching_xyz1(row_srt_dpl(row_strt_dpl_loop)),:)=nan;
                    count_skip_strt=count_skip_strt+1;
                    skip_strt_dlp(count_skip_strt,1)=researching_xyz1(row_srt_dpl(row_strt_dpl_loop));
                end
            end
        end
        
        if(count_skip_strt==0) 
            row_strt_size=size(row_strt); %read size 
            distance_strt_matrix=0;
            for row_strt_loop=1:row_strt_size(1,1)
                distance_strt_matrix(row_strt_loop,1) = distance_strt(row_strt(row_strt_loop));
            end
            distance_strt_min=min(distance_strt_matrix);
            order_dist_strt=find(distance_strt_matrix==distance_strt_min);       
            xyz_1{1,change_xyz_loop}(researching_xyz1(row_strt(order_dist_strt)),:)=nan;
        end      
    end
end
end

%****************************************rest***************************
object_number_count=0;
cont_rest_th=0;
cont_angle_th=0;
fast_total_coor={0};
angle_information=0;

for rest_start_searching= 1:max_frame_size

    fast_obj_start=2*rest_start_searching-1;
   
    if(fast_obj_start == 1)
        maxsize_fast_xyz=size(xyz_1{1,fast_obj_start}); %read size
        change_inside_loop = maxsize_fast_xyz(1,1);
    else
        change_inside_loop = object_number_count;
    end
    
    for inside_loop=1:change_inside_loop
        stop_search=0;
        distance_fast_mat=0; 
        distance_max_hold=0;
        maxsize_fast2_xyz=size(xyz{1,fast_obj_start+1}); %read size    
        
        for fast_second_mt=1:maxsize_fast2_xyz(1,1) %for fast sperm
        
            if(fast_obj_start == 1)
                
                xx_1 = xyz_1{1,fast_obj_start+1}(fast_second_mt,1) - xyz_1{1,fast_obj_start}(inside_loop,1); % comparing x coordinates   
                yy_1 = xyz_1{1,fast_obj_start+1}(fast_second_mt,2) - xyz_1{1,fast_obj_start}(inside_loop,2); % comparing y coordinates
                
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
                [row_dist,col_dist,check_dist]=find(distance_fast_mat < cont_rest_th(inside_loop,1)*threshold_diff_dislocation_shrt & distance_fast_mat(fast_second_mt,1) > cont_rest_th(inside_loop,1)*threshold_diff_dislocation_shrt2);            
            else
                [row_dist,col_dist,check_dist]=find(distance_fast_mat < cont_rest_th(inside_loop,1)*threshold_diff_dislocation & distance_fast_mat(fast_second_mt,1) > cont_rest_th(inside_loop,1)*threshold_diff_dislocation2);
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
                    priority_angle_order(try_loop,1)=try_loop;
                    priority_angle_order(try_loop,2)=try_loop;%give priority for smaller angle!!
                end
            else 
                for try_loop = 1: row_dist_size(1,1)
                    angle_fast_test = atan((xyz{1,fast_obj_start+1}(row_dist(try_loop),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))/(xyz{1,fast_obj_start+1}(row_dist(try_loop),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi; % unit is degree   
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
                    angle_fast = atan((xyz_1{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2)-xyz_1{1,fast_obj_start}(inside_loop,2))/(xyz_1{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)-xyz_1{1,fast_obj_start}(inside_loop,1)))*180/pi; % unit is degree
                    dist_thresh_sec = distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)));
                
                    if(angle_fast <= 45 & angle_fast >= -45) %vetorc
                        vector=1; % x vetorc
                        diff_vector = xyz_1{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1) - xyz_1{1,fast_obj_start}(inside_loop,1);
                        posit_nega=sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                    else
                        vector=0; % y vetorc
                        diff_vector = xyz_1{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2) - xyz_1{1,fast_obj_start}(inside_loop,2); 
                        posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                    end
                    %first xyz{1,fast_obj_start}(inside_loop,1), xyz{1,fast_obj_start}(inside_loop,2)      
                    %second xyz(row_dist(try_loop),1),xyz(row_dist(try_loop),2)
                    
                    posit_nega2 = posit_nega;
                    pass_numfil=1;
                
                else

                    if(cont_angle_th(inside_loop,1) <= 45 & cont_angle_th(inside_loop,1) >= -45) %vetorc
                        vector=1; % x vetorc
                        diff_vector = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1) - fast_total_coor{1,inside_loop}(fast_obj_start,1);
                        posit_nega=sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                    else 
                        vector=0; % y vetorc
                        diff_vector = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2) - fast_total_coor{1,inside_loop}(fast_obj_start,2); 
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
                    
                    angle_fast_cmp=atan((xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))/(xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi;
                    
                    if( cont_rest_th(inside_loop,1) > decide_slow)           
                        if(abs(cont_angle_th(inside_loop,1)) > tangent_control)              
                            if(abs(cont_angle_th(inside_loop,1))+control_vecter_angle > abs(angle_fast_cmp) & abs(cont_angle_th(inside_loop,1))-control_vecter_angle < abs(angle_fast_cmp)) % unit is degree
                                pass_numfil=1;
                            end
                            
                        else
                            
                            if(cont_angle_th(inside_loop,1)+control_vecter_angle > angle_fast_cmp & cont_angle_th(inside_loop,1)-control_vecter_angle < angle_fast_cmp) % unit is degree
                            pass_numfil=1;
                            end
                        end
                        
                    else
                        if(abs(cont_angle_th(inside_loop,1))+control_vecter_angle*multilple_short_angle > abs(angle_fast_cmp) & abs(cont_angle_th(inside_loop,1))-control_vecter_angle*multilple_short_angle < abs(angle_fast_cmp)) % unit is degree
                            pass_numfil=1;
                        end
                    end
                end
                    if(posit_nega == posit_nega2)
                        pass_value=1; 
                        if(fast_obj_start == 1)                         
                            pass_value=0; 
                        else
                            
                            if((distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)))+cont_rest_th(inside_loop,1))*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start-1,1))^2+(xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start-1,2))^2 ))          
                                pass_value=0; 
                            else   
                                pass_value=1;
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
                                    [row_secd,col_secd,check_secd]=find(distance_repeat_rest < dist_thresh_sec*threshold_diff_dislocation_shrt & distance_repeat_rest > dist_thresh_sec*threshold_diff_dislocation_shrt2);
                                else
                                    [row_secd,col_secd,check_secd]=find(distance_repeat_rest < dist_thresh_sec*threshold_diff_dislocation & distance_repeat_rest > dist_thresh_sec*threshold_diff_dislocation2);
                                end
                                
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
                                            
                                            if(fast_obj_start == 1)
                                                if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz_1{1,fast_obj_start}(inside_loop,1))^2+(xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz_1{1,fast_obj_start}(inside_loop,2))^2 ))
                                                    pass_value=0;
                                                else
                                                    pass_value=1;
                                                end
                                                
                                            else
                                                if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1))^2+(xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))^2 ))
                                                    pass_value=0;
                                                else
                                                    pass_value=1;
                                                end
                                            end

                                            if(pass_value == 0)                                       
                                                angle_secd_fast_cmp=atan((xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2))/(xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1)))*180/pi; 

                                                if(dist_thresh_sec > decide_slow)
                                                    if(abs(angle_fast)>tangent_control)                                                     
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
                                                    if( (dist_thresh_third/dist_thresh_sec) > increase_max_limit )
                                                        distance_max_hold=dist_thresh_third;
                                                    elseif( (dist_thresh_third/dist_thresh_sec) < (1/increase_max_limit) )
                                                        distance_max_hold=dist_thresh_sec;
                                                    end
                                                    % finish third
                                                    distance_thrd_mat=0;
                                                    maxsize_fast_xyz3=size(xyz{1,fast_obj_start+3});

                                                    for dist_thrd_loop = 1: maxsize_fast_xyz3(1,1)                                   
                                                        xx_3 = xyz{1,fast_obj_start+3}(dist_thrd_loop,1) - xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1);
                                                        yy_3 = xyz{1,fast_obj_start+3}(dist_thrd_loop,2) - xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2);
                                                        distance_thrd_mat(dist_thrd_loop,1) = sqrt((xx_3^2)+(yy_3^2));                           
                                                    end
  
                                                    if(dist_thresh_third <= decide_slow)
                                                        if( distance_max_hold ~= 0)
                                                            [row_third,col_third,check_third]=find(distance_thrd_mat < distance_max_hold*threshold_diff_dislocation);
                                                        else
                                                        [row_third,col_third,check_third]=find(distance_thrd_mat < dist_thresh_third*threshold_diff_dislocation_shrt & distance_thrd_mat > dist_thresh_third*threshold_diff_dislocation_shrt2);
                                                        end
                                                    else 
                                                        [row_third,col_third,check_third]=find(distance_thrd_mat < dist_thresh_third*threshold_diff_dislocation & distance_thrd_mat > dist_thresh_third*threshold_diff_dislocation2);
                                                    end
     
                                                    if (check_third == 1)

                                                        row_third_size=size(row_third);
                                                        angle_rest_th=0;
                                                        dist_rest_th=0;                         
                                                        priority_angle_order3=0;
                                
                                                        for find_rest_loop = 1: row_third_size(1,1)
                                                            angle_rest_th_test=atan((xyz{1,fast_obj_start+3}(row_third(find_rest_loop),2)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2))/(xyz{1,fast_obj_start+3}(row_third(find_rest_loop),1)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)))*180/pi;                                                           
                                                            priority_angle_order3(find_rest_loop,1)=abs(angle_secd_fast-angle_rest_th_test);
                                                            priority_angle_order3(find_rest_loop,2)=find_rest_loop;%give priority for smaller angle!!
                                                            
                                                            if(abs(angle_secd_fast)>=75)
                                                                priority_angle_order3(find_rest_loop,1)=abs(abs(angle_secd_fast)-abs(angle_rest_th_test));
                                                            end
                                                       
                                                            [temp_pr1,ord_pr1] = sort(priority_angle_order3(:,1));
                                                            priority_angle_order3 = priority_angle_order3(ord_pr1,:);
                                                        end
                                                        
                                                        for find_rest_loop2 = 1: row_third_size(1,1)
                                                            pass_numfil=0;  
                                                            if(stop_search == 0)
                                                                diff_vector2=0;
                                                                posit_nega2=nan;
                                                                if(vector == 1)       
                                                                    diff_vector2= xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1) - xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1);
                                                                    posit_nega2=sign(diff_vector2);                                           
                                                                elseif(vector == 0)                                            
                                                                    diff_vector2= xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2) - xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2);                                            
                                                                    posit_nega2=sign(diff_vector2);
                                                                end

                                                                if(posit_nega == posit_nega2)
                                                                    pass_value=1;
                                                                    if(fast_obj_start == 1)    
                                                                        if((distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)))+dist_thresh_third)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1))^2+(xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2))^2 ))
                                                                            pass_value=0;                                                     
                                                                        else            
                                                                            pass_value=1;
                                                                        end         
                                                                    else
                                                                        
                                                                        if((distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)))+dist_thresh_third)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1))^2+(xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2)-xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2))^2 ))                                                                                 
                                                                            pass_value=0;                                                     
                                                                        else            
                                                                            pass_value=1;
                                                                        end  
                                                                    end

                                                                    if(pass_value == 0)                                                      
                                                                        angle_rest_th_cmp=atan((xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2))/(xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)))*180/pi;
                                                                        
                                                                        if( dist_thresh_third > decide_slow)
                                                                            if(abs(angle_secd_fast)>tangent_control)    
                                                                                if(abs(angle_secd_fast)+control_vecter_angle > abs(angle_rest_th_cmp) & abs(angle_secd_fast)-control_vecter_angle < abs(angle_rest_th_cmp))  
                                                                                    pass_numfil=1;
                                                                                end
                                                                                
                                                                            else
                                                                                
                                                                                if(angle_secd_fast+control_vecter_angle > angle_rest_th_cmp & angle_secd_fast-control_vecter_angle < angle_rest_th_cmp)  
                                                                                    pass_numfil=1;
                                                                                end
                                                                            end
                                                                        else
                                                                            if(abs(angle_secd_fast)+control_vecter_angle*multilple_short_angle > abs(angle_rest_th_cmp) & abs(angle_secd_fast)-control_vecter_angle*multilple_short_angle < abs(angle_rest_th_cmp))  
                                                                                pass_numfil=1;
                                                                            end  
                                                                        end
                                                                        
                                                                        if( pass_numfil == 1)
                                                                            angle_rest_th=atan((xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2))/(xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)))*180/pi;                 
                                                                            dist_rest_th=distance_thrd_mat(row_third(priority_angle_order3(find_rest_loop2,2)));
                                                                            %fourth xyz(row_third(find_rest_loop),1),xyz(row_third(find_rest_loop),2)
                                                                            
                                                                            if( distance_max_hold == 0)
                                                                                if( (dist_rest_th/dist_thresh_third) > increase_max_limit )
                                                                                    distance_max_hold=dist_rest_th;
                                                                                elseif( (dist_rest_th/dist_thresh_third) < (1/increase_max_limit) )
                                                                                    distance_max_hold=dist_thresh_third;
                                                                                end
                                                                            end
                                                                            
                                                                            distance_four_mat=0;
                                                                            maxsize_fast_xyz4=size(xyz{1,fast_obj_start+4});
                                                                            
                                                                            for dist_four_loop = 1: maxsize_fast_xyz4(1,1)                                                       
                                                                                xx_4 = xyz{1,fast_obj_start+4}(dist_four_loop,1) - xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1);    
                                                                                yy_4 = xyz{1,fast_obj_start+4}(dist_four_loop,2) - xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2);     
                                                                                distance_four_mat(dist_four_loop,1) = sqrt((xx_4^2)+(yy_4^2));                           
                                                                            end

                                                                            if(dist_rest_th <= decide_slow)
                                                                                if( distance_max_hold ~= 0)
                                                                                    [row_four,col_four,check_four]=find(distance_four_mat < distance_max_hold*threshold_diff_dislocation);
                                                                                else
                                                                                    [row_four,col_four,check_four]=find(distance_four_mat < dist_rest_th*threshold_diff_dislocation_shrt & distance_four_mat > dist_rest_th*threshold_diff_dislocation_shrt2); 
                                                                                end
                                                                            else
                                                                                [row_four,col_four,check_four]=find(distance_four_mat < dist_rest_th*threshold_diff_dislocation & distance_four_mat > dist_rest_th*threshold_diff_dislocation2);   
                                                                            end

                                                                                if (check_four == 1)
                                                                                    row_four_size=size(row_four);
                                                                                    angle_four_th=0;
                                                                                    dist_four_th=0; 
                                                                                    priority_angle_order4=0;
                                                                                    
                                                                                    for find_four_loop = 1: row_four_size(1,1)
                                                                                        angle_final_th_test=atan((xyz{1,fast_obj_start+4}(row_four(find_four_loop),2)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2))/(xyz{1,fast_obj_start+4}(row_four(find_four_loop),1)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)))*180/pi;
                                                                                        
                                                                                        priority_angle_order4(find_four_loop,1)=abs(angle_rest_th-angle_final_th_test);
                                                                                        priority_angle_order4(find_four_loop,2)=find_four_loop;%give priority for smaller angle!!
                                                                                        if(abs(angle_rest_th)>=75)
                                                                                            priority_angle_order4(find_four_loop,1)=abs(abs(angle_rest_th)-abs(angle_final_th_test));
                                                                                        end

                                                                                        [temp_pr1,ord_pr1] = sort(priority_angle_order4(:,1));
                                                                                        priority_angle_order4 = priority_angle_order4(ord_pr1,:);
                                                                                    end
                                                                                    
                                                                                    for find_four_loop2 = 1: row_four_size(1,1)                                                          
                                                                                      if(stop_search == 0)  
                                                                                        pass_numfil=0;
                                                                                        diff_vector2=0;
                                                                                        posit_nega2=nan;
                                                                                        
                                                                                        if(vector == 1)       
                                                                                            diff_vector2= xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1) - xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1);      
                                                                                            posit_nega2=sign(diff_vector2);        
                                                                                        elseif(vector == 0)                                            
                                                                                            diff_vector2= xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2) - xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2);                                            
                                                                                            posit_nega2=sign(diff_vector2);
                                                                                        end
                                                                                        
                                                                                        if(posit_nega == posit_nega2)
                                                                                            
                                                                                            pass_value=1;
                                                                                            if((distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)))+dist_rest_th)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1))^2+(xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2)-xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2))^2 ))       
                                                                                                pass_value=0;    
                                                                                            else
                                                                                                pass_value=1;      
                                                                                            end
    
                                                                                            if(pass_value == 0)                                                                                             
                                                                                                angle_final_th_cmp=atan((xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2))/(xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)))*180/pi;

                                                                                                if( dist_rest_th > decide_slow)
                                                                                                    if(abs(angle_rest_th)>tangent_control)
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
                                                                                                    angle_four_th=atan((xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),2)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2))/(xyz{1,fast_obj_start+4}(row_four(priority_angle_order4(find_four_loop2,2)),1)-xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1)))*180/pi;                                    
                                                                                                    dist_four_th=distance_four_mat(row_four(priority_angle_order4(find_four_loop2,2)));
                                                                                                    
                                                                                                    if(fast_obj_start == 1)
                                                                                                        object_number_count = object_number_count+1;
                                                                                                        
                                                                                                        fst_mv_oj_coor=0;
                                        
                                                                                                        fst_mv_oj_coor(1,1) = xyz_1{1,fast_obj_start}(inside_loop,1);                              
                                                                                                        fst_mv_oj_coor(1,2) = xyz_1{1,fast_obj_start}(inside_loop,2);          
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
                                                                                                        
                                                                                                        angle_information(object_number_count,1)=abs(abs(angle_fast-angle_secd_fast)+abs(angle_secd_fast-angle_rest_th)+abs(angle_rest_th-angle_four_th));
                                                                                                        
                                                                                                        if( fast_total_coor{1,1} == 0)
                                                                                                            fast_total_coor(1,object_number_count)={fst_mv_oj_coor};
                                                                                                            cont_rest_th(object_number_count,1)=dist_thresh_third;
                                                                                                            cont_angle_th(object_number_count,1)=angle_secd_fast;
                                                                                                            stop_search=1; %prevent duplicate trajectories
                                                                                                        end                 
                                                                                                        %------prevent duplication--------
                                                                                                        maxsize_fstc=size(fast_total_coor);
                                                                                                        if (fast_total_coor{1,1} ~= 0)
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
                                                                                                                    object_number_count = object_number_count-1; %duplication happen
                                                                                                                else
                                                                                                                    fast_total_coor{1,duplication_number} = {0};
                                                                                                                    fast_total_coor{1,duplication_number} = fst_mv_oj_coor; %smaller angle is right trajectory
                                                                                                                    object_number_count = object_number_count-1;
                                                                                                                end
                                                                                                            elseif(duplication_number == 0)
                                                                                                                fast_total_coor(1,object_number_count)={fst_mv_oj_coor};
                                                                                                                cont_rest_th(object_number_count,1)=dist_thresh_third;
                                                                                                                cont_angle_th(object_number_count,1)=angle_secd_fast;
                                                                                                                stop_search=1; %prevent duplicate trajectories
                                                                                                            end
                                                                                                        end
         
                                                                                                    else

                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start+1,1) = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),1);
                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start+1,2) = xyz{1,fast_obj_start+1}(row_dist(priority_angle_order(try_loop2,2)),2);       
                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start+1,3) = inside_loop;
                                                             
                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start+2,1) = xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),1);
                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start+2,2) = xyz{1,fast_obj_start+2}(row_secd(priority_angle_order2(find_secd_loop2,2)),2);          
                                                                                                        fast_total_coor{1,inside_loop}(fast_obj_start+2,3) = inside_loop;

                                                                                                        if( dist_rest_th > decide_slow)
                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+3,1) = xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),1);
                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+3,2) = xyz{1,fast_obj_start+3}(row_third(priority_angle_order3(find_rest_loop2,2)),2);          
                                                                                                            fast_total_coor{1,inside_loop}(fast_obj_start+3,3) = inside_loop;
                                                                                                        end

                                                                                                        cont_rest_th(inside_loop,1)=dist_thresh_third;
                                                                                                        cont_angle_th(inside_loop,1)=angle_secd_fast;
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

%--------------------------------------------------------------------------
prevent_duplication_strt2=0;
maxsize_fst_coor= size(fast_total_coor);

for delete_strt_loop= 1:maxsize_fst_coor(1,2)
prevent_duplication_strt2(delete_strt_loop,1)=fast_total_coor{1,delete_strt_loop}(1,1);
prevent_duplication_strt2(delete_strt_loop,2)=fast_total_coor{1,delete_strt_loop}(1,2);
end

xyz_11=0;
comp_prevent_strt=0;
rsh_number=0;
    
xyz_11=xyz_1{1,fast_obj_start};
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

new_number_count=object_number_count;

for rest_start_searching= 1:max_frame_size

    fast_obj_start=2*rest_start_searching-1;

    if(fast_obj_start == 1)
        maxsize_fast_xyz=size(xyz_11); %read size
        change_inside_loop = maxsize_fast_xyz(1,1);
        start_number_skip=1;
    else
        change_inside_loop = new_number_count; 
        start_number_skip = object_number_count;
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
                if( estm_ft_size(1,1) >= fast_obj_start-1 )
                    
                    xx_1 = xyz{1,fast_obj_start+2}(fast_second_mt,1) - fast_total_coor{1,inside_loop}(fast_obj_start,1); % comparing x coordinates   
                    yy_1 = xyz{1,fast_obj_start+2}(fast_second_mt,2) - fast_total_coor{1,inside_loop}(fast_obj_start,2); % comparing y coordinates
                    
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
                    angle_fast_test = atan((xyz{1,fast_obj_start+2}(row_dist(try_loop),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))/(xyz{1,fast_obj_start+2}(row_dist(try_loop),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi; % unit is degree   
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
                        diff_vector = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1) - fast_total_coor{1,inside_loop}(fast_obj_start,1);
                        posit_nega=sign(diff_vector); % x_vector>0 : 1 , x_vector<0 : -1 , 0=0
                    else 
                        vector=0; % y vetorc
                        diff_vector = xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2) - fast_total_coor{1,inside_loop}(fast_obj_start,2); 
                        posit_nega=sign(diff_vector); % y_vector>0 : 1 , y_vector<0 : -1 , 0=0
                    end

                    diff_vector2=0;
                    posit_nega2=nan;  

                    if(vector == 1)         
                        diff_vector2= xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1);    
                        posit_nega2=sign(diff_vector2);
                    elseif (vector == 0)            
                        diff_vector2= xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2);
                        posit_nega2=sign(diff_vector2);          
                    end
                    
                    angle_fast_cmp=atan((xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))/(xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi;

                    if(abs(cont_angle_th(inside_loop,1))+control_vecter_angle*multilple_short_angle > abs(angle_fast_cmp) & abs(cont_angle_th(inside_loop,1))-control_vecter_angle*multilple_short_angle < abs(angle_fast_cmp)) % unit is degree
                        pass_numfil=1;
                    end                   
                end
                
                if(posit_nega == posit_nega2) 
                    pass_value=1;   
                    if(fast_obj_start == 1)                                
                        pass_value=0;   
                    else 
                        if((distance_fast_mat(row_dist(priority_angle_order(try_loop2,2)))+cont_rest_th(inside_loop,1))*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start-2,1))^2+(xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start-2,2))^2 ))    
                            pass_value=0; 
                        else 
                            pass_value=1; 
                        end 
                    end
                    if(pass_value == 0) 
                        if(pass_numfil == 1) 
                            if(fast_obj_start > 1)      
                                angle_fast = atan((xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))/(xyz{1,fast_obj_start+2}(row_dist(priority_angle_order(try_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1)))*180/pi;   
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
                                                if((distance_repeat_rest(row_secd(priority_angle_order2(find_secd_loop2,2)))+dist_thresh_sec)*increase_distance_threshold < sqrt( (xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),1)-fast_total_coor{1,inside_loop}(fast_obj_start,1))^2+(xyz{1,fast_obj_start+4}(row_secd(priority_angle_order2(find_secd_loop2,2)),2)-fast_total_coor{1,inside_loop}(fast_obj_start,2))^2 ))
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
ans=0;
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
imshow(abs(recon_stack(:,:,6)),[])
hold on
plot(xyz{1,6}(:,1),xyz{1,6}(:,2),'g.')
set(gca, 'XLim', [2250, 2350], 'YLim', [1100,1200])
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
set(gca, 'XLim', [2250, 2350], 'YLim', [1100,1200])
axis on

set(gca, 'XLim', [900,960], 'YLim', [0,60])
axis on

%-----------check
set(gca, 'XLim', [580,620], 'YLim', [1050,1150])
axis on

set(gca, 'XLim', [230,280], 'YLim', [300,360])
axis on

%next process increse for slow --> 400 % // incease rate 0.65 

set(gca, 'XLim', [2350,2550], 'YLim', [1360,1420])
axis on