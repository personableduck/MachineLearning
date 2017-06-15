%New trajectories3 by DUCK-HA HWANG

%------solve diff holo dislocation, change xyz coordinate information.

threshold_diff_dislocation= 1.3; %*****check it!! important value (1+1*30%) for combining close coordinates

fps=3; %fps means frame per second
distance_threshold_fast = (120/fps)/pixelsize; %******threshold object tracking. need to control. maximum ditance to distinguish as a same object
%unit is um. // I limited that the fasted sperm's speed is 110 um/sec.
%fps means frame per second

for combine_compare_loop=1:iter-2
    
    behind_compare=combine_compare_loop+1;
    
    maxsize_cmb_cp=size(xyz{1,combine_compare_loop}); %read size
    
    maxsize_cmb_cp2=size(xyz{1,behind_compare}); %read size

    
    for inft_loop=1:maxsize_cmb_cp(1,1)
        for bhid_loop=1:maxsize_cmb_cp2(1,1)

            xx_shrt = xyz{1,combine_compare_loop}(inft_loop,1) - xyz{1,behind_compare}(bhid_loop,1); % comparing x coordinates
            yy_shrt = xyz{1,combine_compare_loop}(inft_loop,2) - xyz{1,behind_compare}(bhid_loop,2); % comparing y coordinates
        
            distance_shrt_limit= sqrt((xx_shrt^2)+(yy_shrt^2));

            if (distance_shrt_limit <= threshold_diff_dislocation)         
                xyz{1,behind_compare}(bhid_loop,1) = xyz{1,combine_compare_loop}(inft_loop,1);
                xyz{1,behind_compare}(bhid_loop,2) = xyz{1,combine_compare_loop}(inft_loop,2);
                
            end
        end
    end
    
end

%------change xyz

object_number_count=0;

fast_obj_start=1; % if you need to find starting point on the middle of frame, then you need to change the value as a loop(for function)
maxsize_fast_xyz=size(xyz{1,fast_obj_start}); %read size  
    
   for inside_loop=1:maxsize_fast_xyz(1,1)-1
       
       count_distance_num=0;
       distance_fast_mat=0;
       
       last_distance_comp={0};
       find_alaram=0; %check for finding
       vector_info_fast=0;
       
       for inside_bhd=inside_loop+1:maxsize_fast_xyz(1,1)
           
           count_distance_num=count_distance_num+1;
    
            xx_1 = xyz{1,fast_obj_start}(inside_bhd,1) - xyz{1,fast_obj_start}(inside_loop,1); % comparing x coordinates
            yy_1 = xyz{1,fast_obj_start}(inside_bhd,2) - xyz{1,fast_obj_start}(inside_loop,2); % comparing y coordinates
    
            distance_fast_mat(count_distance_num,1) = sqrt((xx_1^2)+(yy_1^2));
            
       end
       
       last_distance_comp(1,inside_loop)={distance_fast_mat};
       
       [temp,ord] = sort(last_distance_comp{1,inside_loop}(:,1)); %rearrangement for distance from the smallest to the largest
       find_distance_order = last_distance_comp{1,inside_loop}(ord,:);
           
       %******************************************need to improve(check
       %threschold first
       
       for try_time_loop = 1:5;
           %add loop here 5 time 
                  
           cont_third_chk=0;
           result_value=0;
           
            %---size problem after inside_loop be max       
                  if(inside_loop > maxsize_fast_xyz(1,1)-5)
                      try_time_loop=1;           
                  end
           
           [row_fdo]=find(last_distance_comp{1,inside_loop} == find_distance_order(try_time_loop,1),1);      
           
           distance_fast=last_distance_comp{1,inside_loop}(row_fdo);
      
            if(distance_fast <= distance_threshold_fast)                
                                      
                angle_fast=asin((xyz{1,fast_obj_start}(inside_loop+row_fdo,2)-xyz{1,fast_obj_start}(inside_loop,2))/distance_fast); %unit is Radians(rad)
                increase_angle=angle_fast+(pi/180*45); % ****important value for erro (sum 45 dgree)  
                decrease_angle=angle_fast-(pi/180*45); % ****important value for erro (minus 45 dgree)
                
                if (min_intensities_xyz{1,fast_obj_start}(inside_loop) > min_intensities_xyz{1,fast_obj_start}(inside_loop+row_fdo)) % compare intensity to decide before one.         
                 
                    if(xyz{1,fast_obj_start}(inside_loop,1) > xyz{1,fast_obj_start}(inside_loop+row_fdo,1)) % distinguish for 1,2 quaters 
                        
                        angle_fast=(180*pi/180)-angle_fast;
                        increase_angle=angle_fast+(pi/180*45); % ****important value for erro (sum 45 dgree)  
                        decrease_angle=angle_fast-(pi/180*45); % ****important value for erro (minus 45 dgree)     

                        vector_info_fast=2;
                        
                        [row_ang,col_ang,check_ang]=find(xyz{1,fast_obj_start+1} < (xyz{1,fast_obj_start}(inside_loop+row_fdo,1)+(cos(decrease_angle)*(distance_fast*0.8))) & xyz{1,fast_obj_start+1} > (xyz{1,fast_obj_start}(inside_loop+row_fdo,1)+(cos(increase_angle)*(distance_fast*1.2))),1); %find x coordinate
                        if (check_ang == 1)                          
                            if ( xyz{1,fast_obj_start+1}(row_ang,2) > (xyz{1,fast_obj_start}(inside_loop+row_fdo,2)+(sin(increase_angle)*(distance_fast*0.8))) & xyz{1,fast_obj_start+1}(row_ang,2) < (xyz{1,fast_obj_start}(inside_loop+row_fdo,2)+(sin(decrease_angle)*(distance_fast*1.2))))          
                                % my second moving coordinats: xyz{1,fast_obj_loop+1}(row_ang,1),xyz{1,fast_obj_loop+1}(row_ang,2)
                                %------check two frame
                                %improve accuracy check three frame
                                
                                cont_third_chk=1;
                            end
                        end
                                
                    else
                        
                        vector_info_fast=1;
                        
                        [row_ang,col_ang,check_ang]=find(xyz{1,fast_obj_start+1} < (xyz{1,fast_obj_start}(inside_loop+row_fdo,1)+(cos(decrease_angle)*(distance_fast*1.2))) & xyz{1,fast_obj_start+1} > (xyz{1,fast_obj_start}(inside_loop+row_fdo,1)+(cos(increase_angle)*(distance_fast*0.8))),1); %find x coordinate
                        
                        if (check_ang == 1)                          
                            if ( xyz{1,fast_obj_start+1}(row_ang,2) > (xyz{1,fast_obj_start}(inside_loop+row_fdo,2)+(sin(decrease_angle)*(distance_fast*0.8))) & xyz{1,fast_obj_start+1}(row_ang,2) < (xyz{1,fast_obj_start}(inside_loop+row_fdo,2)+(sin(increase_angle)*(distance_fast*1.2))))          
                                % my second moving coordinats: xyz{1,fast_obj_loop+1}(row_ang,1),xyz{1,fast_obj_loop+1}(row_ang,2)
                                %------check two frame
                                %improve accuracy check three frame
                                
                                cont_third_chk=1;
                            end
                        end
                                
                    end
                    
                           if(cont_third_chk ==1) 
                               
                                if(xyz{1,fast_obj_start}(inside_loop,1) > xyz{1,fast_obj_start}(inside_loop+row_fdo,1)) % distinguish for 1,2 quaters                        
                                    [row_ang_3,col_ang_3,check_ang_3]= find(xyz{1,fast_obj_start+2} > (xyz{1,fast_obj_start+1}(row_ang,1)+(cos(increase_angle)*(distance_fast*1.2))) & xyz{1,fast_obj_start+2} < (xyz{1,fast_obj_start+1}(row_ang,1)+(cos(decrease_angle)*(distance_fast*0.8))),1); % find x coordinate
                                    if (check_ang_3 == 1)
                                        if ( xyz{1,fast_obj_start+2}(row_ang_3,2) < (xyz{1,fast_obj_start+1}(row_ang,2)+(sin(decrease_angle)*(distance_fast*1.2))) & xyz{1,fast_obj_start+2}(row_ang_3,2) > (xyz{1,fast_obj_start+1}(row_ang,2)+(sin(increase_angle)*(distance_fast*0.8))))
                                            result_value=1;
                                        end
                                    end
                                else
                                    [row_ang_3,col_ang_3,check_ang_3]= find(xyz{1,fast_obj_start+2} < (xyz{1,fast_obj_start+1}(row_ang,1)+(cos(decrease_angle)*(distance_fast*1.2))) & xyz{1,fast_obj_start+2} > (xyz{1,fast_obj_start+1}(row_ang,1)+(cos(increase_angle)*(distance_fast*0.8))),1); % find x coordinate
                                    if (check_ang_3 == 1)
                                        if ( xyz{1,fast_obj_start+2}(row_ang_3,2) > (xyz{1,fast_obj_start+1}(row_ang,2)+(sin(decrease_angle)*(distance_fast*0.8))) & xyz{1,fast_obj_start+2}(row_ang_3,2) < (xyz{1,fast_obj_start+1}(row_ang,2)+(sin(increase_angle)*(distance_fast*1.2))))
                                            result_value=1;
                                        end
                                    end
                                end
  
                                if(result_value ==1)
                                
                                    %my third moving object:
                                    %xyz{1,fast_obj_loop+2}(row_ang_3,1),xyz{1,fast_obj_loop+2}(row_ang_3,2)
                                    object_number_count = object_number_count+1;
                                   
                                    fst_mv_oj_coor(1,1) = xyz{1,fast_obj_start}(inside_loop,1);                                
                                    fst_mv_oj_coor(1,2) = xyz{1,fast_obj_start}(inside_loop,2);            
                                    fst_mv_oj_coor(1,3) = 2; %fast moving value is 2
                                    fst_mv_oj_coor(1,4) = object_number_count;
                                   
                                    fst_mv_oj_coor(2,1) = xyz{1,fast_obj_start+1}(row_ang,1);
                                    fst_mv_oj_coor(2,2) = xyz{1,fast_obj_start+1}(row_ang,2);           
                                    fst_mv_oj_coor(2,3) = 2; %fast moving value is 2
                                    fst_mv_oj_coor(2,4) = object_number_count;
                                   
                                    fst_mv_oj_coor(3,1) = xyz{1,fast_obj_start+2}(row_ang_3,1);
                                    fst_mv_oj_coor(3,2) = xyz{1,fast_obj_start+2}(row_ang_3,2);          
                                    fst_mv_oj_coor(3,3) = 2; %fast moving value is 2
                                    fst_mv_oj_coor(3,4) = object_number_count;
                                   
                                    find_alaram=1; %check for finding
                                    continue   
                                end
                           end
        
                else % for intensity opposite direction
                    
                    angle_fast=-angle_fast;
                    
                    if(xyz{1,fast_obj_start}(inside_loop,1) > xyz{1,fast_obj_start}(inside_loop+row_fdo,1)) % distinguish for 1,2 quaters 
                        
                        angle_fast=-(180*pi/180)-angle_fast;
                        increase_angle=angle_fast+(pi/180*45); % ****important value for erro (sum 45 dgree)  
                        decrease_angle=angle_fast-(pi/180*45); % ****important value for erro (minus 45 dgree)     

                        vector_info_fast=4;
                        [row_ang,col_ang,check_ang]=find(xyz{1,fast_obj_start+1} < (xyz{1,fast_obj_start}(inside_loop,1)+(cos(increase_angle)*(distance_fast*1.2))) & xyz{1,fast_obj_start+1} > (xyz{1,fast_obj_start}(inside_loop,1)+(cos(decrease_angle)*(distance_fast*0.8))),1); % find x coordinate
                    
                        if (check_ang == 1)                          
                            if ( xyz{1,fast_obj_start+1}(row_ang,2) > (xyz{1,fast_obj_start}(inside_loop,2)+(sin(decrease_angle)*(distance_fast*1.2))) & xyz{1,fast_obj_start+1}(row_ang,2) < (xyz{1,fast_obj_start}(inside_loop,2)+(sin(increase_angle)*(distance_fast*0.8))))
                         
                                cont_third_chk=1;
                            end
                        end
                        
                        
                    else
                        vector_info_fast=3;
                        [row_ang,col_ang,check_ang]=find(xyz{1,fast_obj_start+1} > (xyz{1,fast_obj_start}(inside_loop,1)+(cos(decrease_angle)*(distance_fast*1.2))) & xyz{1,fast_obj_start+1} < (xyz{1,fast_obj_start}(inside_loop,1)+(cos(increase_angle)*(distance_fast*0.8))),1); % find x coordinate
                  
                        if (check_ang == 1)                          
                            if ( xyz{1,fast_obj_start+1}(row_ang,2) < (xyz{1,fast_obj_start}(inside_loop,2)+(sin(decrease_angle)*(distance_fast*0.8))) & xyz{1,fast_obj_start+1}(row_ang,2) > (xyz{1,fast_obj_start}(inside_loop,2)+(sin(increase_angle)*(distance_fast*1.2))))
                         
                                cont_third_chk=1;
                            end
                        end
                        
                    end
                    
                    %%%%%
                    if(cont_third_chk ==1) 

                        if(xyz{1,fast_obj_start}(inside_loop,1) > xyz{1,fast_obj_start}(inside_loop+row_fdo,1)) % distinguish for 1,2 quaters                        
                            [row_ang_3,col_ang_3,check_ang_3]= find(xyz{1,fast_obj_start+2} > (xyz{1,fast_obj_start+1}(row_ang,1)+(cos(decrease_angle)*(distance_fast*0.8))) & xyz{1,fast_obj_start+2} < (xyz{1,fast_obj_start+1}(row_ang,1)+(cos(increase_angle)*(distance_fast*1.2))),1); % find x coordinate
                            if (check_ang_3 == 1)
                                if ( xyz{1,fast_obj_start+2}(row_ang_3,2) > (xyz{1,fast_obj_start+1}(row_ang,2)+(sin(decrease_angle)*(distance_fast*1.2))) & xyz{1,fast_obj_start+2}(row_ang_3,2) < (xyz{1,fast_obj_start+1}(row_ang,2)+(sin(increase_angle)*(distance_fast*0.8))))
                                    result_value=1;     
                                end
                            end
                        else
                            [row_ang_3,col_ang_3,check_ang_3]= find(xyz{1,fast_obj_start+2} > (xyz{1,fast_obj_start+1}(row_ang,1)+(cos(decrease_angle)*(distance_fast*1.2))) & xyz{1,fast_obj_start+2} < (xyz{1,fast_obj_start+1}(row_ang,1)+(cos(increase_angle)*(distance_fast*0.8))),1); % find x coordinate
                            if (check_ang_3 == 1)
                                if ( xyz{1,fast_obj_start+2}(row_ang_3,2) < (xyz{1,fast_obj_start+1}(row_ang,2)+(sin(decrease_angle)*(distance_fast*0.8))) & xyz{1,fast_obj_start+2}(row_ang_3,2) > (xyz{1,fast_obj_start+1}(row_ang,2)+(sin(increase_angle)*(distance_fast*1.2))))                                                  
                                    result_value=1;     
                                end
                            end
                        end
                        
                        if(result_value ==1)
     
                            object_number_count = object_number_count+1;
                                   
                            fst_mv_oj_coor(1,1) = xyz{1,fast_obj_start}(inside_loop,1);
                            fst_mv_oj_coor(1,2) = xyz{1,fast_obj_start}(inside_loop,2);            
                            fst_mv_oj_coor(1,3) = 2; %fast moving value is 2
                            fst_mv_oj_coor(1,4) = object_number_count;
                                   
                            fst_mv_oj_coor(2,1) = xyz{1,fast_obj_start+1}(row_ang,1);
                            fst_mv_oj_coor(2,2) = xyz{1,fast_obj_start+1}(row_ang,2);           
                            fst_mv_oj_coor(2,3) = 2; %fast moving value is 2
                            fst_mv_oj_coor(2,4) = object_number_count;
                                   
                            fst_mv_oj_coor(3,1) = xyz{1,fast_obj_start+2}(row_ang_3,1);
                            fst_mv_oj_coor(3,2) = xyz{1,fast_obj_start+2}(row_ang_3,2);          
                            fst_mv_oj_coor(3,3) = 2; %fast moving value is 2
                            fst_mv_oj_coor(3,4) = object_number_count;
                                   
                            find_alaram=1; %check for finding     
                            continue   
                        end
                    end
                end% end for intensity
            end%end for distance Threshold  
      
       end  %end for try loop
            
       if (find_alaram == 0)           
           continue
       end
           
%------- find 1 object
%----- we already get 1,2,3 frame moving object
       
       for rest_tj_loop=1:(iter-4)%keep going to find rest of trajectories
           
           rest_frame_num = rest_tj_loop+3;
           
           if( vector_info_fast== 1)
  
               [row_ang_rest,col_ang_rest,check_ang_rest]= find(xyz{1,rest_frame_num} < (fst_mv_oj_coor(3,1)+(cos(decrease_angle)*(distance_fast*1.2))) & xyz{1,rest_frame_num} > (fst_mv_oj_coor(3,1)+(cos(increase_angle)*(distance_fast*0.8))),1); % find x coordinate     
                    if (check_ang_rest == 1)                          
                        if ( fst_mv_oj_coor(3,2) > (fst_mv_oj_coor(3,2)+(sin(decrease_angle)*(distance_fast*0.8))) & fst_mv_oj_coor(3,2) < (fst_mv_oj_coor(3,2)+(sin(increase_angle)*(distance_fast*1.2))))   
                            
                            fst_mv_oj_coor(rest_frame_num,1)=xyz{1,rest_frame_num}(row_ang_rest,1);
                            fst_mv_oj_coor(rest_frame_num,2)=xyz{1,rest_frame_num}(row_ang_rest,2);
                            fst_mv_oj_coor(rest_frame_num,3) = 2; %fast moving value is 2
                            fst_mv_oj_coor(rest_frame_num,4) = object_number_count;
                            
                            %result of finding rest trajectories      
                        end
                    end
           end
     
           if( vector_info_fast== 2)
   
               [row_ang_rest,col_ang_rest,check_ang_rest]= find(xyz{1,rest_frame_num} < (fst_mv_oj_coor(3,1)+(cos(decrease_angle)*(distance_fast*1.8))) & xyz{1,rest_frame_num} > (fst_mv_oj_coor(3,1)+(cos(increase_angle)*(distance_fast*1.2))),1); % find x coordinate     
               
               if (check_ang_rest == 1)                            
                   if ( fst_mv_oj_coor(3,2) > (fst_mv_oj_coor(3,2)+(sin(increase_angle)*(distance_fast*0.8))) & fst_mv_oj_coor(3,2) < (fst_mv_oj_coor(3,2)+(sin(decrease_angle)*(distance_fast*1.2))))
                             
                       fst_mv_oj_coor(rest_frame_num,1)=xyz{1,rest_frame_num}(row_ang_rest,1);   
                       fst_mv_oj_coor(rest_frame_num,2)=xyz{1,rest_frame_num}(row_ang_rest,2);
                       fst_mv_oj_coor(rest_frame_num,3) = 2; %fast moving value is 2
                       fst_mv_oj_coor(rest_frame_num,4) = object_number_count;
                       %result of finding rest trajectories      
                        
                   end
               end
               
           end
              
           if( vector_info_fast== 3)  
                        
               [row_ang_rest,col_ang_rest,check_ang_rest]= find(xyz{1,rest_frame_num} > (fst_mv_oj_coor(3,1)+(cos(decrease_angle)*(distance_fast*1.2))) & xyz{1,rest_frame_num} < (fst_mv_oj_coor(3,1)+(cos(increase_angle)*(distance_fast*0.8))),1); % find x coordinate     
               
               if (check_ang_rest == 1)                            
                   if ( fst_mv_oj_coor(3,2) < (fst_mv_oj_coor(3,2)+(sin(decrease_angle)*(distance_fast*0.8))) & fst_mv_oj_coor(3,2) > (fst_mv_oj_coor(3,2)+(sin(increase_angle)*(distance_fast*1.2))))
                              
                       fst_mv_oj_coor(rest_frame_num,1)=xyz{1,rest_frame_num}(row_ang_rest,1);   
                       fst_mv_oj_coor(rest_frame_num,2)=xyz{1,rest_frame_num}(row_ang_rest,2);   
                       fst_mv_oj_coor(rest_frame_num,3) = 2; %fast moving value is 2
                       fst_mv_oj_coor(rest_frame_num,4) = object_number_count;
                       %result of finding rest trajectories      
                        
                   end
               end
               
           end

           if( vector_info_fast== 4)
               
               [row_ang_rest,col_ang_rest,check_ang_rest]= find(xyz{1,rest_frame_num} < (fst_mv_oj_coor(3,1)+(cos(increase_angle)*(distance_fast*1.2))) & xyz{1,rest_frame_num} > (fst_mv_oj_coor(3,1)+(cos(decrease_angle)*(distance_fast*0.8))),1); % find x coordinate     
               
               if (check_ang_rest == 1)                            
                   if ( fst_mv_oj_coor(3,2) > (fst_mv_oj_coor(3,2)+(sin(decrease_angle)*(distance_fast*1.2))) & fst_mv_oj_coor(3,2) < (fst_mv_oj_coor(3,2)+(sin(increase_angle)*(distance_fast*0.8))))
                       
                       fst_mv_oj_coor(rest_frame_num,1)=xyz{1,rest_frame_num}(row_ang_rest,1);   
                       fst_mv_oj_coor(rest_frame_num,2)=xyz{1,rest_frame_num}(row_ang_rest,2);   
                       fst_mv_oj_coor(rest_frame_num,3) = 2; %fast moving value is 2
                       fst_mv_oj_coor(rest_frame_num,4) = object_number_count;
                       %result of finding rest trajectories      
                        
                   end
               end
               
           end
               
               
              
       end%end for rest loop

       fast_trj_coordinates(1,object_number_count)={fst_mv_oj_coor}; %actual coordinate information
                 
   end %end for distance
                

                    
                   
                               
                               
% figure, imshow(abs(recon_stack(:,:,1)),[])
% for test_image = 1:iter-1
%     hold on
%     plot(xyz{1,test_image}(:,1),xyz{1,test_image}(:,2),'g.')
% end
%     
%          
% set(gca, 'XLim', [0, 15], 'YLim', [880, 910])
% axis on
                        
          
 figure, imshow(recon_stack_diff(:,:,1),[])
 hold on
 plot(xyz{1,1}(:,1),xyz{1,1}(:,2),'g.')
 
set(gca, 'XLim', [0, 15], 'YLim', [880, 910])
axis on


% sin(30*pi/180)
% asin(1/2)*180/pi
