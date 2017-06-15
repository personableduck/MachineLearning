copy_ct=coor_to;
copy_cf=coor_from;
number_max_matrix=0; %decide threshold

  for max_routine = 1:(iter-1) %To find maximum difference between coor_from and coor_to
  threshold_difference = max(abs(copy_ct{1,max_routine} - copy_cf{1,max_routine}));
  threshold_figure = max(threshold_difference);
  number_max_matrix=number_max_matrix+1;
  computation_threshold(number_max_matrix,1) = threshold_figure;
  end
  
positive_threshold = max(computation_threshold)+3; %threshold for rearrangement figures of coor_from and coor_to // figure 3 is about 30% of max difference.
negative_threshold=-(positive_threshold);

copy_ct_matrix_number=0;

for copy_ct_routine_num =1:(iter-2) %iter is images number
    copy_cf_routine_num = copy_ct_routine_num+1;
    
    copy_ct_sizef= size(copy_ct{1,copy_ct_routine_num}); %size calculate
    copy_ct_size= copy_ct_sizef(1,1);
    
    copy_cf_sizef= size(copy_cf{1,copy_cf_routine_num});
    copy_cf_size= copy_cf_sizef(1,1);
    
    for copy_cf_figure=1:copy_cf_size
        for copy_ct_figure=1:copy_ct_size

        if((negative_threshold <= (copy_cf{1,copy_cf_routine_num}(copy_cf_figure,1) - copy_ct{1,copy_ct_routine_num}(copy_ct_figure,1)) & (copy_cf{1,copy_cf_routine_num}(copy_cf_figure,1) - copy_ct{1,copy_ct_routine_num}(copy_ct_figure,1)) <= positive_threshold) & (negative_threshold <= (copy_cf{1,copy_cf_routine_num}(copy_cf_figure,2) - copy_ct{1,copy_ct_routine_num}(copy_ct_figure,2)) & (copy_cf{1,copy_cf_routine_num}(copy_cf_figure,2) - copy_ct{1,copy_ct_routine_num}(copy_ct_figure,2)) <= positive_threshold))
            
            copy_ct{1,copy_ct_routine_num}(copy_cf_figure,1)=0;
            copy_ct{1,copy_ct_routine_num}(copy_cf_figure,2)=0;
            
        end
        
        end
    end
    
    
    copy_change_matrix_num=0;
    copy_ct_sizef2= size(copy_ct{1,copy_ct_routine_num}); %size calculate
    copy_ct_size2= copy_ct_sizef2(1,1);
    
    copy_ct_size_next2= size(copy_ct{1,(copy_ct_routine_num+1)}); %size calculate
    copy_ct_size_next2= copy_ct_size_next2(1,1);
    
    copy_cf_sizef11= size(copy_cf{1,1}); %size calculate
    copy_cf_sizef11_add= copy_cf_sizef11(1,1);
    
    for copy_change_num = 1:copy_ct_size2
    
    if copy_ct{1,copy_ct_routine_num}(copy_change_num,1) > 0
        copy_change_matrix_num=copy_change_matrix_num+1;
        copy_ct{1,(copy_ct_routine_num+1)}(copy_ct_size_next2+copy_change_matrix_num,1)=copy_ct{1,copy_ct_routine_num}(copy_change_num,1);
        copy_ct{1,(copy_ct_routine_num+1)}(copy_ct_size_next2+copy_change_matrix_num,2)=copy_ct{1,copy_ct_routine_num}(copy_change_num,2);
        
        copy_cf{1,(copy_ct_routine_num+1)}(copy_ct_size_next2+copy_change_matrix_num,1)=copy_cf{1,copy_ct_routine_num}(copy_change_num,1);
        copy_cf{1,(copy_ct_routine_num+1)}(copy_ct_size_next2+copy_change_matrix_num,2)=copy_cf{1,copy_ct_routine_num}(copy_change_num,2);
    end
    
    end
end
            
            
            
            