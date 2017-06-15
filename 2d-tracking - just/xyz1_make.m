comp_xyz_new=0;
comp_xyz_new2=0;
rsh_number=0;
    
comp_xyz_new=xyz{1,1};
comp_xyz_new2=xyz{1,2};

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
maxsize_behinding=size(xyz{1,1});

count_change_num=0;
change_nan=0;

for find_start_loop=1:maxsize_rsch_xyz1(1,1)
    
    distance_strt=0;
    
    for behinding_loop= 1:maxsize_behinding(1,1)
        
        srt_xx=xyz{1,1}(behinding_loop,1)-xyz{1,1}(researching_xyz1(find_start_loop,1),1);
        srt_yy=xyz{1,1}(behinding_loop,2)-xyz{1,1}(researching_xyz1(find_start_loop,1),2);
        
        distance_strt(behinding_loop,1)=sqrt(srt_xx^2+srt_yy^2);
        
    end

    [row_strt,col_strt,check_strt]=find(distance_strt <= distance_threshold_fast & distance_strt ~= 0);
    
    if(check_strt == 1)
        
        row_strt_size=size(row_strt); %read size 
        
        for row_strt_loop=1:row_strt_size(1,1)
            if(min_intensities_xyz{1,1}(researching_xyz1(find_start_loop,1)) < min_intensities_xyz{1,1}(row_strt(row_strt_loop)) )
               count_change_num=count_change_num+1;
               change_nan(count_change_num,1)=researching_xyz1(find_start_loop,1);
            end
        end
    end
end

change_nan_size=size(change_nan); %read size

xyz_1=xyz{1,1};

for nan_put1= 1: change_nan_size(1,1)
    xyz_1(change_nan(nan_put1,1),:)=nan;
end



