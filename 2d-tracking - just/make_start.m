%--------------------------------make starting point


xyz_1=0;
xyz_1=xyz;

for change_xyz_loop=1:iter-2

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