%--------------------------------make starting point
xyz_1=0;
xyz_1=xyz;
for xyz_change_loop=1:iter-2
    
    comp_xyz_new=0;
    comp_xyz_new2=0;
    rsh_number=0;
    researching_special=0;

    comp_xyz_new=xyz{1,xyz_change_loop};
    comp_xyz_new2=xyz{1,xyz_change_loop+1};

    size_cmp_n=size(comp_xyz_new2);

    for loop_cmp_del=1:size_cmp_n(1,1)

        [row_del_x]=find(comp_xyz_new(:,1) == comp_xyz_new2(loop_cmp_del,1));
        [row_del_y]=find(comp_xyz_new(:,2) == comp_xyz_new2(loop_cmp_del,2));

        rr_ss=size(row_del_x);
        for check_loop_dlt=1:rr_ss
            [row_del_f,col_del_f,check_del_f]=find(row_del_y == row_del_x(check_loop_dlt));

            if(check_del_f ~= 0)
                rsh_number=rsh_number+1;
                researching_special(rsh_number,1)=xyz{1,xyz_change_loop}(row_del_x(check_loop_dlt),1);
                researching_special(rsh_number,2)=xyz{1,xyz_change_loop}(row_del_x(check_loop_dlt),2);
                researching_special(rsh_number,3)=row_del_x(check_loop_dlt);
            end
        end
    end

    maxsize_rsch_xyz1=size(researching_special); %read size

    for find_start_loop=1:maxsize_rsch_xyz1(1,1)-1
        [row_search_strt,col_search_strt,check_search_strt]=find( researching_special(:,1) < researching_special(find_start_loop,1) + distance_threshold_fast & researching_special(:,1) > researching_special(find_start_loop,1) - distance_threshold_fast); 
        maxsize_check_search_strt=size(check_search_strt);
        if(maxsize_check_search_strt(1,1) > 1)
            counting_strt=0;
            distance_strt=0;
            for row_search_strt_loop=1:maxsize_check_search_strt(1,1)
                srt_xx=researching_special((row_search_strt(row_search_strt_loop)),1)-researching_special(find_start_loop,1);
                srt_yy=researching_special((row_search_strt(row_search_strt_loop)),2)-researching_special(find_start_loop,2);
                if(srt_xx == 0)
                    srt_xx=nan;
                end
                if(srt_yy == 0)
                    srt_yy=nan;
                end
                counting_strt=counting_strt+1;
                distance_strt(counting_strt,1)=sqrt(srt_xx^2+srt_yy^2);
                distance_strt(counting_strt,2)=row_search_strt_loop;
            end
            [row_distrt,col_distrt,check_distrt]= find(distance_strt(:,1) <= distance_threshold_fast);
            if(check_distrt == 1)           
                left_only=nan;
                left_one_moving=0;
                left_one_moving(1,1)=min_intensities_xyz{1,xyz_change_loop}(researching_special(find_start_loop,3));
                left_one_moving(1,2)=researching_special(find_start_loop,3);
                maxsize_row_distrt=size(row_distrt);
                for row_distrt_loop=1:maxsize_row_distrt(1,1)
                    left_one_moving(row_distrt_loop+1,1)=min_intensities_xyz{1,xyz_change_loop}(researching_special(row_search_strt(distance_strt(row_distrt(row_distrt_loop),2)),3));
                    left_one_moving(row_distrt_loop+1,2)=researching_special(row_search_strt(distance_strt(row_distrt(row_distrt_loop),2)),3);
                end
                left_only=max(left_one_moving(:,1));
                maxsize_left_one_moving=size(left_one_moving);
                for left_one_moving_loop=1:maxsize_left_one_moving(1,1)
                    if(left_one_moving(left_one_moving_loop,1) ~= left_only)
                        xyz_1{1,xyz_change_loop}(left_one_moving(left_one_moving_loop,2),:) = nan;
                    end
                end
            else
                [row_search_strt2,col_search_strt2,check_search_strt2]=find( xyz{1,xyz_change_loop}(:,1) < xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1) + distance_threshold_fast & xyz{1,xyz_change_loop}(:,1) > xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1) - distance_threshold_fast);
                maxsize_check_search_strt2=size(check_search_strt2);
                if(maxsize_check_search_strt2(1,1) > 1)
                    counting_strt=0;
                    distance_strt=0;
                    check_distrt=0;
                    row_distrt=0;
                    col_distrt=0;
                    for row_search_strt_loop=1:maxsize_check_search_strt2(1,1)
                        srt_xx=xyz{1,xyz_change_loop}(row_search_strt2(row_search_strt_loop),1)-xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1);
                        srt_yy=xyz{1,xyz_change_loop}(row_search_strt2(row_search_strt_loop),2)-xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),2);
                        if(srt_xx == 0)
                            srt_xx=nan;
                        end
                        if(srt_yy == 0)
                            srt_yy=nan;
                        end
                        counting_strt=counting_strt+1;
                        distance_strt(counting_strt,1)=sqrt(srt_xx^2+srt_yy^2);
                        distance_strt(counting_strt,2)=row_search_strt_loop;
                    end
                    [row_distrt,col_distrt,check_distrt]= find(distance_strt(:,1) <= distance_threshold_fast);
                    if(check_distrt == 1)
                        xyz_1{1,xyz_change_loop}(researching_special(find_start_loop,3),:)=nan;
                    end
                end
            end
        else
            [row_search_strt3,col_search_strt3,check_search_strt3]=find( xyz{1,xyz_change_loop}(:,1) < xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1) + distance_threshold_fast & xyz{1,xyz_change_loop}(:,1) > xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1) - distance_threshold_fast);
            maxsize_check_search_strt3=size(check_search_strt3);
            if(maxsize_check_search_strt3(1,1) > 1)
                counting_strt=0;
                distance_strt=0;
                check_distrt=0;
                row_distrt=0;
                col_distrt=0;
                for row_search_strt_loop=1:maxsize_check_search_strt3(1,1)
                    srt_xx=xyz{1,xyz_change_loop}(row_search_strt3(row_search_strt_loop),1)-xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),1);
                    srt_yy=xyz{1,xyz_change_loop}(row_search_strt3(row_search_strt_loop),2)-xyz{1,xyz_change_loop}(researching_special(find_start_loop,3),2);
                    if(srt_xx == 0)
                        srt_xx=nan;
                    end
                    if(srt_yy == 0)
                        srt_yy=nan;
                    end
                    counting_strt=counting_strt+1;
                    distance_strt(counting_strt,1)=sqrt(srt_xx^2+srt_yy^2);
                    distance_strt(counting_strt,2)=row_search_strt_loop;
                end
                [row_distrt,col_distrt,check_distrt]= find(distance_strt(:,1) <= distance_threshold_fast);
                if(check_distrt == 1)
                    xyz_1{1,xyz_change_loop}(researching_special(find_start_loop,3),:)=nan;
                end
            end
        end
    end
end