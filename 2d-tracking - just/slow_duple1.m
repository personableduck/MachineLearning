function [fast_total_coor,new_search1,stop_search,new_number_count]=slow_duple1(xyz_11,fst_mv_oj_coor,fast_total_coor,new_number_count,stop_search,angle_information,object_number_count,fast_obj_start)
new_search1=nan;
duplication_number=0;
prevent_check_size=size(fast_total_coor);
%-----prevent duplication----------
if( prevent_check_size(1,2) == object_number_count)
    fast_total_coor(1,new_number_count)={fst_mv_oj_coor};
    stop_search=1; %prevent duplicate trajectories
else
    %------prevent duplication--------
    maxsize_fstc=size(fast_total_coor);     
    maxsize_cmp_fst=size(fst_mv_oj_coor);
    duplehappen=0;
    for prevent_duplication_loop=object_number_count+1:maxsize_fstc(1,2)	
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
        if(duplication_happen >= 1) 
            duplication_number=prevent_duplication_loop;
        end
    end
    if(duplication_number ~= 0)                                                                                  
        if(angle_information(duplication_number,1) < angle_information(new_number_count,1))
            duplehappen=1;
        end
        if( duplehappen == 1)
            new_number_count= new_number_count-1; %duplication happen
        else
            [row_srh_x,col_srh_x,check_srh_x]=find(xyz_11{1,fast_obj_start}(:,1) == fast_total_coor{1,duplication_number}(fast_obj_start,1));
            [row_srh_y,col_srh_y,check_srh_y]=find(xyz_11{1,fast_obj_start}(:,2) == fast_total_coor{1,duplication_number}(fast_obj_start,2));
            rr_srh=length(row_srh_x);
            for check_loop_srh=1:rr_srh
                [row_srh_f,col_srh_f,check_srh_f]=find(row_srh_y == row_srh_x(check_loop_srh));
                if(check_srh_f ~= 0)
                    new_search1=row_srh_x(check_loop_srh);
                end
            end    
            
            fast_total_coor(1,duplication_number) = {0};
            new_number_count = new_number_count-1;
            maxsize_fst_check=size(fst_mv_oj_coor);
            for thirdnumberchange=1:maxsize_fst_check(1,1)
                fst_mv_oj_coor(thirdnumberchange,3) = duplication_number;
            end
            fast_total_coor(1,duplication_number) = {fst_mv_oj_coor}; %smaller angle is right trajectory
            stop_search=1;
        end
    else
        fast_total_coor(1,new_number_count)={fst_mv_oj_coor};
        stop_search=1; %prevent duplicate trajectories
    end
end


end