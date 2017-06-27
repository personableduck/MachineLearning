function [fast_total_coor,stop_search,new_number_count]=slow_last(fst_mv_oj_coor,fast_total_coor,new_number_count,stop_search)

 %-----prevent duplication----
prevent_check_size=size(fast_total_coor);
%------prevent duplication--------
maxsize_fstc=size(fast_total_coor);     
maxsize_cmp_fst=size(fst_mv_oj_coor);
duplecount=0;
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
    if(duplication_happen >= 1) 
        duplecount=duplecount+1;
    end
end

if(duplecount == 0)
    fast_total_coor(1,new_number_count)={fst_mv_oj_coor};
    stop_search=1; %prevent duplicate trajectories
end

end
