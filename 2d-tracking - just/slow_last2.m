function [fast_total_coor,stop_search,prevent_rest_check]=slow_last2(origin_size_fst,fst_mv_oj_coor,fast_total_coor,prevent_rest_check,duple_mat,ans_xyz3,jump_number)

%------prevent duplication--------
maxsize_fstc=size(prevent_rest_check);     
maxsize_cmp_fst=size(fst_mv_oj_coor);
duplication_number=0;
duplecount=0;
duplehappen=0;
prevent_duplication_loop=0;
for prevent_duplication_loop=1:maxsize_fstc(1,2)	
    duplication_happen=0;
    for fst_check_loop =1:maxsize_cmp_fst(1,1)
        [row_dplc_x,col_dplc_x,check_dplc_x]=find(prevent_rest_check{1,prevent_duplication_loop}(:,1) == fst_mv_oj_coor(fst_check_loop,1));
        [row_dplc_y,col_dplc_y,check_dplc_y]=find(prevent_rest_check{1,prevent_duplication_loop}(:,2) == fst_mv_oj_coor(fst_check_loop,2));
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

if(duplication_number == 0)
    prevent_rest_check{1,duple_mat}(1,1) = fast_total_coor{1,duple_mat}(origin_size_fst,1);
    prevent_rest_check{1,duple_mat}(1,2) = fast_total_coor{1,duple_mat}(origin_size_fst,2);
    prevent_rest_check{1,duple_mat}(1,3) = fast_total_coor{1,duple_mat}(origin_size_fst,3);
    prevent_rest_check{1,duple_mat}(2,1) = fst_mv_oj_coor(1,1);
    prevent_rest_check{1,duple_mat}(2,2) = fst_mv_oj_coor(1,2);
    prevent_rest_check{1,duple_mat}(2,3) = fst_mv_oj_coor(1,3); 
    prevent_rest_check{1,duple_mat}(3,1) = fst_mv_oj_coor(2,1);
    prevent_rest_check{1,duple_mat}(3,2) = fst_mv_oj_coor(2,2);
    prevent_rest_check{1,duple_mat}(3,3) = fst_mv_oj_coor(2,3);

    fast_total_coor{1,duple_mat}(origin_size_fst+1,1) = fst_mv_oj_coor(1,1);
    fast_total_coor{1,duple_mat}(origin_size_fst+1,2) = fst_mv_oj_coor(1,2);       
    fast_total_coor{1,duple_mat}(origin_size_fst+1,3) = duple_mat;
    fast_total_coor{1,duple_mat}(origin_size_fst+1,4) = jump_number;

    fast_total_coor{1,duple_mat}(origin_size_fst+2,1) = fst_mv_oj_coor(2,1);
    fast_total_coor{1,duple_mat}(origin_size_fst+2,2) = fst_mv_oj_coor(2,2);         
    fast_total_coor{1,duple_mat}(origin_size_fst+2,3) = duple_mat;
    fast_total_coor{1,duple_mat}(origin_size_fst+2,4) = jump_number;

    if(~isnan(ans_xyz3))
        prevent_rest_check{1,duple_mat}(4,1)= fst_mv_oj_coor(3,1);           
        prevent_rest_check{1,duple_mat}(4,2)= fst_mv_oj_coor(3,2);
        prevent_rest_check{1,duple_mat}(4,3)= fst_mv_oj_coor(3,3);

        fast_total_coor{1,duple_mat}(origin_size_fst+3,1) = fst_mv_oj_coor(3,1);
        fast_total_coor{1,duple_mat}(origin_size_fst+3,2) = fst_mv_oj_coor(3,2);
        fast_total_coor{1,duple_mat}(origin_size_fst+3,3) = duple_mat;
        fast_total_coor{1,duple_mat}(origin_size_fst+3,4) = jump_number;
    end
    stop_search=1; %prevent duplicate trajectories
end


end