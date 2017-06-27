function [fast_total_coor,prevent_rest_check,duple_mat,stop_search]=prevent_duple2(fst_mv_oj_coor,fast_total_coor,inside_loop,prevent_rest_check,ans_xyz5,origin_size_fst,jump_number,stop_search,angle_information)
%-----prevent duplication----
duplication_number=0;
duple_mat=nan;

if( prevent_rest_check{1,1} == 0)

    prevent_rest_check{1,inside_loop}(1,1) = fast_total_coor{1,inside_loop}(origin_size_fst,1);
    prevent_rest_check{1,inside_loop}(1,2) = fast_total_coor{1,inside_loop}(origin_size_fst,2);
    prevent_rest_check{1,inside_loop}(1,3) = fast_total_coor{1,inside_loop}(origin_size_fst,3);
    prevent_rest_check{1,inside_loop}(2,1) = fst_mv_oj_coor(1,1);
    prevent_rest_check{1,inside_loop}(2,2) = fst_mv_oj_coor(1,2);
    prevent_rest_check{1,inside_loop}(2,3) = fst_mv_oj_coor(1,3); 
    prevent_rest_check{1,inside_loop}(3,1) = fst_mv_oj_coor(2,1);
    prevent_rest_check{1,inside_loop}(3,2) = fst_mv_oj_coor(2,2); 
    prevent_rest_check{1,inside_loop}(3,3) = fst_mv_oj_coor(2,3);
    prevent_rest_check{1,inside_loop}(4,1) = fst_mv_oj_coor(3,1);
    prevent_rest_check{1,inside_loop}(4,2) = fst_mv_oj_coor(3,2);
    prevent_rest_check{1,inside_loop}(4,3) = fst_mv_oj_coor(3,3);
    prevent_rest_check{1,inside_loop}(5,1) = fst_mv_oj_coor(4,1);
    prevent_rest_check{1,inside_loop}(5,2) = fst_mv_oj_coor(4,2);
    prevent_rest_check{1,inside_loop}(5,3) = fst_mv_oj_coor(4,3);

    fast_total_coor{1,inside_loop}(origin_size_fst+1,1) = fst_mv_oj_coor(1,1);
    fast_total_coor{1,inside_loop}(origin_size_fst+1,2) = fst_mv_oj_coor(1,2);       
    fast_total_coor{1,inside_loop}(origin_size_fst+1,3) = inside_loop;
    fast_total_coor{1,inside_loop}(origin_size_fst+1,4) = jump_number;

    fast_total_coor{1,inside_loop}(origin_size_fst+2,1) = fst_mv_oj_coor(2,1);
    fast_total_coor{1,inside_loop}(origin_size_fst+2,2) = fst_mv_oj_coor(2,2);         
    fast_total_coor{1,inside_loop}(origin_size_fst+2,3) = inside_loop;
    fast_total_coor{1,inside_loop}(origin_size_fst+2,4) = jump_number;

    fast_total_coor{1,inside_loop}(origin_size_fst+3,1) = fst_mv_oj_coor(3,1);
    fast_total_coor{1,inside_loop}(origin_size_fst+3,2) = fst_mv_oj_coor(3,2);       
    fast_total_coor{1,inside_loop}(origin_size_fst+3,3) = inside_loop;
    fast_total_coor{1,inside_loop}(origin_size_fst+3,4) = jump_number;

    fast_total_coor{1,inside_loop}(origin_size_fst+4,1) = fst_mv_oj_coor(4,1);
    fast_total_coor{1,inside_loop}(origin_size_fst+4,2) = fst_mv_oj_coor(4,2);
    fast_total_coor{1,inside_loop}(origin_size_fst+4,3) = inside_loop; 
    fast_total_coor{1,inside_loop}(origin_size_fst+4,4) = jump_number;

    if(~isnan(ans_xyz5))
        fast_total_coor{1,inside_loop}(origin_size_fst+5,1) = fst_mv_oj_coor(5,1);
        fast_total_coor{1,inside_loop}(origin_size_fst+5,2) = fst_mv_oj_coor(5,2);
        fast_total_coor{1,inside_loop}(origin_size_fst+5,3) = inside_loop;
        fast_total_coor{1,inside_loop}(origin_size_fst+5,4) = jump_number;

        prevent_rest_check{1,inside_loop}(6,1)= fst_mv_oj_coor(5,1);           
        prevent_rest_check{1,inside_loop}(6,2)= fst_mv_oj_coor(5,2);
        prevent_rest_check{1,inside_loop}(6,3)= fst_mv_oj_coor(5,3);
    end
    stop_search=1; %prevent duplicate trajectories
else
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
        if(duplication_happen >= 2) 
            duplication_number=prevent_duplication_loop;
        end
    end
end
    
if(duplication_number ~= 0)
    if(angle_information(duplication_number,1) < angle_information(inside_loop,1))
        duplehappen=1;
    end
    if( duplehappen == 1)
        %duplication happen
    else %%%------------------------------check-------------------------------------------
        prevent_rest_check{1,inside_loop}(1,1) = fast_total_coor{1,inside_loop}(origin_size_fst,1);
        prevent_rest_check{1,inside_loop}(1,2) = fast_total_coor{1,inside_loop}(origin_size_fst,2);
        prevent_rest_check{1,inside_loop}(1,3) = fast_total_coor{1,inside_loop}(origin_size_fst,3); 
        prevent_rest_check{1,inside_loop}(2,1) = fst_mv_oj_coor(1,1);
        prevent_rest_check{1,inside_loop}(2,2) = fst_mv_oj_coor(1,2);
        prevent_rest_check{1,inside_loop}(2,3) = fst_mv_oj_coor(1,3);
        prevent_rest_check{1,inside_loop}(3,1) = fst_mv_oj_coor(2,1);
        prevent_rest_check{1,inside_loop}(3,2) = fst_mv_oj_coor(2,2); 
        prevent_rest_check{1,inside_loop}(3,3) = fst_mv_oj_coor(2,3); 
        prevent_rest_check{1,inside_loop}(4,1) = fst_mv_oj_coor(3,1);
        prevent_rest_check{1,inside_loop}(4,2) = fst_mv_oj_coor(3,2);
        prevent_rest_check{1,inside_loop}(4,3) = fst_mv_oj_coor(3,3);
        prevent_rest_check{1,inside_loop}(5,1) = fst_mv_oj_coor(4,1);
        prevent_rest_check{1,inside_loop}(5,2) = fst_mv_oj_coor(4,2);
        prevent_rest_check{1,inside_loop}(5,3) = fst_mv_oj_coor(4,3);

        fast_total_coor{1,inside_loop}(origin_size_fst+1,1) = fst_mv_oj_coor(1,1);
        fast_total_coor{1,inside_loop}(origin_size_fst+1,2) = fst_mv_oj_coor(1,2);       
        fast_total_coor{1,inside_loop}(origin_size_fst+1,3) = fst_mv_oj_coor(1,3); 
        fast_total_coor{1,inside_loop}(origin_size_fst+1,4) = jump_number;

        fast_total_coor{1,inside_loop}(origin_size_fst+2,1) = fst_mv_oj_coor(2,1);
        fast_total_coor{1,inside_loop}(origin_size_fst+2,2) = fst_mv_oj_coor(2,2); 
        fast_total_coor{1,inside_loop}(origin_size_fst+2,3) = fst_mv_oj_coor(2,3);
        fast_total_coor{1,inside_loop}(origin_size_fst+2,4) = jump_number;

        fast_total_coor{1,inside_loop}(origin_size_fst+3,1) = fst_mv_oj_coor(3,1);
        fast_total_coor{1,inside_loop}(origin_size_fst+3,2) = fst_mv_oj_coor(3,2); 
        fast_total_coor{1,inside_loop}(origin_size_fst+3,3) = fst_mv_oj_coor(3,3);
        fast_total_coor{1,inside_loop}(origin_size_fst+3,4) = jump_number;

        fast_total_coor{1,inside_loop}(origin_size_fst+4,1) = fst_mv_oj_coor(4,1);
        fast_total_coor{1,inside_loop}(origin_size_fst+4,2) = fst_mv_oj_coor(4,2);
        fast_total_coor{1,inside_loop}(origin_size_fst+4,3) = fst_mv_oj_coor(4,3);
        fast_total_coor{1,inside_loop}(origin_size_fst+4,4) = jump_number;

        if(~isnan(ans_xyz5))
            prevent_rest_check{1,inside_loop}(6,1)= fst_mv_oj_coor(5,1);           
            prevent_rest_check{1,inside_loop}(6,2)= fst_mv_oj_coor(5,2);
            prevent_rest_check{1,inside_loop}(6,3)= fst_mv_oj_coor(5,3);

            fast_total_coor{1,inside_loop}(origin_size_fst+5,1) = fst_mv_oj_coor(5,1);
            fast_total_coor{1,inside_loop}(origin_size_fst+5,2) = fst_mv_oj_coor(5,2);
            fast_total_coor{1,inside_loop}(origin_size_fst+5,3) = fst_mv_oj_coor(5,3);
            fast_total_coor{1,inside_loop}(origin_size_fst+5,4) = jump_number;
        end
        
        duple_mat=prevent_rest_check{1,duplication_number}(1,3);

        %delete duplication before                                                                                                                                    
        fast_total_coor{1,duplication_number}(origin_size_fst+1,1) = nan;
        fast_total_coor{1,duplication_number}(origin_size_fst+1,2) = nan;      
        fast_total_coor{1,duplication_number}(origin_size_fst+1,3) = nan;
        fast_total_coor{1,duplication_number}(origin_size_fst+1,4) = nan;

        fast_total_coor{1,duplication_number}(origin_size_fst+2,1) = nan;
        fast_total_coor{1,duplication_number}(origin_size_fst+2,2) = nan;        
        fast_total_coor{1,duplication_number}(origin_size_fst+2,3) = nan;
        fast_total_coor{1,duplication_number}(origin_size_fst+2,4) = nan;

        fast_total_coor{1,duplication_number}(origin_size_fst+3,1) = nan;
        fast_total_coor{1,duplication_number}(origin_size_fst+3,2) = nan;       
        fast_total_coor{1,duplication_number}(origin_size_fst+3,3) = nan;
        fast_total_coor{1,duplication_number}(origin_size_fst+3,4) = nan;

        fast_total_coor{1,duplication_number}(origin_size_fst+4,1) = nan;
        fast_total_coor{1,duplication_number}(origin_size_fst+4,2) = nan;
        fast_total_coor{1,duplication_number}(origin_size_fst+4,3) = nan; 
        fast_total_coor{1,duplication_number}(origin_size_fst+4,4) = nan;

        fast_total_coor{1,duplication_number}(origin_size_fst+5,1) = nan;
        fast_total_coor{1,duplication_number}(origin_size_fst+5,2) = nan;
        fast_total_coor{1,duplication_number}(origin_size_fst+5,3) = nan;
        fast_total_coor{1,duplication_number}(origin_size_fst+5,4) = nan;
        
        stop_search=1; %prevent duplicate trajectories
    end 
else
    prevent_rest_check{1,inside_loop}(1,1) = fast_total_coor{1,inside_loop}(origin_size_fst,1);
    prevent_rest_check{1,inside_loop}(1,2) = fast_total_coor{1,inside_loop}(origin_size_fst,2);
    prevent_rest_check{1,inside_loop}(1,3) = fast_total_coor{1,inside_loop}(origin_size_fst,3);
    prevent_rest_check{1,inside_loop}(2,1) = fst_mv_oj_coor(1,1);
    prevent_rest_check{1,inside_loop}(2,2) = fst_mv_oj_coor(1,2);
    prevent_rest_check{1,inside_loop}(2,3) = fst_mv_oj_coor(1,3); 
    prevent_rest_check{1,inside_loop}(3,1) = fst_mv_oj_coor(2,1);
    prevent_rest_check{1,inside_loop}(3,2) = fst_mv_oj_coor(2,2); 
    prevent_rest_check{1,inside_loop}(3,3) = fst_mv_oj_coor(2,3);
    prevent_rest_check{1,inside_loop}(4,1) = fst_mv_oj_coor(3,1);
    prevent_rest_check{1,inside_loop}(4,2) = fst_mv_oj_coor(3,2);
    prevent_rest_check{1,inside_loop}(4,3) = fst_mv_oj_coor(3,3);
    prevent_rest_check{1,inside_loop}(5,1) = fst_mv_oj_coor(4,1);
    prevent_rest_check{1,inside_loop}(5,2) = fst_mv_oj_coor(4,2);
    prevent_rest_check{1,inside_loop}(5,3) = fst_mv_oj_coor(4,3);

    fast_total_coor{1,inside_loop}(origin_size_fst+1,1) = fst_mv_oj_coor(1,1);
    fast_total_coor{1,inside_loop}(origin_size_fst+1,2) = fst_mv_oj_coor(1,2);       
    fast_total_coor{1,inside_loop}(origin_size_fst+1,3) = inside_loop;
    fast_total_coor{1,inside_loop}(origin_size_fst+1,4) = jump_number;

    fast_total_coor{1,inside_loop}(origin_size_fst+2,1) = fst_mv_oj_coor(2,1);
    fast_total_coor{1,inside_loop}(origin_size_fst+2,2) = fst_mv_oj_coor(2,2);         
    fast_total_coor{1,inside_loop}(origin_size_fst+2,3) = inside_loop;
    fast_total_coor{1,inside_loop}(origin_size_fst+2,4) = jump_number;

    fast_total_coor{1,inside_loop}(origin_size_fst+3,1) = fst_mv_oj_coor(3,1);
    fast_total_coor{1,inside_loop}(origin_size_fst+3,2) = fst_mv_oj_coor(3,2);       
    fast_total_coor{1,inside_loop}(origin_size_fst+3,3) = inside_loop;
    fast_total_coor{1,inside_loop}(origin_size_fst+3,4) = jump_number;

    fast_total_coor{1,inside_loop}(origin_size_fst+4,1) = fst_mv_oj_coor(4,1);
    fast_total_coor{1,inside_loop}(origin_size_fst+4,2) = fst_mv_oj_coor(4,2);
    fast_total_coor{1,inside_loop}(origin_size_fst+4,3) = inside_loop;  
    fast_total_coor{1,inside_loop}(origin_size_fst+4,4) = jump_number;

    if(~isnan(ans_xyz5))
        prevent_rest_check{1,inside_loop}(6,1)= fst_mv_oj_coor(5,1);           
        prevent_rest_check{1,inside_loop}(6,2)= fst_mv_oj_coor(5,2);
        prevent_rest_check{1,inside_loop}(6,3)= fst_mv_oj_coor(5,3);

        fast_total_coor{1,inside_loop}(origin_size_fst+5,1) = fst_mv_oj_coor(5,1);
        fast_total_coor{1,inside_loop}(origin_size_fst+5,2) = fst_mv_oj_coor(5,2);
        fast_total_coor{1,inside_loop}(origin_size_fst+5,3) = inside_loop;
        fast_total_coor{1,inside_loop}(origin_size_fst+5,4) = jump_number;
    end
    stop_search=1; %prevent duplicate trajectories
end  

if( stop_search == 0 )
    prevent_rest_check{1,inside_loop}(1,1) = nan;
    prevent_rest_check{1,inside_loop}(1,2) = nan;
    prevent_rest_check{1,inside_loop}(1,3) = nan;
end

if(duplication_number == 0)
    duple_mat=nan;
end

end