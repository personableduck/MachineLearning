ppp=0;

  for vv = 1:(iter-1)
  th_vx = max(abs(coor_to{1,vv} - coor_from{1,vv}));
  th_v1 = max(th_vx);
  ppp=ppp+1;
  cmp_th(ppp,1) = th_v1;
  end
  
threshold3 = max(cmp_th)+5;
threshold4=-(threshold3);
jjj=0;
ff_num= iter;
ff_n= iter-2;


for n2 = 1:ff_n %fame number
    
    nn=n2+1;
    
    tc= size(coor_to{1,n2});
    ii= tc(1,1);
    
    fc= size(coor_from{1,nn});
    jj= fc(1,1);
    
    for r1=1:ii
        for r2=1:jj
            xx = coor_to{1,n2}(r1,1) - coor_from{1,nn}(r2,1);
            yy = coor_to{1,n2}(r1,2) - coor_from{1,nn}(r2,2); 
            
            if((threshold4 <= xx & xx <= threshold3) & (threshold4 <= yy & yy <= threshold3))
                xs=coor_from{1,n2}(r1,1);
                ys=coor_from{1,n2}(r1,2);
                jjj = jjj+1;
                
                object_num(jjj,1)=xs; % x coordinate of moving object
                object_num(jjj,2)=ys; % y coordinate of moving object
                object_num(jjj,3)=n2; % matrix cool's number
                object_num(jjj,4)=jjj; % number of object
                object_num(jjj,5)=r1; % number of order for coor_to
                object_num(jjj,6)=r2; % number of order for coor_from
                
                object_count=object_num;

            end
            
            jkj = jjj;
            jkj2= jjj-1;

            if (jkj>1) & (object_num(jkj,1) == object_num(jkj2,1)) & (object_num(jkj,2) == object_num(jkj2,2))
                jjj= jjj-1;
            end
        end
    end
    
    object_mat(n2,1)={object_num};
    jjj=0;
    object_num=0;
    
end 


