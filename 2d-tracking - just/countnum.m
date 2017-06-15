cc_f2 = coor_from;
cc_t2 = coor_to;
ppp=0;
ccfs= size(coor_from{1,1});
enclfs= ccfs(1,1);
adadd=0;

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
    
    tc= size(cc_t2{1,n2});
    ii= tc(1,1);
    
    fc= size(cc_f2{1,nn});
    jj= fc(1,1);
 
    for r1=1:ii
        for r2=1:jj
            xx1 = cc_t2{1,n2}(r1,1) - cc_f2{1,nn}(r2,1);
            yy1 = cc_t2{1,n2}(r1,2) - cc_f2{1,nn}(r2,2); 
            
            kcom(jj,1) = 0;
            
            if((threshold4 <= xx1 & xx1 <= threshold3) & (threshold4 <= yy1 & yy1 <= threshold3))

            kcom(r2,1) = 1;
            
            end
            
        end
    end
    
    
    for compare_kcom=1:jj
        
        if kcom(compare_kcom,1) == 0
           
            adadd=adadd+1;
            cc_f2{1,1}(enclfs+adadd,1) = coor_from{1,nn}(compare_kcom,1);
            cc_f2{1,1}(enclfs+adadd,2) = coor_from{1,nn}(compare_kcom,2);
            cc_t2{1,1}(enclfs+adadd,1) = coor_from{1,nn}(compare_kcom,1);
            cc_t2{1,1}(enclfs+adadd,2) = coor_from{1,nn}(compare_kcom,2);
            
            ccts= size(cc_t2{1,n2+1});
            enclts= ccts(1,1);
            cc_t2{1,1}(enclts+adadd,1) = coor_from{1,nn}(compare_kcom,1);
            cc_t2{1,1}(enclts+adadd,2) = coor_from{1,nn}(compare_kcom,2);
        
        end
    end
    
    kcom=0;
end
