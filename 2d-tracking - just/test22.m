cff=coor_from
ctt=coor_to


for randomd=1:iter-2

jjj=0;

tc= size(ctt{1,1});
    ii= tc(1,1);
    
    fc= size(cff{1,2});
    jj= fc(1,1);
    
    for r1=1:ii
        for r2=1:jj
            xx = ctt{1,randomd}(r1,1) - cff{1,randomd+1}(r2,1);
            yy = ctt{1,randomd}(r1,2) - cff{1,randomd+1}(r2,2); 
            
            if((threshold4 <= xx & xx <= threshold3) & (threshold4 <= yy & yy <= threshold3))
                cff{1,randomd+1}(r2,1)=0;
                cff{1,randomd+1}(r2,2)=0;     
            end
            end
    end
 ouou=0;
 for ahah2 =1:205
 
 if cff{1,2}(ahah2,1) > 0
 ouou=ouou+1;
 cff{1,1}(205+ouou,1)=cff{1,2}(ahah2,1);
 cff{1,1}(205+ouou,2)=cff{1,2}(ahah2,2);
 end
 end