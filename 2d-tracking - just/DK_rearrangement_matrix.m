
for ctct= 2:iter-1

dknm= size(coor_from{1,1});
dkjj= dknm(1,1);

for dksel= 1: dkjj

obj_cor((iter-1)*2,2)=0;
    
obj_cor(1,1) = coor_from{1,1}(dksel,1); %for consider
obj_cor(1,2) = coor_from{1,1}(dksel,2);
obj_cor(2,1) = coor_to{1,1}(dksel,1);
obj_cor(2,2) = coor_to{1,1}(dksel,2);

ppp=0;

  for vv = 1:(iter-1)
  th_vx = max(abs(coor_to{1,vv} - coor_from{1,vv}));
  th_v1 = max(th_vx);
  ppp=ppp+1;
  cmp_th(ppp,1) = th_v1;
  end
  
threshold3 = max(cmp_th)+5; %threshold object tracking
threshold4=-(threshold3);
jjj=0;
ff_num= iter;
ff_n= iter-2;

pkpk=0;


for nkk= 2:(iter-1)

pkpk=pkpk+1;

    
fc2= size(coor_from{1,nkk});
jj2= fc2(1,1);

for jjk= 1:jj2
    
if((threshold4 < (obj_cor((pkpk*2),1) - coor_from{1,nkk}(jjk,1)) & (obj_cor((pkpk*2),1) - coor_from{1,nkk}(jjk,1)) < threshold3) & (threshold4 < (obj_cor((pkpk*2),2) - coor_from{1,nkk}(jjk,2)) & (obj_cor((pkpk*2),2) - coor_from{1,nkk}(jjk,2)) < threshold3))
    
    obj_cor((2*nkk)-1,1) = coor_from{1,nkk}(jjk,1);
    obj_cor((2*nkk)-1,2) = coor_from{1,nkk}(jjk,2);
    obj_cor((2*nkk),1) = coor_to{1,nkk}(jjk,1);
    obj_cor((2*nkk),2) = coor_to{1,nkk}(jjk,2);
    
 
end 

for detec= 1: (iter-1)*2
    
    if obj_cor(detec,1)==0 & obj_cor(detec,2) ==0
        obj_cor(detec,1)=obj_cor((detec-1),1);
        obj_cor(detec,2)=obj_cor((detec-1),2);
    end
end

end
end




obj_mat_cor(dksel,1)={obj_cor};
obj_cor=0;

end


%%%%%%%%number search
ktkt=0;
nush1=0;
nush2=0;
nush3=0;
nush4=0;


    
    numcot= size(coor_from{1,ctct});
    mnumc= numcot(1,1);
    
    dknm2= size(coor_from{1,1});
    dkjj2= dknm2(1,1);

for ctct2= 1: mnumc
for rvsc= 1:dkjj2
for rvsc2= 1:(iter-1)*2
        
if ( threshold4 <= coor_from{1,ctct}(ctct2,1) - obj_mat_cor{rvsc,1}(rvsc2,1) & coor_from{1,ctct}(ctct2,1) - obj_mat_cor{rvsc,1}(rvsc2,1) <= threshold3) & (threshold4 <= coor_from{1,ctct}(ctct2,2) - obj_mat_cor{rvsc,1}(rvsc2,2) & coor_from{1,ctct}(ctct2,2) - obj_mat_cor{rvsc,1}(rvsc2,2) <= threshold3)

ktkt=ktkt+1;
ktnl=1;
kbnl(ktkt,1)= ktnl;

else
    
    ktkt=ktkt+1;
    ktnl=0;
    kbnl(ktkt,1)= ktnl;
    
end

end
end

if kbnl == 0
    
nush1 = coor_from{1,ctct}(ctct2,1);
nush2 = coor_from{1,ctct}(ctct2,2);
nush3 = coor_to{1,ctct}(ctct2,1);
nush4 = coor_to{1,ctct}(ctct2,2);

end

ktkt=0;
kbnl=0;

if nush1 > 0 
    
numffc= size(coor_from{1,1});
mnuff= numffc(1,1);
            
mknumf= mnuff+1;
            
coor_from{1,1}(mknumf,1)= nush1;
coor_from{1,1}(mknumf,2)= nush2;
coor_to{1,1}(mknumf,1)= nush3;
coor_to{1,1}(mknumf,2)= nush4;


nush1 = 0;
nush2 = 0;
nush3 = 0;
nush4 = 0;

end
end
end


%%%%% onceagiagn 

dknm2= size(coor_from{1,1});
dkjj2= dknm2(1,1);

for dksel2= 1: dkjj2

obj_cor((iter-1)*2,2)=0;
    
obj_cor(1,1) = coor_from{1,1}(dksel2,1); %for consider
obj_cor(1,2) = coor_from{1,1}(dksel2,2);
obj_cor(2,1) = coor_to{1,1}(dksel2,1);
obj_cor(2,2) = coor_to{1,1}(dksel2,2);

ppp=0;

  for vv = 1:(iter-1)
  th_vx = max(abs(coor_to{1,vv} - coor_from{1,vv}));
  th_v1 = max(th_vx);
  ppp=ppp+1;
  cmp_th(ppp,1) = th_v1;
  end
  
threshold3 = max(cmp_th)+2; %threshold object tracking
threshold4=-(threshold3);
jjj=0;
ff_num= iter;
ff_n= iter-2;

pkpk=0;


for nkk= 2:(iter-1)

pkpk=pkpk+1;

    
fc22= size(coor_from{1,nkk});
jj22= fc22(1,1);

for jjk= 1:jj22
    
if((threshold4 < (obj_cor((pkpk*2),1) - coor_from{1,nkk}(jjk,1)) & (obj_cor((pkpk*2),1) - coor_from{1,nkk}(jjk,1)) < threshold3) & (threshold4 < (obj_cor((pkpk*2),2) - coor_from{1,nkk}(jjk,2)) & (obj_cor((pkpk*2),2) - coor_from{1,nkk}(jjk,2)) < threshold3))
    
    obj_cor((2*nkk)-1,1) = coor_from{1,nkk}(jjk,1);
    obj_cor((2*nkk)-1,2) = coor_from{1,nkk}(jjk,2);
    obj_cor((2*nkk),1) = coor_to{1,nkk}(jjk,1);
    obj_cor((2*nkk),2) = coor_to{1,nkk}(jjk,2);
    
 
end 

for detec= 1: (iter-1)*2
    
    if obj_cor(detec,1)==0 & obj_cor(detec,2) ==0
        obj_cor(detec,1)=obj_cor((detec-1),1);
        obj_cor(detec,2)=obj_cor((detec-1),2);
    end
end

end
end




obj_mat_cor(dksel2,1)={obj_cor};
obj_cor=0;

end

%%%%%%SPEED

disfps= 9; %fps
timess= 3.72;  %iter/fps;

dkcor= size(obj_mat_cor);
dkcormax= dkcor(1,1);


spsp=0;

for ababs = 1: dkcormax
spsp=spsp+1;
sumdis=0;

for objds = 1:((iter-1)*2-1)
    
    obdistances= sqrt((obj_mat_cor{ababs,1}(objds+1,1)-obj_mat_cor{ababs,1}(objds,1))^2 + (obj_mat_cor{ababs,1}(objds+1,2)-obj_mat_cor{ababs,1}(objds,2))^2); % um
    
    sumdis = obdistances + sumdis;
    
end

speedss = (sumdis*pixelsize)/timess; %um/s

obj_spd(spsp,1) = speedss;
end

figure
hist(obj_spd,25)
xlabel('Speed(um/sec)')
ylabel('Counts')
median(obj_spd)
mean(obj_spd)
dkcormax %objectnumber