
ttt1=onefr(1:820 , 1:1060);
ttt2=onefr(1233:2464 , 1:1640);
ttt3=onefr(1:1232 , 1641:3280);
ttt4=onefr(1233:2464 , 1641:3280);
%ttt=onefr;

close all
count=0;
standr=[]
meani=[]
noiser=[]
chvk=[]
z2_list=[]
dff=[]

tic
for z=360:390
count=count+1;
z2_list(count)=z;
rec_r= Prop_SSA(ttt1,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list(count),lambda_V/1000/n_glass);
%rec_r=Prop_SSA_v2(ttt1,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list(count),lambda_V/1000/n_glass);
figure,imshow(abs(rec_r),[])

ch=abs(rec_r);
ftv=  mean(ch(:))-2*std(ch(:)) > ch;
result=(ch).*ftv;
pixel_number(count)=sum(ftv(:));
GG = result(result~=0);

ta{count}=GG;

ch=abs(rec_r);
ftv=  mean(ch(:))-2*std(ch(:)) < ch;
result=(ch).*ftv;
GG = result(result~=0);

% figure,imshow(result,[])

tb{count}=GG;
end
avg_meana=[]
avg_meanb=[]
dff=[]

for ccct=1:count
qr=ta{ccct};
avg_meana(ccct)=mean(abs(qr(:)));
end

for ccct=1:count
qr=tb{ccct};
avg_meanb(ccct)=mean(abs(qr(:)));
end

dff=avg_meana-avg_meanb;
figure,plot(z2_list,dff)
hold on
hline = refline([0 max(dff)]);
hline.Color = 'r';
drawnow;

for h=1:count
if max(dff) == dff(h)
    z2_list(h)
end
end


for h=1:count
if min(abs(real(noiser3(:)))) == abs(real(noiser3(h)))
    z2_list(h)
    break
end
end

title(['Optimized Z2 is ', num2str(z2_list(h))]);
%figure(12234),scatter(z2_list,noiser)
hold on
hline = refline([0 min(abs(real(noiser3(:))))]);
hline.Color = 'r';
drawnow;


standr=[]
meani=[]
noiser=[]
chvk=[]
count=0;
GG=[]

for z=z2_list(h)-10:z2_list(h)+10
count=count+1;
z3_list(count)=z;
rec_r= Prop_SSA(ttt,pixelsize/interp_factor,pixelsize/interp_factor,-z3_list(count),lambda_V/1000/n_glass);

ch=((rec_r));
ftv=  mean(ch(:))-2*std(ch(:)) > ch;
result=abs(ch).*ftv;

figure,imshow(result,[])
title(['depth is ', num2str(z3_list(count))]);

ch=((rec_r));
ftv=  mean(ch(:))-2*std(ch(:)) < ch;
result=(ch).*ftv;
GG = result(result~=0);

noiser3(count)= (std(result(:)) / mean(result(:)));

% figure,imshow(result,[])
% title(['depth is ', num2str(z3_list(count))]);
% %xlabel(['Noise is:' ,num2str(noiser3(count))]); 
% xlabel(['Std: ',num2str(std(GG)),' Mean: ',num2str(mean(GG)), ' MAX: ',num2str(max(GG)) ,' Min: ',num2str(min(GG)) ,' Noise: ',num2str( std(GG) / mean(GG)) ]); 

end

for h=1:count
if max(noiser3(:)) == noiser3(h)
    h
    break
end
end

hold off
figure(12234),scatter(z3_list,noiser)
hold on
hline = refline([0 max(noiser(:))]);
hline.Color = 'r';
drawnow;

for h=1:count
if min(noiser(:)) == noiser(h)
    h
    break
end
end

standr=[]
meani=[]
noiser=[]
chvk=[]
count=0;

for z=1:20
    
count=count+1;
z4_list(count)=z3_list(h)-1+z*0.1;
rec_r= Prop_SSA(ttt,pixelsize/interp_factor,pixelsize/interp_factor,-z3_list(count),lambda_V/1000/n_glass);

standr(count)=std(rec_r(:));
meani(count)=mean(rec_r(:));
% noiser(count)= abs(real(standr(count) / meani(count)));
noiser(count)= (standr(count) / meani(count));
end

for h=1:count
if max(noiser(:)) == noiser(h)
    h
    break
end
end

hold off
figure(12234),scatter(z4_list,noiser)
hold on
hline = refline([0 max(noiser(:))]);
hline.Color = 'r';
drawnow;

standr=[]
meani=[]
noiser=[]
chvk=[]
count=0;

for z=1:20
    
count=count+1;
z5_list(count)=z4_list(h)-0.1+z*0.01;
rec_r= Prop_SSA(ttt,pixelsize/interp_factor,pixelsize/interp_factor,-z3_list(count),lambda_V/1000/n_glass);

standr(count)=std(rec_r(:));
meani(count)=mean(rec_r(:));
% noiser(count)= abs(real(standr(count) / meani(count)));
noiser(count)= (standr(count) / meani(count));
end

for h=1:count
if max(noiser(:)) == noiser(h)
    h
    break
end
end

hold off
figure(12234),scatter(z5_list,noiser)
hold on
hline = refline([0 max(noiser(:))]);
hline.Color = 'r';
drawnow;

z2_list=z5_list(h)

rec_r= Prop_SSA(ttt,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list,lambda_V/1000/n_glass);
figure,imshow(abs(real(rec_r)),[])  