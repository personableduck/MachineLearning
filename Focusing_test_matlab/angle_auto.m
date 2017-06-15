bandwidth_bak = 20;    %16 <--before basler sensor// 90<--pi sensor// bandwidth for lowpassfilter background // you need to change
dark_bias = min(If{1,1}(:));
ori_level = mean2(If{1,1});
[onefr,myimg_filtered] = intensitycompensation(If{1,1},bandwidth_bak,bandwidth_bak, dark_bias);
onefr = onefr/mean2(onefr) * ori_level; %apply law pass filter to average hologram image

ttt1=onefr(1:1232 , 1:1640);
ttt2=onefr(1233:2464 , 1:1640);
ttt3=onefr(1:1232 , 1641:3280);
ttt4=onefr(1233:2464 , 1641:3280);
%ttt=onefr; avg_holo

set(gca, 'XLim', [1,1640], 'YLim', [1,1232])

%set(gca, 'XLim', [1600,2100], 'YLim', [350,750])
ttt3=ttt1;

close all
tic
z1_list=[]
count=0;
for z=30:50
count=count+1;
z1_list(count)=z*10;
%rec_r= Prop_SSA(ttt1,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list(count),lambda_V/1000/n_glass);
rec_r=Prop_SSA_v2(ttt3,pixelsize/interp_factor,pixelsize/interp_factor,-z1_list(count),lambda_V/1000/n_glass);

ch=(angle(rec_r));
ftv=  mean(ch(:))-2*std(ch(:)) > ch;
result=(ch).*ftv;
pixel_number(count)=sum(ftv(:));
GG = result(result~=0);

ta{count}=GG;

ch=(angle(rec_r));
ftv=  mean(ch(:))-2*std(ch(:)) < ch;
result=(ch).*ftv;
GG = result(result~=0);

tb{count}=GG;
end
avg_meana=[]
avg_meanb=[]
dff=[]

for ccct=1:count
qr=ta{ccct};
avg_meana(ccct)=mean(qr(:));
end

for ccct=1:count
qr=tb{ccct};
avg_meanb(ccct)=mean(qr(:));
end

dff=abs((avg_meana)-(avg_meanb));
figure,plot(z1_list,dff)
hold on
hline = refline([0 max(dff)]);
hline.Color = 'r';
drawnow;

for h=1:count
if max(dff) == dff(h)
    z1_list(h)
    break
end
end

z2_list=[]
count=0;
for z=1:20
count=count+1;
z2_list(count)=(z1_list(h)-10)+z;
%z2_list(count)=381.5823;
%rec_r= Prop_SSA(ttt1,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list(count),lambda_V/1000/n_glass);
rec_r=Prop_SSA_v2(ttt3,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list(count),lambda_V/1000/n_glass);
% figure,imshow(real(rec_r),[])
%figure,imshow(angle(rec_r),[])
%figure,imshow(abs(angle(rec_r)),[])

ch=(angle(rec_r));
ftv=  mean(ch(:))-2*std(ch(:)) > ch;
result=(ch).*ftv;
pixel_number(count)=sum(ftv(:));
GG = result(result~=0);

ta{count}=GG;

% figure,imshow(result,[])
% title(['depth is ', num2str(z2_list(count))]);
% %xlabel(['Noise is:' ,num2str(noiser3(count))]); 
% xlabel(['Std: ',num2str(std(GG)),' Mean: ',num2str(mean(GG)), ' MAX: ',num2str(max(GG)) ,' Min: ',num2str(min(GG)) ,' Noise: ',num2str( std(GG) / mean(GG)) ]); 

ch=(angle(rec_r));
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
strdara(ccct)=std(qr(:));
avg_meana(ccct)=mean(qr(:));
end

for ccct=1:count
qr=tb{ccct};
strdarb(ccct)=std(qr(:));
avg_meanb(ccct)=mean(qr(:));
end

dff=abs((avg_meana)-(avg_meanb));
figure,plot(z2_list,dff)
hold on
hline = refline([0 max(dff)]);
hline.Color = 'r';
drawnow;

for h=1:count
if max(dff) == dff(h)
    z2_list(h)
    break
end
end

count=0;

for z=1:20
count=count+1;
z3_list(count)=(z2_list(h)-1)+0.1*z;

%rec_r= Prop_SSA(ttt1,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list(count),lambda_V/1000/n_glass);
rec_r=Prop_SSA_v2(ttt3,pixelsize/interp_factor,pixelsize/interp_factor,-z3_list(count),lambda_V/1000/n_glass);
% figure,imshow(real(rec_r),[])
%figure,imshow(abs(rec_r),[])
% figure,imshow((angle(rec_r)),[])

ch=(angle(rec_r));
ftv=  mean(ch(:))-2*std(ch(:)) > ch;
result=(ch).*ftv;
pixel_number(count)=sum(ftv(:));
GG = result(result~=0);

ta{count}=GG;

% figure,imshow(result,[])
% title(['depth is ', num2str(z3_list(count))]);
% %xlabel(['Noise is:' ,num2str(noiser3(count))]); 
% xlabel(['Std: ',num2str(std(GG)),' Mean: ',num2str(mean(GG)), ' MAX: ',num2str(max(GG)) ,' Min: ',num2str(min(GG)) ,' Noise: ',num2str( std(GG) / mean(GG)) ]); 

ch=(angle(rec_r));
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
strdara(ccct)=std(qr(:));
avg_meana(ccct)=mean(qr(:));
% maxi_vala(ccct)=max(qr(:));
% minim_vala(ccct)= min(qr(:));
end

for ccct=1:count
qr=tb{ccct};
strdarb(ccct)=std(qr(:));
avg_meanb(ccct)=mean(qr(:));
% maxi_valb(ccct)=max(qr(:));
% minim_valb(ccct)= min(qr(:));
end

dff=abs((avg_meana)-(avg_meanb));
figure,plot(z3_list,dff)
hold on
hline = refline([0 max(dff)]);
hline.Color = 'r';
drawnow;

for h=1:count
if max(dff) == dff(h)
    z3_list(h)
    break
end
end

count=0;

for z=1:20
count=count+1;
z4_list(count)=(z3_list(h)-0.1)+0.01*z;

% z4_list(count)=370.5;
%rec_r= Prop_SSA(ttt1,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list(count),lambda_V/1000/n_glass);
rec_r=Prop_SSA_v2(ttt3,pixelsize/interp_factor,pixelsize/interp_factor,-z4_list(count),lambda_V/1000/n_glass);
% figure,imshow(real(rec_r),[])
%figure,imshow(abs(rec_r),[])
%figure,imshow((angle(rec_r)),[])

ch=(angle(rec_r));
ftv=  mean(ch(:))-2*std(ch(:)) > ch;
result=(ch).*ftv;
pixel_number(count)=sum(ftv(:));
GG = result(result~=0);

ta{count}=GG;

% figure,imshow(result,[])
% 
% 
% title(['depth is ', num2str(z4_list(count))]);
% %xlabel(['Noise is:' ,num2str(noiser3(count))]); 
% xlabel(['Std: ',num2str(std(GG)),' Mean: ',num2str(mean(GG)), ' MAX: ',num2str(max(GG)) ,' Min: ',num2str(min(GG)) ,' Noise: ',num2str( std(GG) / mean(GG)) ]); 

ch=(angle(rec_r));
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
strdara(ccct)=std(qr(:));
avg_meana(ccct)=mean(qr(:));
% maxi_vala(ccct)=max(qr(:));
% minim_vala(ccct)= min(qr(:));
end

for ccct=1:count
qr=tb{ccct};
strdarb(ccct)=std(qr(:));
avg_meanb(ccct)=mean(qr(:));
% maxi_valb(ccct)=max(qr(:));
% minim_valb(ccct)= min(qr(:));
end

dff=abs((avg_meana)-(avg_meanb));
figure,plot(z4_list,dff)
hold on
hline = refline([0 max(dff)]);
hline.Color = 'r';
drawnow;

for h=1:count
if max(dff) == dff(h)
    z4_list(h)
    break
end
end

count=0;

for z=1:20
count=count+1;
z5_list(count)=(z4_list(h)-0.01)+0.001*z;

% z4_list(count)=370.5;
%rec_r= Prop_SSA(ttt1,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list(count),lambda_V/1000/n_glass);
rec_r=Prop_SSA_v2(ttt3,pixelsize/interp_factor,pixelsize/interp_factor,-z5_list(count),lambda_V/1000/n_glass);
% figure,imshow(real(rec_r),[])
%figure,imshow(abs(rec_r),[])
%figure,imshow((angle(rec_r)),[])

ch=(angle(rec_r));
ftv=  mean(ch(:))-2*std(ch(:)) > ch;
result=(ch).*ftv;
pixel_number(count)=sum(ftv(:));
GG = result(result~=0);

ta{count}=GG;

% figure,imshow(result,[])


% title(['depth is ', num2str(z5_list(count))]);
% %xlabel(['Noise is:' ,num2str(noiser3(count))]); 
% xlabel(['Std: ',num2str(std(GG)),' Mean: ',num2str(mean(GG)), ' MAX: ',num2str(max(GG)) ,' Min: ',num2str(min(GG)) ,' Noise: ',num2str( std(GG) / mean(GG)) ]); 

ch=(angle(rec_r));
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
strdara(ccct)=std(qr(:));
avg_meana(ccct)=mean(qr(:));
% maxi_vala(ccct)=max(qr(:));
% minim_vala(ccct)= min(qr(:));
end

for ccct=1:count
qr=tb{ccct};
strdarb(ccct)=std(qr(:));
avg_meanb(ccct)=mean(qr(:));
% maxi_valb(ccct)=max(qr(:));
% minim_valb(ccct)= min(qr(:));
end

dff=abs((avg_meana)-(avg_meanb));
figure,plot(z5_list,dff)
hold on
hline = refline([0 max(dff)]);
hline.Color = 'r';
drawnow;

for h=1:count
if max(dff) == dff(h)
    z5_list(h)
    break
end
end

count=0;

for z=1:20
count=count+1;
z6_list(count)=(z5_list(h)-0.001)+0.0001*z;

% z4_list(count)=370.5;
%rec_r= Prop_SSA(ttt1,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list(count),lambda_V/1000/n_glass);
rec_r=Prop_SSA_v2(ttt3,pixelsize/interp_factor,pixelsize/interp_factor,-z6_list(count),lambda_V/1000/n_glass);
% figure,imshow(real(rec_r),[])
%figure,imshow(abs(rec_r),[])
%figure,imshow((angle(rec_r)),[])

ch=(angle(rec_r));
ftv=  mean(ch(:))-2*std(ch(:)) > ch;
result=(ch).*ftv;
pixel_number(count)=sum(ftv(:));
GG = result(result~=0);

ta{count}=GG;

% figure,imshow(result,[])


% title(['depth is ', num2str(z6_list(count))]);
% %xlabel(['Noise is:' ,num2str(noiser3(count))]); 
% xlabel(['Std: ',num2str(std(GG)),' Mean: ',num2str(mean(GG)), ' MAX: ',num2str(max(GG)) ,' Min: ',num2str(min(GG)) ,' Noise: ',num2str( std(GG) / mean(GG)) ]); 

ch=(angle(rec_r));
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
strdara(ccct)=std(qr(:));
avg_meana(ccct)=mean(qr(:));
% maxi_vala(ccct)=max(qr(:));
% minim_vala(ccct)= min(qr(:));
end

for ccct=1:count
qr=tb{ccct};
strdarb(ccct)=std(qr(:));
avg_meanb(ccct)=mean(qr(:));
% maxi_valb(ccct)=max(qr(:));
% minim_valb(ccct)= min(qr(:));
end

dff=abs((avg_meana)-(avg_meanb));
figure,plot(z6_list,dff)
hold on
hline = refline([0 max(dff)]);
hline.Color = 'r';
drawnow;

for h=1:count
if max(dff) == dff(h)
    z6_list(h)
    break
end
end
toc

rec_r=Prop_SSA_v2(ttt3,pixelsize/interp_factor,pixelsize/interp_factor,-z6_list(h),lambda_V/1000/n_glass);
figure,imshow(angle(rec_r),[])
% rec_r= Prop_SSA(ttt3,pixelsize/interp_factor,pixelsize/interp_factor,-z6_list(h),lambda_V/1000/n_glass);
% figure,imshow(abs(angle(rec_r)),[])