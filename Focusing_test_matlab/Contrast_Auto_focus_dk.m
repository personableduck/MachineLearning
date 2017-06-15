% Made by DUCK-HA HWANG at March 15 2017 


function [z2_list] = Contrast_Auto_focus_dk(img_file,focus_type,lowest_depth,highest_depth,pixelsize,lambda_V)
% img_file // one channel image file 
% focus_type // 0:amplitude focus, 1:phase focus
% lowest_depth // lowest z2 depth :300~500
% highest_depth // highest z2 depth :300~500
% pixelsize // pixelsize um: 1.12um
% lambda_V // wavelength: 639 nm

%for example: z2_list = Contrast_Auto_focus_dk(img_file,0,300,500,1.12,639)

tic

%% 1 step
lowest_depth=lowest_depth/10;
highest_depth=highest_depth/10;

z1=[];
count=0;
for z=lowest_depth:highest_depth %start from 10um z2 gap
count=count+1;
z1(count)=z*10;
rec_r=Prop_SSA_v2(img_file,pixelsize,pixelsize,-z1(count),lambda_V/1000); %reconstruction
    if focus_type == 0
        focus_img=abs(rec_r); %amplitude 
    elseif focus_type == 1  
        focus_img=angle(rec_r); %phase 
    else
        display 'Set focus_type 0(amplitude) or 1(phase)'
        break
    end    
object_binary=  mean(focus_img(:))-2*std(focus_img(:)) > focus_img; %object filter
intensity_obj=(focus_img).*object_binary; %restore object intensity
object_result{count} = intensity_obj(intensity_obj~=0); %extract intensity value except 0

background_binary=  mean(focus_img(:))-2*std(focus_img(:)) < focus_img; %background filter
intensity_bkg=(focus_img).*background_binary; %restore background intensity
background_result{count} = intensity_bkg(intensity_bkg~=0); %extract intensity value except 0
end

average_obj=[]; 
average_bck=[];
contrast_result=[];

for i=1:count
    average_calculation=object_result{i};
    average_obj(i)=mean(average_calculation(:));
end

for j=1:count
    average_calculation=background_result{j};
    average_bck(j)=mean(average_calculation(:));
end

contrast_result=abs((average_obj)-(average_bck)); % contrast result from difference between objects' average intensity and background average intensity

for h=1:count
    if max(contrast_result) == contrast_result(h)
        display(['Process 16.6667 %.......Depth result: ', num2str(z1(h))]) 
        break
    end
end

hold off
figure(07071),plot(z1,contrast_result)
hold on
hline = refline([0 max(contrast_result)]);
hline.Color = 'r';
drawnow;
title(['Depth is ', num2str(z1(h)) ]);
xlabel('Step 1'); 
set(gcf,'Color',[1 1 1]);
hold off

%% 2 step
z2=[];
count=0;

for z=1:20 %test 1um of z2 gap
count=count+1; %
z2(count)=(z1(h)-10)+z; 
rec_r=Prop_SSA_v2(img_file,pixelsize,pixelsize,-z2(count),lambda_V/1000);
    if focus_type == 0
        focus_img=abs(rec_r); %amplitude 
    elseif focus_type == 1  
        focus_img=angle(rec_r); %phase 
    end   
object_binary=  mean(focus_img(:))-2*std(focus_img(:)) > focus_img; %object filter
intensity_obj=(focus_img).*object_binary; %restore object intensity
object_result{count} = intensity_obj(intensity_obj~=0); %extract intensity value except 0

background_binary=  mean(focus_img(:))-2*std(focus_img(:)) < focus_img; %background filter
intensity_bkg=(focus_img).*background_binary; %restore background intensity
background_result{count} = intensity_bkg(intensity_bkg~=0); %extract intensity value except 0
end

average_obj=[]; 
average_bck=[];
contrast_result=[];

for i=1:count
    average_calculation=object_result{i};
    average_obj(i)=mean(average_calculation(:));
end

for j=1:count
    average_calculation=background_result{j};
    average_bck(j)=mean(average_calculation(:));
end

contrast_result=abs((average_obj)-(average_bck)); % contrast result from difference between objects' average intensity and background average intensity

for h=1:count
    if max(contrast_result) == contrast_result(h)
        display(['Process 33.3333 %.......Depth result: ', num2str(z2(h))]); 
        break
    end
end

figure(07071),plot(z2,contrast_result)
hold on
hline = refline([0 max(contrast_result)]);
hline.Color = 'r';
drawnow;
title(['Depth is ', num2str(z2(h)) ]);
xlabel('Step 2'); 
set(gcf,'Color',[1 1 1]);
hold off

%% 3 step

z3=[];
count=0;

for z=1:20 %test 1um of z2 gap
count=count+1; %
z3(count)=(z2(h)-1)+z*0.1; 
rec_r=Prop_SSA_v2(img_file,pixelsize,pixelsize,-z3(count),lambda_V/1000);
    if focus_type == 0
        focus_img=abs(rec_r); %amplitude 
    elseif focus_type == 1  
        focus_img=angle(rec_r); %phase 
    end    
object_binary=  mean(focus_img(:))-2*std(focus_img(:)) > focus_img; %object filter
intensity_obj=(focus_img).*object_binary; %restore object intensity
object_result{count} = intensity_obj(intensity_obj~=0); %extract intensity value except 0

background_binary=  mean(focus_img(:))-2*std(focus_img(:)) < focus_img; %background filter
intensity_bkg=(focus_img).*background_binary; %restore background intensity
background_result{count} = intensity_bkg(intensity_bkg~=0); %extract intensity value except 0
end

average_obj=[]; 
average_bck=[];
contrast_result=[];

for i=1:count
    average_calculation=object_result{i};
    average_obj(i)=mean(average_calculation(:));
end

for j=1:count
    average_calculation=background_result{j};
    average_bck(j)=mean(average_calculation(:));
end

contrast_result=abs((average_obj)-(average_bck)); % contrast result from difference between objects' average intensity and background average intensity

for h=1:count
    if max(contrast_result) == contrast_result(h)
        display(['Process 50 %.......Depth result: ', num2str(z3(h))]); 
        break
    end
end

figure(07071),plot(z3,contrast_result)
hold on
hline = refline([0 max(contrast_result)]);
hline.Color = 'r';
drawnow;
title(['Depth is ', num2str(z3(h)) ]);
xlabel('Step 3'); 
set(gcf,'Color',[1 1 1]);
hold off

%% 4 step

z4=[];
count=0;

for z=1:20
count=count+1;
z4(count)=(z3(h)-0.1)+0.01*z;
rec_r=Prop_SSA_v2(img_file,pixelsize,pixelsize,-z4(count),lambda_V/1000);
    if focus_type == 0
        focus_img=abs(rec_r); %amplitude 
    elseif focus_type == 1  
        focus_img=angle(rec_r); %phase 
    end    
object_binary=  mean(focus_img(:))-2*std(focus_img(:)) > focus_img; %object filter
intensity_obj=(focus_img).*object_binary; %restore object intensity
object_result{count} = intensity_obj(intensity_obj~=0); %extract intensity value except 0

background_binary=  mean(focus_img(:))-2*std(focus_img(:)) < focus_img; %background filter
intensity_bkg=(focus_img).*background_binary; %restore background intensity
background_result{count} = intensity_bkg(intensity_bkg~=0); %extract intensity value except 0
end

average_obj=[]; 
average_bck=[];
contrast_result=[];

for i=1:count
    average_calculation=object_result{i};
    average_obj(i)=mean(average_calculation(:));
end

for j=1:count
    average_calculation=background_result{j};
    average_bck(j)=mean(average_calculation(:));
end

contrast_result=abs((average_obj)-(average_bck)); % contrast result from difference between objects' average intensity and background average intensity

for h=1:count
    if max(contrast_result) == contrast_result(h)
        display(['Process 66.6667 %.......Depth result: ', num2str(z4(h))]); 
        break
    end
end

figure(07071),plot(z4,contrast_result)
hold on
hline = refline([0 max(contrast_result)]);
hline.Color = 'r';
drawnow;
title(['Depth is ', num2str(z4(h)) ]);
xlabel('Step 4'); 
set(gcf,'Color',[1 1 1]);
hold off

%% step 5
z5=[];
count=0;

for z=1:20
count=count+1;
z5(count)=(z4(h)-0.01)+0.001*z;
rec_r=Prop_SSA_v2(img_file,pixelsize,pixelsize,-z5(count),lambda_V/1000);
    if focus_type == 0
        focus_img=abs(rec_r); %amplitude 
    elseif focus_type == 1  
        focus_img=angle(rec_r); %phase 
    end    
object_binary=  mean(focus_img(:))-2*std(focus_img(:)) > focus_img; %object filter
intensity_obj=(focus_img).*object_binary; %restore object intensity
object_result{count} = intensity_obj(intensity_obj~=0); %extract intensity value except 0

background_binary=  mean(focus_img(:))-2*std(focus_img(:)) < focus_img; %background filter
intensity_bkg=(focus_img).*background_binary; %restore background intensity
background_result{count} = intensity_bkg(intensity_bkg~=0); %extract intensity value except 0
end

average_obj=[]; 
average_bck=[];
contrast_result=[];

for i=1:count
    average_calculation=object_result{i};
    average_obj(i)=mean(average_calculation(:));
end

for j=1:count
    average_calculation=background_result{j};
    average_bck(j)=mean(average_calculation(:));
end

contrast_result=abs((average_obj)-(average_bck)); % contrast result from difference between objects' average intensity and background average intensity

for h=1:count
    if max(contrast_result) == contrast_result(h)
        display(['Process 83.3333 %.......Depth result: ', num2str(z5(h))]);
        break
    end
end

figure(07071),plot(z5,contrast_result)
hold on
hline = refline([0 max(contrast_result)]);
hline.Color = 'r';
drawnow;
title(['Depth is ', num2str(z5(h)) ]);
xlabel('Step 5'); 
set(gcf,'Color',[1 1 1]);
hold off

%% step 6
z6=[];
count=0;

for z=1:20
count=count+1;
z6(count)=(z5(h)-0.001)+0.0001*z;
rec_r=Prop_SSA_v2(img_file,pixelsize,pixelsize,-z6(count),lambda_V/1000);
    if focus_type == 0
        focus_img=abs(rec_r); %amplitude 
    elseif focus_type == 1  
        focus_img=angle(rec_r); %phase 
    end    
object_binary=  mean(focus_img(:))-2*std(focus_img(:)) > focus_img; %object filter
intensity_obj=(focus_img).*object_binary; %restore object intensity
object_result{count} = intensity_obj(intensity_obj~=0); %extract intensity value except 0

background_binary=  mean(focus_img(:))-2*std(focus_img(:)) < focus_img; %background filter
intensity_bkg=(focus_img).*background_binary; %restore background intensity
background_result{count} = intensity_bkg(intensity_bkg~=0); %extract intensity value except 0
end

average_obj=[]; 
average_bck=[];
contrast_result=[];

for i=1:count
    average_calculation=object_result{i};
    average_obj(i)=mean(average_calculation(:));
end

for j=1:count
    average_calculation=background_result{j};
    average_bck(j)=mean(average_calculation(:));
end

contrast_result=abs((average_obj)-(average_bck)); % contrast result from difference between objects' average intensity and background average intensity

for h=1:count
    if max(contrast_result) == contrast_result(h)
        display(['Process 100 %.......Depth result: ', num2str(z6(h))]);
        break
    end
end

figure(07071),plot(z6,contrast_result)
hold on
hline = refline([0 max(contrast_result)]);
hline.Color = 'r';
drawnow;
title(['Depth is ', num2str(z6(h)) ]);
xlabel('Step 6'); 
set(gcf,'Color',[1 1 1]);
hold off

toc

rec_r=Prop_SSA_v2(img_file,pixelsize,pixelsize,-z6(h),lambda_V/1000);
figure(27270),imshow(abs(rec_r),[])
set(gcf,'Color',[1 1 1]);
drawnow;

z2_list=z6(h);
end