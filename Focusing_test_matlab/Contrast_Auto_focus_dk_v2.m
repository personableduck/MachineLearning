% Made by DUCK-HA HWANG at March 15 2017 


function [z2_list] = Contrast_Auto_focus_dk(img_file,focus_type,lowest_depth,highest_depth,noise_crop,pixelsize,lambda_V)
% img_file // one channel image file 
% focus_type // 0:amplitude focus, 1:phase focus
% lowest_depth // lowest z2 depth :300~500
% highest_depth // highest z2 depth :300~500
% pixelsize // pixelsize um: 1.12um
% lambda_V // wavelength: 639 nm
% noise_crop // crop the boundary area because of reconstruction noise

%for example: z2_list = Contrast_Auto_focus_dk(img_file,0,300,500,150,1.12,639)

tic

%% 1 step
lowest_depth=lowest_depth/10;
highest_depth=highest_depth/10;

z1=[];
count=0;
object_number=[];
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

[height,width]=size(focus_img);
focus_img=focus_img(noise_crop:height-noise_crop,noise_crop:width-noise_crop);    
    
object_binary=  mean(focus_img(:)) + (sign(mean(focus_img(:)))*-1)*2*std(focus_img(:)) > focus_img; %object filter
filter_size = bwareafilt(object_binary,[8 36]);

stn_obj=regionprops(filter_size,'Centroid','MajorAxisLength','MinorAxisLength');
stn_minor=[stn_obj.MinorAxisLength] > 2;
stn_count=[stn_obj.MajorAxisLength] < 8;

stn_mid=stn_minor .* stn_count;
stn_mid=stn_mid.';

stn_coor=cell2mat({stn_obj.Centroid}');

stn_finding=[stn_coor stn_mid;];
[stn_row]= find(stn_finding(:,3) ~= 0);
stn_final{count}=stn_coor(stn_row,:);

object_number(count)=length(stn_final{count});
end


figure(07071),plot(z1,object_number)
hold on
hline = refline([0 mean(object_number)]);
hline.Color = 'r';
title('Searching Range');
xlabel('Step 1'); 
set(gcf,'Color',[1 1 1]);
drawnow;
hold off

find_start=0;
find_peak=0;
find_finish=0;
for h=1:count
    if find_start==0
        if mean(object_number) < object_number(h)
            lowest_depth=z1(h);
            find_start=1;
            display(['Object detection start: ', num2str(lowest_depth), 'um']) 
        end
    end
    if mean(object_number) > object_number(h)
        find_start=0;
    end
    if max(object_number) == object_number(h)
        find_peak=1;
    end
    if find_peak ==1
        if mean(object_number) > object_number(h)
            highest_depth=z1(h-1);
            display(['Object detection finish: ', num2str(highest_depth), 'um']) 
            find_finish=1;
            break
        end
    end
end

if find_finish==0
    display('Error: Put higher highest_depth value')
    z2_list=[];
    return
end

display(['Process 14.2857 %.......Depth range: ', num2str(lowest_depth), 'um and ',num2str(highest_depth), 'um'  ]) 
%% 2 step

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
    end    
    
[height,width]=size(focus_img);
focus_img=focus_img(noise_crop:height-noise_crop,noise_crop:width-noise_crop);    

object_binary=  mean(focus_img(:)) + (sign(mean(focus_img(:)))*-1)*2*std(focus_img(:)) > focus_img; %object filter
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
        display(['Process 28.5714 %.......Depth result: ', num2str(z1(h)),'um']) 
        break
    end
end

hold off
figure(07071),plot(z1,contrast_result)
hold on
hline = refline([0 max(contrast_result)]);
hline.Color = 'r';
title(['Depth is ', num2str(z1(h)),'um' ]);
xlabel('Step 2'); 
set(gcf,'Color',[1 1 1]);
drawnow;
hold off

%% 3 step
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

[height,width]=size(focus_img);
focus_img=focus_img(noise_crop:height-noise_crop,noise_crop:width-noise_crop);        
    
object_binary=  mean(focus_img(:)) + (sign(mean(focus_img(:)))*-1)*2*std(focus_img(:)) > focus_img; %object filter
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
        display(['Process 42.8571 %.......Depth result: ', num2str(z2(h)),'um']); 
        break
    end
end

hold off
figure(07071),plot(z2,contrast_result)
hold on
hline = refline([0 max(contrast_result)]);
hline.Color = 'r';
title(['Depth is ', num2str(z2(h)),'um' ]);
xlabel('Step 3'); 
set(gcf,'Color',[1 1 1]);
drawnow;
hold off

%% 4 step

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
    
[height,width]=size(focus_img);
focus_img=focus_img(noise_crop:height-noise_crop,noise_crop:width-noise_crop);    
 
object_binary=  mean(focus_img(:)) + (sign(mean(focus_img(:)))*-1)*2*std(focus_img(:)) > focus_img; %object filter
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
        display(['Process 57.1429 %.......Depth result: ', num2str(z3(h)),'um']); 
        break
    end
end

hold off
figure(07071),plot(z3,contrast_result)
hold on
hline = refline([0 max(contrast_result)]);
hline.Color = 'r';
title(['Depth is ', num2str(z3(h)),'um' ]);
xlabel('Step 4'); 
set(gcf,'Color',[1 1 1]);
drawnow;
hold off

%% 5 step

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
    
[height,width]=size(focus_img);
focus_img=focus_img(noise_crop:height-noise_crop,noise_crop:width-noise_crop);        
    
object_binary=  mean(focus_img(:)) + (sign(mean(focus_img(:)))*-1)*2*std(focus_img(:)) > focus_img; %object filter
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
        display(['Process 71.4286 %.......Depth result: ', num2str(z4(h)),'um']); 
        break
    end
end

hold off
figure(07071),plot(z4,contrast_result)
hold on
hline = refline([0 max(contrast_result)]);
hline.Color = 'r';
title(['Depth is ', num2str(z4(h)),'um' ]);
xlabel('Step 5'); 
set(gcf,'Color',[1 1 1]);
drawnow;
hold off

%% step 6
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
    
[height,width]=size(focus_img);
focus_img=focus_img(noise_crop:height-noise_crop,noise_crop:width-noise_crop);        
    
object_binary=  mean(focus_img(:)) + (sign(mean(focus_img(:)))*-1)*2*std(focus_img(:)) > focus_img; %object filter
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
        display(['Process 85.7143 %.......Depth result: ', num2str(z5(h)),'um']);
        break
    end
end

hold off
figure(07071),plot(z5,contrast_result)
hold on
hline = refline([0 max(contrast_result)]);
hline.Color = 'r';
title(['Depth is ', num2str(z5(h)),'um' ]);
xlabel('Step 6'); 
set(gcf,'Color',[1 1 1]);
drawnow;
hold off

%% step 7
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
    
[height,width]=size(focus_img);
focus_img=focus_img(noise_crop:height-noise_crop,noise_crop:width-noise_crop);        
      
object_binary=  mean(focus_img(:)) + (sign(mean(focus_img(:)))*-1)*2*std(focus_img(:)) > focus_img; %object filter
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
        display(['Process 100 %.......Depth result: ', num2str(z6(h)),'um']);
        break
    end
end

hold off
figure(07071),plot(z6,contrast_result)
hold on
hline = refline([0 max(contrast_result)]);
hline.Color = 'r';
title(['Depth is ', num2str(z6(h)),'um']);
xlabel('Step 7'); 
set(gcf,'Color',[1 1 1]);
drawnow;
hold off

toc

rec_r=Prop_SSA_v2(img_file,pixelsize,pixelsize,-z6(h),lambda_V/1000);
figure(27270),imshow(abs(rec_r),[])
set(gcf,'Color',[1 1 1]);
title(['Depth is ', num2str(z6(h)),'um']);
if focus_type == 0
    xlabel('Amplitude reconstructed image'); 
elseif focus_type == 1  
    xlabel('Phase reconstructed image');  
end  
drawnow;

z2_list=z6(h);
end