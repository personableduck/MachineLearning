wave=625;
lowest_depth=250/10;
highest_depth=450/10;

noise_crop=150;

z1=[];
count=0;
object_number=[];
for z=lowest_depth:highest_depth %start from 10um z2 gap
count=count+1;
z1(count)=z*10;
rec_r=Prop_SSA_v2(test_img,1.12,1.12,-z1(count),wave/1000); %reconstruction
    if focus_type == 0
        focus_img=abs(rec_r); %amplitude 
    elseif focus_type == 1  
        focus_img=angle(rec_r); %phase 
    else
%         display 'Set focus_type 0(amplitude) or 1(phase)'
%         break
    end    
    
[height,width]=size(focus_img);
focus_img=focus_img(noise_crop:height-noise_crop,noise_crop:width-noise_crop);

figure,imshow(focus_img,[])  
object_binary=  mean(focus_img(:)) + (sign(mean(focus_img(:)))*-1)*2*std(focus_img(:)) > focus_img;
figure,imshow(object_binary)

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

hold on
plot(stn_final{count}(:,1),stn_final{count}(:,2),'mo')
hold off

object_number(count)=length(stn_final{count});
end

figure(07075),plot(z1,object_number)
hold on
hline = refline([0 mean(object_number)]);
hline.Color = 'r';
xlabel('Step 1'); 
set(gcf,'Color',[1 1 1]);
hold off

find_start=0;
find_peak=0;
find_finish=0;
for h=1:count
    if find_start==0
        if mean(object_number) < object_number(h)
            lowest_depth=z1(h);
            find_start=1;
            display(lowest_depth) 
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
            display(highest_depth)
            find_finish=1;
            break
        end
    end
end

if find_finish==0
    display('Put higher highest_depth value')
end
    
lowest_depth=lowest_depth/10;
highest_depth=highest_depth/10;

z1=[];
count=0;
for z=lowest_depth:highest_depth %start from 10um z2 gap
count=count+1;
z1(count)=z*10;
% rec_r=Prop_SSA_v2(img_file,pixelsize,pixelsize,-z1(count),lambda_V/1000); %reconstruction
rec_r=Prop_SSA_v2(test_img,1.12,1.12,-z1(count),wave/1000); 
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

rec_r=Prop_SSA_v2(test_img,1.12,1.12,-z1(h),wave/1000); %reconstruction
if focus_type == 0
    focus_img=abs(rec_r); %amplitude 
elseif focus_type == 1  
    focus_img=angle(rec_r); %phase 
end    

figure,imshow(focus_img,[])