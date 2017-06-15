
selectFiles_tiff;

read_dir = path_files;
length_frame_stack = length(file_1);

for iter=1:length_frame_stack
    
    filename=char(file_1(:,iter));
    [Ref_img] = read_img_wrapper_dk(read_dir,filename); %for color other sensor
    [img_width,img_height]=size(Ref_img);
    
    Ic{1,iter} = InterpSonyImage(Ref_img); 
    
    r = Ref_img(1:2:end,2:2:end) ./ Ref_img(2:2:end,1:2:end); 
    r_ = [r,flipud(r),r;fliplr(r),r,fliplr(r);r,flipud(r),r];
    h = fspecial('gaussian', [60 60], 30);
    rf = imfilter(r_, h);

    rf = rf(img_width/2+1:img_width, img_height/2+1:img_height); %this figures for half of image size// need to check
    
    Ref_img(2:2:end,1:2:end) =  Ref_img(2:2:end,1:2:end) .* rf; 

    If{1,iter} = InterpSonyImage(Ref_img); 
 
%     filename=char(file_1(:,iter));
%     [Ref_img] = read_img_wrapper_dk2(read_dir,filename); %for color other sensor
%     If{1,iter}=Ref_img;
end

If{1,1} = sperm_stack_avg_dk2(If,length_frame_stack);

bandwidth_bak = 20;    %16 <--before basler sensor// 90<--pi sensor// bandwidth for lowpassfilter background // you need to change
dark_bias = min(If{1,1}(:));
ori_level = mean2(If{1,1});
[onefr,myimg_filtered] = intensitycompensation(If{1,1},bandwidth_bak,bandwidth_bak, dark_bias);
onefr = onefr/mean2(onefr) * ori_level; %apply law pass filter to average hologram image

test_img=onefr(1:1000, 1:1500);
z2_list = Contrast_Auto_focus_dk_v2(test_img,0,250,450,150,1.12,625); %laser 639, chris:624 max:625

% close all
% 
% for z=30:50
%        
%     z2_list=z*10;
%     rec_r= Prop_SSA(test_img,1.12,1.12,-z2_list,639/1000);
%   
%     figure,imshow(abs(rec_r),[]);       
%     title(['depth is ', num2str(z2_list)]);
% end
