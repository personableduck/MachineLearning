FLAG_select=0;
FLAG_pic_extrac=0;
FLAG_manual=1;

if FLAG_select
selectFiles_tiff;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Image Extraction
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FLAG_pic_extrac

read_dir = path_files;
length_frame_stack = length(file_1);

for iter=1:length_frame_stack
    
    filename=char(file_1(:,iter));
    [Ref_img] = read_img_wrapper_dk2(read_dir,filename); %for color other sensor

    holo{1,iter} = Ref_img; 

end


avg_holo = sperm_stack_avg_dk2(holo,length_frame_stack);
figure,imshow(avg_holo,[])

bandwidth_bak = 18;    %16 <--before basler sensor// 90<--pi sensor// bandwidth for lowpassfilter background // you need to change
dark_bias = min(avg_holo(:));
ori_level = mean2(avg_holo);
[avg_holo,myimg_filtered] = intensitycompensation(avg_holo,bandwidth_bak,bandwidth_bak, dark_bias);
avg_holo = avg_holo/mean2(avg_holo) * ori_level; %apply law pass filter to average hologram image
figure,imshow(avg_holo,[])

end

if FLAG_manual
    
lambda = 625e-9; % m // Justin's LED: 633nm //New fibered LED (Chris'): 624nm
pixel_size = 1.12e-6; % m
refidx = 1;

pixelsize=pixel_size*1e6; % um
interp_factor=1; %shirinking images 
lambda_V=lambda*1e9; %nm
n_glass=1;

BP = @(x,z) Propagate(x, 1.12, 1, 625, z, true, false, true);

close all

for z=320:330

z2_list=z; %373.39

rec_r= Prop_SSA(avg_holo,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list,lambda_V/1000/n_glass);

figure,imshow(abs(real(rec_r)),[])       
set(gca, 'XLim', [1000,1500], 'YLim', [1150,1500])
% set(gca, 'XLim', [2150,2650], 'YLim', [1100,1500])
title(['depth is ', num2str(z2_list)]);


% figure,imshow(abs(real(rec_r)),[])       
% set(gca, 'XLim', [1600,2100], 'YLim', [350,750])
% title(['Second field depth is ', num2str(z2_list)]);

% foo=imresize(avg_holo(350:750,1600:2100),3);
% Prop_SSA(foo,pixelsize/interp_factor/3,pixelsize/interp_factor/3,-z2_list,lambda_V/1000/n_glass);
% figure,imshow(abs(real(ans)),[]);
% title(['depth is ', num2str(z2_list)]);

% figure,imagesc(abs(real(rec_r))); 
% axis equal;
% set(gca, 'XLim', [1000,1500], 'YLim', [1150,1500])
% title(['depth is ', num2str(z2_list)]);

% rec_2 = BP(avg_holo, z2_list);
% 
% figure,imshow((real(rec_2)),[])       
% set(gca, 'XLim', [1000,1500], 'YLim', [1150,1500])
% title(['depth is ', num2str(z2_list)]);

end

% z2_list=381;
% rec_r= Prop_SSA(avg_holo,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list,lambda_V/1000/n_glass);
% % close all
% figure,imshow(abs(real(rec_r)),[]) 
% % set(gca, 'XLim', [1600,2100], 'YLim', [350,750])
% set(gca, 'XLim', [1000,1500], 'YLim', [1150,1500])
% 

end