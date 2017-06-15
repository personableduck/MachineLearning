selectFiles_jpg;

read_dir = path_files;
length_frame_stack = length(file_1);

for iter=1:length_frame_stack
    
    filename=char(file_1(:,iter));
    [Ref_img] = read_img_wrapper_dk2(read_dir,filename); 
    raw_stack{1,iter} = Ref_img;
      
%     %look later
%     bandwidth_bak = 18;    %16 <--before basler sensor// 90<--pi sensor// bandwidth for lowpassfilter background // you need to change
%     dark_bias = min(raw_stack{1,iter}(:));
%     ori_level = mean2(raw_stack{1,iter});
%     [r01{1,iter},myimg_filtered] = intensitycompensation(raw_stack{1,iter},bandwidth_bak,bandwidth_bak, dark_bias);
%     result{1,iter} = r01{1,iter}/mean2(r01{1,iter}) * ori_level; %apply law pass filter to average hologram image
%     
% %     ROI_mask = roipoly(r0);
%     r0_crop = r0(930:1120,600:1000);
%     noise_rat = std(r0_crop(:))/mean(r0_crop(:))
    
 
%     figure,
%     imshow(If{1,iter},[])
%     title(['Green interporation ', num2str(iter)])
%     axis equal;
%     colorbar;
%     caxis([0 5])
    

%     hFig = figure(100+iter);
%     set(hFig, 'Position', [1 1 1700 700])
%     imagesc(r0{1,iter});
%     axis equal;
%     caxis([0 1023]);
%     colorbar;
%     set(gcf,'Color',[1 1 1]);
%     title(filename,'interpreter','none');

%     figure,
%     imshow(r0{1,iter},[]);
%     title([filename,'  Only green']);
    
%     figure,
%     imshow(avg_holo1{1,iter},[])
%     title(['interporation filtered image',num2str(iter)])
%     axis equal;
%     colorbar;
%     caxis([0 1000])
end

% avg_holo= sperm_stack_avg_dk2(result,length_frame_stack);
% figure,
% imshow(avg_holo,[]);
% title([filename,'  Only green']);
% 
% pixel_size = 1.12e-6; % m    
% lambda = 633e-9; 
% pixelsize=pixel_size*1e6; % um
% interp_factor=1; %shirinking images 
% lambda_V=lambda*1e9; %nm
% n_glass=1;    
% 
% for z=360:380
%     
%     z2_list=z;
%     
%     rec_r= Prop_SSA(avg_holo,pixelsize,pixelsize,-z2_list,lambda_V/1000/n_glass);
%     figure,imshow(abs(real(rec_r)),[])       
%     %set(gca, 'XLim', [600,1000], 'YLim', [450,750])
%     title(['depth is ', num2str(z2_list)]);
% 
% end