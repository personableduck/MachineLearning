%% Max get rid of grid pattern
close all

FLAG_select=0;
FLAG_pic_extrac=0;
FLAG_AF=1;
FLAG_recons=1;
FLAG_manual=0;

crop_number=3; %input for crop number (exmaple): 2 -> 4 field , 3-> 9 field , 4-> 16 field

if FLAG_select
selectFiles_tiff;
end

if FLAG_pic_extrac

read_dir = path_files;
length_frame_stack = length(file_1);

for iter=1:length_frame_stack
    
    filename=char(file_1(:,iter));
    [Ref_img] = read_img_wrapper_dk(read_dir,filename); %for color other sensor
%     raw_img_without_inter{1,iter}=Ref_img; %raw hologram image from dcraw
%     green_one1 = Ref_img((1:2:end),(2:2:end)); %first green
%     green_one_dundle{1,iter}=green_one1; %group of green channle
%     green_one2 = Ref_img((2:2:end),(1:2:end)); %second green
%     green_second_dundle{1,iter}=green_one2; %group of green channle
   
    [img_height,img_width]=size(Ref_img);
%     [green_height,green_width]=size(green_one1); % one green channel's size

    r = Ref_img(1:2:end,2:2:end) ./ Ref_img(2:2:end,1:2:end); 
    r_ = [r,flipud(r),r;fliplr(r),r,fliplr(r);r,flipud(r),r];
    h = fspecial('gaussian', [60 60], 30);
    rf = imfilter(r_, h);
    rf = rf(img_height/2+1:img_height, img_width/2+1:img_width); %this figures for half of image size// need to check
    Ref_img(2:2:end,1:2:end) =  Ref_img(2:2:end,1:2:end) .* rf; 

    If{1,iter} = InterpSonyImage(Ref_img); 
end

avg_holo = sperm_stack_avg_dk2(If,length_frame_stack);
% avg_raw_holo= sperm_stack_avg_dk2(raw_img_without_inter,length_frame_stack);
% avg_holo_one_green=sperm_stack_avg_dk2(green_one_dundle,length_frame_stack);
% avg_holo_second_green=sperm_stack_avg_dk2(green_second_dundle,length_frame_stack);

%for interporation average hologram
bandwidth_bak = 20;    %16 <--before basler sensor// 90<--pi sensor// bandwidth for lowpassfilter background // you need to change
dark_bias = min(avg_holo(:));
ori_level = mean2(avg_holo);
[avg_holo,myimg_filtered] = intensitycompensation(avg_holo,bandwidth_bak,bandwidth_bak, dark_bias);
avg_holo = avg_holo/mean2(avg_holo) * ori_level; %apply law pass filter to average hologram image

figure(1568),imagesc(myimg_filtered)
title('background image (low pass filtered)')
set(gcf,'Color',[1 1 1]);

% %for raw average hologram
% bandwidth_bak = 20;    %16 <--before basler sensor// 90<--pi sensor// bandwidth for lowpassfilter background // you need to change
% dark_bias = min(avg_raw_holo(:));
% ori_level = mean2(avg_raw_holo);
% [avg_raw_holo_filter,myimg_filtered] = intensitycompensation(avg_raw_holo,bandwidth_bak,bandwidth_bak, dark_bias);
% avg_raw_holo_filter = avg_raw_holo_filter/mean2(avg_raw_holo_filter) * ori_level; %apply law pass filter to average hologram image
% 
% % figure(1569),imagesc(myimg_filtered)
% % title('background image (low pass filtered)')
% % set(gcf,'Color',[1 1 1]);
% 
% %for one green channel average hologram
% bandwidth_bak = 20;    %16 <--before basler sensor// 90<--pi sensor// bandwidth for lowpassfilter background // you need to change
% dark_bias = min(avg_holo_one_green(:));
% ori_level = mean2(avg_holo_one_green);
% [avg_holo_one_green_filter,myimg_filtered] = intensitycompensation(avg_holo_one_green,bandwidth_bak,bandwidth_bak, dark_bias);
% avg_holo_one_green_filter = avg_holo_one_green_filter/mean2(avg_holo_one_green_filter) * ori_level; %apply law pass filter to average hologram image
% 
% % figure(1570),imagesc(myimg_filtered)
% % title('background image (low pass filtered)')
% % set(gcf,'Color',[1 1 1]);
% 
% %for one green channel average hologram
% bandwidth_bak = 20;    %16 <--before basler sensor// 90<--pi sensor// bandwidth for lowpassfilter background // you need to change
% dark_bias = min(avg_holo_second_green(:));
% ori_level = mean2(avg_holo_second_green);
% [avg_holo_second_green_filter,myimg_filtered] = intensitycompensation(avg_holo_second_green,bandwidth_bak,bandwidth_bak, dark_bias);
% avg_holo_second_green_filter = avg_holo_second_green_filter/mean2(avg_holo_second_green_filter) * ori_level; %apply law pass filter to average hologram image
% 
% % figure(1571),imagesc(myimg_filtered)
% % title('background image (low pass filtered)')
% % set(gcf,'Color',[1 1 1]);

end

if FLAG_AF

%% Auto_focus Max

% pixelsize=1.12;
% wavelength = 633;
% BP = @(x,z) Propagate(x, pixelsize, 1, wavelength, z, true, false, true);
% for ff = 1
%     I = avg_holo;
%     z = af_quick_v5(I,[300 600], pixelsize, 1,wavelength,0, 'GiniOfGrad', 1, -1, [], 10, 0.1, 40, true, true );
%     If = BP(I, z);
%     figure(1); clf; imshow(angle(If),[]); title(['Phase reconstruction, depth is ', num2str(z)]); drawnow;
%     disp(z);
% end


% exra_height=ceil(img_height*0.061); %150 pixel is extra region for delet wierd line afer reconstruction
% extra_width=ceil(img_width*0.046); %150 pixel is extra region for delet wierd line afer reconstruction

% disp('------------Using raw hologram-----------') 
% 
% crop_img1=avg_raw_holo_filter(1:img_height/2+exra_height,1:img_width/2+extra_width); %120 pixel is extra region for delet wierd line afer reconstruction
% z1 = Auto_Focus(crop_img1,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1); % m
% z1p=['First crop field depth is ',num2str(z1*1e6) 'um'];
% disp(z1p)
% 
% crop_img2=avg_raw_holo_filter(1:img_height/2+exra_height,img_width/2+1-extra_width:img_width);
% z2 = Auto_Focus(crop_img2,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1); % m
% z2p=['Second crop field depth is ',num2str(z2*1e6) 'um'];
% disp(z2p)
% 
% crop_img3=avg_raw_holo_filter(img_height/2+1-exra_height:img_height,1:img_width/2+extra_width);
% z3 = Auto_Focus(crop_img3,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1); % m
% z3p=['Third crop field depth is ',num2str(z3*1e6) 'um'];
% disp(z3p)
% 
% crop_img4=avg_raw_holo_filter(img_height/2+1-exra_height:img_height,img_width/2+1-extra_width:img_width);
% z4 = Auto_Focus(crop_img4,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1); % m
% z4p=['Fourth crop field depth is ',num2str(z4*1e6) 'um'];
% disp(z4p)

% exra_height=ceil(img_height*0.061); %150 pixel is extra region for delet wierd line afer reconstruction
% extra_width=ceil(img_width*0.046); %150 pixel is extra region for delet wierd line afer reconstruction

disp('------------Using interpolation hologram-----------') 
pixelsize=1.12;
wavelength = 633;

exra_height=200; %ceil( (img_height/crop_number)*0.13 ); %150 pixel is extra region for delet wierd line afer reconstruction
extra_width=200; %ceil( (img_width/crop_number)*0.092 ); %150 pixel is extra region for delet wierd line afer reconstruction

crop_img=[];
image_number=0;

BP = @(x,z) Propagate(x, pixelsize, 1, wavelength, z, true, false, true);
horizon_combine=[];
vertical_combine=[];

crop_img=[];
If=[];
final_crop=[];

for i=1:crop_number
    horizon_combine=[];
    
    for j=1:crop_number
        image_number=image_number+1;
        if i == 1
            row_value=[1:ceil((img_height/crop_number)+exra_height)];
            crop_row=[1:ceil(img_height/crop_number)];
        elseif i == crop_number
            row_value=[ceil((img_height/crop_number)*(i-1)-exra_height):img_height];
            crop_row=[exra_height+1:ceil(img_height/crop_number)+exra_height];
        else
            row_value=[ ceil((img_height/crop_number)*(i-1)-exra_height) : ceil(((img_height/crop_number)*i)+exra_height)];
            crop_row=[exra_height+1:ceil(img_height/crop_number)+exra_height];
        end
        
        if j == 1
            colum_value=[1:ceil((img_width/crop_number)+extra_width)];
            crop_colum=[1:ceil(img_width/crop_number)];
        elseif j == crop_number
            colum_value=[ceil((img_width/crop_number)*(j-1)-extra_width):img_width];
            crop_colum=[extra_width+1:ceil(img_width/crop_number)+extra_width];
        else
            colum_value=[ceil((img_width/crop_number)*(j-1)-extra_width) : ceil(((img_width/crop_number)*j)+extra_width)];
            crop_colum=[extra_width+1:ceil(img_width/crop_number)+extra_width];
        end
        
        crop_img{1,image_number}=avg_holo(row_value,colum_value);
        z = af_quick_v5(crop_img{1,image_number},[300 600], pixelsize, 1,wavelength,0, 'GiniOfGrad', 1, -1, [], 10, 0.1, 40, true, true );
        If{1,image_number} = BP(crop_img{1,image_number}, z);
        final_crop{1,image_number} = If{1,image_number}(crop_row,crop_colum);
        
        horizon_combine=horzcat(horizon_combine,final_crop{1,image_number});
     
        figure; clf; imshow(angle(If{1,image_number}),[]); title([num2str(image_number), 'Phase reconstruction, depth is ', num2str(z)]); drawnow;
        z_p=[num2str(image_number),' Crop field depth is ',num2str(z) 'um'];
        disp(z_p)
    end
    
    vertical_combine=[vertical_combine; horizon_combine;];
end
% 
% entire_field=horzcat(half1,half2);
% 
% figure(9452),imshow(angle(entire_field),[])
% title('Phase reconstructed image')
% set(gcf,'Color',[1 1 1]);

half1=[rec_r1(1:img_height/2,1:img_width/2); rec_r3(exra_height+1:end,1:img_width/2);];
half2=[rec_r2(1:img_height/2,extra_width+1:end); rec_r4(exra_height+1:end,extra_width+1:end);];

crop_img1=avg_holo(1:img_height/2+exra_height,1:img_width/2+extra_width); %120 pixel is extra region for delet wierd line afer reconstruction
z1 = Auto_Focus(crop_img1,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1); % m
title(['depth is ', num2str(z1*1e6)]); 
z1p=['First crop field depth is ',num2str(z1*1e6) 'um'];
disp(z1p)

crop_img2=avg_holo(1:img_height/2+exra_height,img_width/2+1-extra_width:img_width);
z2 = Auto_Focus(crop_img2,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1); % m
z2p=['Second crop field depth is ',num2str(z2*1e6) 'um'];
disp(z2p)

crop_img3=avg_holo(img_height/2+1-exra_height:img_height,1:img_width/2+extra_width);
z3 = Auto_Focus(crop_img3,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1); % m
z3p=['Third crop field depth is ',num2str(z3*1e6) 'um'];
disp(z3p)

crop_img4=avg_holo(img_height/2+1-exra_height:img_height,img_width/2+1-extra_width:img_width);
z4 = Auto_Focus(crop_img4,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1); % m
z4p=['Fourth crop field depth is ',num2str(z4*1e6) 'um'];
disp(z4p)



% exra_height=ceil(green_height*0.061); %150 pixel is extra region for delet wierd line afer reconstruction
% extra_width=ceil(green_width*0.046); %150 pixel is extra region for delet wierd line afer reconstruction

% disp('------------Using one green channel-----------') 
% 
% crop_img1=avg_holo_one_green_filter(1:green_height/2+exra_height,1:green_width/2+extra_width); %120 pixel is extra region for delet wierd line afer reconstruction
% z1 = Auto_Focus_negative(crop_img1,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1); % m
% title(['depth is ', num2str(z1*1e6)]); 
% z1p=['First crop field depth is ',num2str(z1*1e6) 'um'];
% disp(z1p)
% 
% crop_img2=avg_holo_one_green_filter(1:green_height/2+exra_height,green_width/2+1-extra_width:green_width);
% z2 = Auto_Focus_negative(crop_img2,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1); % m
% z2p=['Second crop field depth is ',num2str(z2*1e6) 'um'];
% disp(z2p)
% 
% crop_img3=avg_holo_one_green_filter(green_height/2+1-exra_height:green_height,1:green_width/2+extra_width);
% z3 = Auto_Focus_negative(crop_img3,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1); % m
% z3p=['Third crop field depth is ',num2str(z3*1e6) 'um'];
% disp(z3p)
% 
% crop_img4=avg_holo_one_green_filter(green_height/2+1-exra_height:green_height,green_width/2+1-extra_width:green_width);
% z4 = Auto_Focus_negative(crop_img4,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1); % m
% z4p=['Fourth crop field depth is ',num2str(z4*1e6) 'um'];
% disp(z4p)


% z = AF_chooseROI(avg_holo,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,1,1); % m

%z2 = Auto_Focus(avg_holo,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,1,1);

end

if FLAG_recons

%% reconstruction

pixelsize=pixel_size*1e6; % um
interp_factor=1; %shirinking images 
lambda_V=lambda*1e9; %nm
n_glass=1;

z2_list_1=z1*1e6; % um
z2_list_2=z2*1e6; % um
z2_list_3=z3*1e6; % um
z2_list_4=z4*1e6; % um
% z2_list=745; % um

rec_r1= Prop_SSA(crop_img1,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list_1,lambda_V/1000/n_glass);
rec_r2= Prop_SSA(crop_img2,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list_2,lambda_V/1000/n_glass);
rec_r3= Prop_SSA(crop_img3,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list_3,lambda_V/1000/n_glass);
rec_r4= Prop_SSA(crop_img4,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list_4,lambda_V/1000/n_glass);

half1=[rec_r1(1:img_height/2,1:img_width/2); rec_r3(exra_height+1:end,1:img_width/2);];
half2=[rec_r2(1:img_height/2,extra_width+1:end); rec_r4(exra_height+1:end,extra_width+1:end);];

entire_field=horzcat(half1,half2);

figure(9452),imshow(angle(entire_field),[])
title('Phase reconstructed image')
set(gcf,'Color',[1 1 1]);

figure(6527),imshow(abs(entire_field),[])
title('Amplitude reconstructed image')
set(gcf,'Color',[1 1 1]);
figure(6528),imshow(real(entire_field),[])
title('real field reconstructed image')
set(gcf,'Color',[1 1 1]);

end

if FLAG_manual

    for z=360:400
    
    z2_list=z;
    
    rec_r= Prop_SSA(crop_img2,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list,lambda_V/1000/n_glass);

    figure,imshow(abs(angle(rec_r)),[])       
    set(gca, 'XLim', [1400,1800], 'YLim', [150,550])
    title(['depth is ', num2str(z2_list)]);


    end

z2_list=373;
rec_r= Prop_SSA(crop_img3,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list,lambda_V/1000/n_glass);
figure,imshow(abs(rec_r),[]) 
% 
% a1=rec_r;

% z2_list=341.7616;
% rec_r= Prop_SSA(crop_img2,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list,lambda_V/1000/n_glass);
%     close all
%     figure,imshow((real(rec_r)),[]) 
%     % set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
%     title(['real part depth is ', num2str(z2_list)]);
%     figure,imshow(abs(rec_r),[]) 
%     % set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
%     title(['abs part depth is ', num2str(z2_list)]);
%     figure,imshow(angle(rec_r),[])       
%     % set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
%     title(['angle part depth is ', num2str(z2_list)]);
%     figure,imshow(imag(rec_r),[])  
%     % set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
%     title(['imag part depth is ', num2str(z2_list)]);
%     figure,imshow(real(rec_r),[])  
%     % set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
%     title(['real part depth is ', num2str(z2_list)]);
%     figure,imshow(abs(imag(rec_r)),[])  
%     % set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
%     title(['abs(img) part depth is ', num2str(z2_list)]);
%     figure,imshow(abs(real(rec_r)),[])  
%     % set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
%     title(['abs(real) part depth is ', num2str(z2_list)]);

end

% l=350e-6;
% u=400-6;
% delta_z =10;
% lambda = lambda_V*1e-6/1000/1;  %is 645 in um? or lambda_V/1000/1.5
% rect=[1200,1000,900,700];

% tic
% z2 = Auto_Focus_v5(avg_holo,l,u,delta_z,pixel_size,refidx,lambda,0,'Tamura',1,[],[],[],[]);
% toc

% z3 = Auto_Focus(avg_holo,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1);



