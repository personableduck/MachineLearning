


%% max get rid of grid pattern

FLAG_select=1;
FLAG_pic_extrac=1;
FLAG_AF=0;
FLAG_recons=0;
FLAG_manual=0;

FLAG_plane=1;
K_fac=0.5;
FLAG_mean_plane=0;

divide_ROI = 4; %Number of lines/columns that divide the FOV: corresponds to the following number of ROIs: divide_ROI^2


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Image Selection
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    [Ref_img] = read_img_wrapper_dk(read_dir,filename); %for color other sensor
    [img_width,img_height]=size(Ref_img);
    
    Ic{1,iter} = InterpSonyImage(Ref_img); 
    
    r = Ref_img(1:2:end,2:2:end) ./ Ref_img(2:2:end,1:2:end); 
    r_ = [r,flipud(r),r;fliplr(r),r,fliplr(r);r,flipud(r),r];
    h = fspecial('gaussian', [60 60], 30);
    rf = imfilter(r_, h);
%     rf = rf(1233:1232*2, 1641:1640*2); %this figures for half of image size// need to check
    rf = rf(img_width/2+1:img_width, img_height/2+1:img_height); %this figures for half of image size// need to check
    
%     c=Ref_img(2:2:end,1:2:end);

    Ref_img(2:2:end,1:2:end) =  Ref_img(2:2:end,1:2:end) .* rf; 

    If{1,iter} = InterpSonyImage(Ref_img); 
end


avg_holo = sperm_stack_avg_dk2(If,length_frame_stack);

Nb_ROI = divide_ROI^2;
wid_size = round(img_width/divide_ROI)-1;
hei_size = round(img_height/divide_ROI)-1;

figure(1568)
subplot(1,3,1)
imagesc(avg_holo)
title('Initial raw image')
set(gcf,'Color',[1 1 1]);

bandwidth_bak = 20;    %16 <--before basler sensor// 90<--pi sensor// bandwidth for lowpassfilter background // you need to change
dark_bias = min(avg_holo(:));
ori_level = mean2(avg_holo);
[avg_holo,myimg_filtered] = intensitycompensation(avg_holo,bandwidth_bak,bandwidth_bak, dark_bias);
avg_holo = avg_holo/mean2(avg_holo) * ori_level; %apply law pass filter to average hologram image

avg_holo(1:3,:)=mean(avg_holo(4,:));
avg_holo(end-3:end,:)=mean(avg_holo(end-4,:));
avg_holo(:,1:3)=mean(avg_holo(:,4));
avg_holo(:,end-3:end)=mean(avg_holo(:,end-4));

figure(1568)
subplot(1,3,2)
imagesc(myimg_filtered)
title('background image (low pass filtered)')
subplot(1,3,3)
imagesc(avg_holo)
title('Averaged filtered image')

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auto focus
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FLAG_AF

z2_l = 325e-6; % m 400 ~ 550 before
z2_u = 425e-6; % m 
delta_z2 = 0.1e-6; % m 

tf_phase = 1;
lambda = 633e-9; % m // Justin's LED: 633nm //New fibered LED (Chris'): 624nm
pixel_size = 1.12e-6; % m
refidx = 1;

K=0.12; %fraction of added bands

avg_holo_expan=zeros(img_width+2*round(K*wid_size),img_height+2*round(K*hei_size));
avg_holo_expan(1:img_width,1:img_height)=avg_holo;
avg_holo_expan(1:img_width,end-img_height+1:end)=avg_holo;
avg_holo_expan(end-img_width+1:end,1:img_height)=avg_holo;
avg_holo_expan(end-img_width+1:end,end-img_height+1:end)=avg_holo;

avg_holo_expan(round(K*wid_size)+1:img_width+round(K*wid_size),round(K*hei_size)+1:img_height+round(K*hei_size))=avg_holo;
% figure(458),imagesc(avg_holo_expan)

vec_wid=NaN(divide_ROI,wid_size+2*round(K*wid_size));
vec_hei=NaN(divide_ROI,hei_size+2*round(K*hei_size));
for k_it=1:divide_ROI
vec_wid(k_it,:)=(k_it-1)*wid_size+1:k_it*wid_size+2*round(K*wid_size);
vec_hei(k_it,:)=(k_it-1)*hei_size+1:k_it*hei_size+2*round(K*hei_size);
end


crop_img=NaN(wid_size+2*round(K*wid_size),hei_size+2*round(K*hei_size),Nb_ROI);
z0=NaN(1,Nb_ROI);
incr=0;
for i_c=1:divide_ROI
    for j_c=1:divide_ROI
incr=incr+1;
eval(['crop_img(:,:,incr)=avg_holo_expan(vec_wid(' num2str(i_c) ',:),vec_hei(' num2str(j_c) ',:));']);
z0(incr) = af_quick_v5(crop_img(:,:,incr), [z2_l*1e6 z2_u*1e6], pixel_size*1e6, 1, lambda*1e9, 2, 'Gini', 0, 1, [], 10, 0.01, 40, true, true );
% z(incr) = af_quick_v5(crop_img{1,image_number},[300 600], pixelsize, 1,wavelength,0, 'GiniOfGrad', 1, -1, [], 10, 0.1, 40, true, true );
% z(incr) = Auto_Focus(crop_img(:,:,incr),z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1); % m
title(['depth is ', num2str(z0(incr))]); 
zp=['Crop field ' num2str(incr) '/' num2str(Nb_ROI) ': depth is ',num2str(z0(incr)) 'um'];
disp(zp)
    end
end

% z = AF_chooseROI(avg_holo,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,1,1); % m

%z2 = Auto_Focus(avg_holo,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,1,1);

incr3=0;
mat_z=NaN(divide_ROI);
for i_r=1:divide_ROI
    for j_r=1:divide_ROI
            incr3=incr3+1;
            mat_z(i_r,j_r)=z0(incr3);
    end
end

figure(4444),imshow(mat_z*1e6,[])
title('Distribution of z2 distances (um)')
truesize(gcf,[1000 1000])
colorbar
axis equal
set(gcf,'Color',[1 1 1]);

end

if FLAG_plane

    
hori_zvec=nanmean(mat_z,1);  
vert_zvec=nanmean(mat_z,2);
fit_hori = robustfit(1:divide_ROI,hori_zvec);
fit_vert = robustfit(1:divide_ROI,vert_zvec);
    
[Xin,Yin] = meshgrid(1:divide_ROI);
plane_mat = K_fac*fit_hori(2)*Xin + K_fac*fit_vert(2)*Yin + mean([fit_hori(1) fit_vert(1)]);
mat_z_plot = inpaint_nans(mat_z);

figure(11212),imshow(plane_mat,[])
truesize(gcf,[1000 1000])
colorbar
axis equal
set(gcf,'Color',[1 1 1]);

figure(7788),surf(plane_mat*1e-6)
set(gcf,'Color',[1 1 1]);
set(gca,'NextPlot','add');
figure(7788),surf(mat_z_plot*1e-6)

incr4=0;
z=NaN(1,divide_ROI);
for i_r=1:divide_ROI
    for j_r=1:divide_ROI
            incr4=incr4+1;
            z(incr4)=plane_mat(i_r,j_r);
    end
end

elseif FLAG_mean_plane
    
z(:,:)=mean(z0(:));

else
z=z0;

end






if FLAG_recons

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reconstruction
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pixelsize=pixel_size*1e6; % um
interp_factor=1; %shirinking images 
lambda_V=lambda*1e9; %nm
n_glass=1;

z2_list=z*1e6; % um
% z2_list=745; % um
% z(:)=385;

rec_r=NaN(wid_size+2*round(K*wid_size),hei_size+2*round(K*hei_size));
for k_p=1:Nb_ROI
% rec_r= Prop_SSA(crop_img(:,:,k_p),pixelsize/interp_factor,pixelsize/interp_factor,-z2_list(k_p),lambda_V/1000/n_glass);
% eval(['rec_r(:,:,' num2str(k_p) ')=Prop_SSA(crop_img(:,:,' num2str(k_p) '),pixelsize/interp_factor,pixelsize/interp_factor,-z2_list(' num2str(k_p) '),lambda_V/1000/n_glass);']);
eval(['rec_r(:,:,' num2str(k_p) ') = Propagate(crop_img(:,:,' num2str(k_p) '), pixelsize, 1, lambda_V, z(' num2str(k_p) '), true, false, true);']);
% figure(9452),imshow(angle(rec_r(:,:,k_p)),[])
figure(9452),imshow(abs(rec_r(:,:,k_p)),[])
end

rec_r_crop=rec_r(round(K*wid_size):wid_size+round(K*wid_size),round(K*hei_size):hei_size+round(K*hei_size),:);
% rec_r_crop=rec_r;

%Assembling of the cropped ROIs:
incr2=0;
for i_r=1:divide_ROI
    incr2=incr2+1;
    eval(['conc_mat_' num2str(i_r) '=rec_r_crop(:,:,' num2str((i_r-1)*divide_ROI+1) ');']);
    for j_r=2:divide_ROI
        eval(['conc_mat_' num2str(i_r) '=[conc_mat_' num2str(i_r) '(:,:,1) rec_r_crop(:,:,' num2str(divide_ROI*(i_r-1)+j_r) ')];']);
    end
end

conc_mat_final = conc_mat_1;
for i_r2=2:divide_ROI
eval(['conc_mat_final = [conc_mat_final ; conc_mat_' num2str(i_r2) '];']);
end

figure(9452),imshow(angle(conc_mat_final),[])
title('Phase reconstructed image')
set(gcf,'Color',[1 1 1]);
figure(6527),imshow(abs(conc_mat_final),[])
title('Amplitude reconstructed image')
set(gcf,'Color',[1 1 1]);


end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% manual
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FLAG_manual

    for z=370:390
    
    z2_list=z;
    
    rec_r= Prop_SSA(avg_holo,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list,lambda_V/1000/n_glass);

    figure,imshow(angle(rec_r),[])       
    set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
    title(['depth is ', num2str(z2_list)]);


    end

z2_list=378;
rec_r= Prop_SSA(avg_holo,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list,lambda_V/1000/n_glass);
close all
figure,imshow(angle(real(rec_r)),[]) 
set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
title(['real part depth is ', num2str(z2_list)]);
figure,imshow(abs(rec_r),[]) 
set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
title(['abs part depth is ', num2str(z2_list)]);
figure,imshow(angle(rec_r),[])       
set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
title(['angle part depth is ', num2str(z2_list)]);
figure,imshow(imag(rec_r),[])  
set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
title(['imag part depth is ', num2str(z2_list)]);
figure,imshow(real(rec_r),[])  
set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
title(['real part depth is ', num2str(z2_list)]);
figure,imshow(abs(imag(rec_r)),[])  
set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
title(['abs(img) part depth is ', num2str(z2_list)]);
figure,imshow(abs(real(rec_r)),[])  
set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
title(['abs(real) part depth is ', num2str(z2_list)]);

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








