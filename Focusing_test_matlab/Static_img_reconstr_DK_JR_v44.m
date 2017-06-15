


%% max get rid of grid pattern

FLAG_select=1;
FLAG_pic_extrac=1;
FLAG_AF=0;
FLAG_recons=0;
FLAG_manual=0;

divide_ROI = 3; %Number of lines/columns that divide the FOV: corresponds to the following number of ROIs: divide_ROI^2


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Image Selection
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FLAG_select
selectFiles_tiff;
%selectFiles_jpg;
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
    
    img_G{1,iter} = deBayer_RPi_v2(Ref_img,'G');
    
    raw{1,iter} = Ref_img;
    r0{1,iter} = Ref_img(1:2:end,2:2:end);
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
avg_holog = sperm_stack_avg_dk2(img_G,length_frame_stack);
avg_holo_r = sperm_stack_avg_dk2(raw,length_frame_stack);

Nb_ROI = divide_ROI^2;
wid_size = round(img_width/divide_ROI);
hei_size = round(img_height/divide_ROI);

figure(1568)
subplot(1,3,1)
imagesc(avg_holo)
title('Initial raw image')
set(gcf,'Color',[1 1 1]);

avg_holo = sperm_stack_avg_dk2(If,length_frame_stack);
bandwidth_bak = 20;    %16 <--before basler sensor// 90<--pi sensor// bandwidth for lowpassfilter background // you need to change
dark_bias = min(avg_holo(:));
ori_level = mean2(avg_holo);
[avg_holo,myimg_filtered] = intensitycompensation(avg_holo,bandwidth_bak,bandwidth_bak, dark_bias);
figure,imagesc(myimg_filtered)
avg_holo = avg_holo/mean2(avg_holo) * ori_level; %apply law pass filter to average hologram image

figure(1568)
subplot(1,3,2)
imagesc(myimg_filtered)
title('background image (low pass filtered)')
subplot(1,3,3)
imagesc(avg_holo)
title('Averaged filtered image')

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auto_focus Chris
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FLAG_AF

z2_l = 300e-6; % m 400 ~ 550 before
z2_u = 500e-6; % m 
delta_z2 = 0.1e-6; % m 

tf_phase = 1;
lambda = 633e-9; % m // Justin's LED: 633nm //New fibered LED (Chris'): 624nm
pixel_size = 1.12e-6; % m
refidx = 1;

vec_wid=NaN(divide_ROI,wid_size);
vec_hei=NaN(divide_ROI,hei_size);
for k_it=1:divide_ROI
vec_wid(k_it,:)=(k_it-1)*wid_size+1:k_it*wid_size;
vec_hei(k_it,:)=(k_it-1)*hei_size+1:k_it*hei_size;
end


crop_img=NaN(wid_size,hei_size);
z=NaN(1,Nb_ROI);
incr=0;
pos_mat=NaN(Nb_ROI,2);
for i_c=1:divide_ROI
    for j_c=1:divide_ROI
incr=incr+1;
eval(['crop_img(:,:,incr)=avg_holo(vec_wid(' num2str(i_c) ',:),vec_hei(' num2str(j_c) ',:));']);
pos_mat(incr,:)=[i_c j_c];
% z(incr) = af_quick_v5(crop_img(:,:,incr), [z2_l*1e6 z2_u*1e6], pixel_size*1e6, 1, lambda*1e9, 0, 'Tamura', 1, -1, [], 10, 0.1, 40, true, true );
z(incr) = Auto_Focus(crop_img(:,:,incr),z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,[],1,1); % m
title(['depth is ', num2str(z(incr)*1e6)]); 
zp=['Crop field ' num2str(incr) '/' num2str(Nb_ROI) ': depth is ',num2str(z(incr)*1e6) 'um'];
disp(zp)
    end
end

% z = AF_chooseROI(avg_holo,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,1,1); % m

%z2 = Auto_Focus(avg_holo,z2_l,z2_u,delta_z2,pixel_size,refidx,lambda,tf_phase,1,1);

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

for k_p=1:Nb_ROI
rec_r= Prop_SSA(crop_img(:,:,k_p),pixelsize/interp_factor,pixelsize/interp_factor,-z2_list(k_p),lambda_V/1000/n_glass);
eval(['rec_r(:,:,' num2str(k_p) ')=Prop_SSA(crop_img(:,:,' num2str(k_p) '),pixelsize/interp_factor,pixelsize/interp_factor,-z2_list(' num2str(k_p) '),lambda_V/1000/n_glass);']);
% figure(9452),imshow(angle(rec_r(:,:,k_p)),[])
figure(9452),imshow(abs(rec_r(:,:,k_p)),[])
end


%Assembling of the cropped ROIs:
incr2=0;
for i_r=1:divide_ROI
    incr2=incr2+1;
    eval(['conc_mat_' num2str(i_r) '=rec_r(:,:,' num2str((i_r-1)*divide_ROI+1) ');']);
    for j_r=2:divide_ROI
        eval(['conc_mat_' num2str(i_r) '=[conc_mat_' num2str(i_r) '(:,:,1) rec_r(:,:,' num2str(divide_ROI*(i_r-1)+j_r) ')];']);
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

pixel_size = 1.12e-6; % m    
lambda = 527e-9; 
pixelsize=pixel_size*1e6; % um
interp_factor=1; %shirinking images 
lambda_V=lambda*1e9; %nm
n_glass=1;    

i=1;
count=0;
standr=[];
meani=[];
maxmin=[];
noiser=[];
close all
for z=360:380
    z2_list=z;
    rec_r{1,i}= Prop_SSA(avg_holo,pixelsize,pixelsize,-z2_list,lambda_V/1000/n_glass);
    figure,imshow(angle(rec_r{1,i}),[])
    %figure,imshow(angle(rec_r{1,i}),[])
    set(gca, 'XLim', [500,1200], 'YLim', [300,800])
    title(['depth is ', num2str(z2_list)]);
%     count=count+1;
%     standr(count)=std(c(:));
%     meani(count)=mean(c(:));
%     maxmin(count)=max(c(:))-min(c(:));
%     noiser(count)=standr(count) / meani(count);
%     xlabel(['Std: ' ,num2str(standr(count)),' Mean: ',num2str(meani(count)), ' MAX-MIn: ',num2str(maxmin(count)) ,'Noise: ',num2str(noiser(count))]); 
end

% figure,plot(1:count,standr)
% figure,plot(1:count,meani)
% figure,plot(1:count,maxmin)
% figure,plot(1:count,noiser)


    z2_list=386;
    rec_r{1,1}= Prop_SSA(avg_holo,pixelsize,pixelsize,-z2_list,lambda_V/1000/n_glass);
    figure,imshow(abs(abs(rec_r{1,1})),[]) 
    set(gca, 'XLim', [2170,2380], 'YLim', [830,970])
    
    z2_list=372;
    rec_r{1,1}= Prop_SSA(avg_holo,pixelsize,pixelsize,-z2_list,lambda_V/1000/n_glass);
    figure,imshow(abs(real(rec_r{1,1})),[]) 
    set(gca, 'XLim', [500,1200], 'YLim', [300,800])

hold on
plot([2300,2300 + 44.6429],[950,950],'Color','w','LineWidth',3);
text(2317,955, '50 um','Color','w','FontSize',30);    
    
avgimg=sperm_stack_avg_dk2(img_G,iter);
figure,imshow(img_G,[])   
z2_list=378.39;
recavg= Prop_SSA(avgimg,pixelsize,pixelsize,-z2_list,lambda_V/1000/n_glass);
figure,imshow(abs(real(recavg)),[])     

cc=abs(real(recavg));
r = cc(1:2:end,2:2:end) ./ cc(2:2:end,1:2:end); 
r_ = [r,flipud(r),r;fliplr(r),r,fliplr(r);r,flipud(r),r];
h = fspecial('gaussian', [60 60], 30);
rf = imfilter(r_, h);
rf = rf(img_width/2+1:img_width, img_height/2+1:img_height); %this figures for half of image size// need to check
cc(2:2:end,1:2:end) =  cc(2:2:end,1:2:end) .* rf; 
gg = InterpSonyImage(cc); 
figure,imshow(gg,[]) 


bckimg=sperm_bck_sub(img_G{1,1},avgimg);
z2_list=378.39;
recbck= Prop_SSA(bckimg,pixelsize,pixelsize,-z2_list,lambda_V/1000/n_glass); 
figure,imshow(recbck,[])     

% z2_list=378;
% rec_r= Prop_SSA(avg_holo,pixelsize/interp_factor,pixelsize/interp_factor,-z2_list,lambda_V/1000/n_glass);
% close all
% figure,imshow(angle(real(rec_r)),[]) 
% set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
% title(['real part depth is ', num2str(z2_list)]);
% figure,imshow(abs(rec_r),[]) 
% set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
% title(['abs part depth is ', num2str(z2_list)]);
% figure,imshow(angle(rec_r),[])       
% set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
% title(['angle part depth is ', num2str(z2_list)]);
% figure,imshow(imag(rec_r),[])  
% set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
% title(['imag part depth is ', num2str(z2_list)]);
% figure,imshow(real(rec_r),[])  
% set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
% title(['real part depth is ', num2str(z2_list)]);
% figure,imshow(abs(imag(rec_r)),[])  
% set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
% title(['abs(img) part depth is ', num2str(z2_list)]);
% figure,imshow(abs(real(rec_r)),[])  
% set(gca, 'XLim', [500,1400], 'YLim', [700,1400])
% title(['abs(real) part depth is ', num2str(z2_list)]);

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








