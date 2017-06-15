selectFiles_tiff;

read_dir = path_files;
length_frame_stack = length(file_1);

for iter=1:length_frame_stack
    
    filename=char(file_1(:,iter));
    [Ref_img] = read_img_wrapper_dk(read_dir,filename); %for color other sensor
    raw_stack{1,iter} = Ref_img;
    [img_width,img_height]=size(Ref_img);

    img_G{1,iter} = deBayer_RPi_v2(Ref_img,'G');
    
    Ic{1,iter} = InterpSonyImage(Ref_img); 
    
    r0{1,iter} = Ref_img(1:2:end,2:2:end);
    r = Ref_img(1:2:end,2:2:end) ./ Ref_img(2:2:end,1:2:end); 
    r_ = [r,flipud(r),r;fliplr(r),r,fliplr(r);r,flipud(r),r];
    h = fspecial('gaussian', [60 60], 30);
    rf = imfilter(r_, h);
    rf = rf(img_width/2+1:img_width, img_height/2+1:img_height); %this figures for half of image size// need to check
    
%     %look later
    bandwidth_bak = 18;    %16 <--before basler sensor// 90<--pi sensor// bandwidth for lowpassfilter background // you need to change
    dark_bias = min(img_G{1,iter}(:));
    ori_level = mean2(img_G{1,iter});
    [img_G{1,iter},myimg_filtered] = intensitycompensation(img_G{1,iter},bandwidth_bak,bandwidth_bak, dark_bias);
    img_G{1,iter} = img_G{1,iter}/mean2(img_G{1,iter}) * ori_level; %apply law pass filter to average hologram image
    
%     ROI_mask = roipoly(r0);
%     r0_crop = r0(930:1120,600:1000);
%     noise_rat = std(r0_crop(:))/mean(r0_crop(:))
    
    
    Ref_img(2:2:end,1:2:end) =  Ref_img(2:2:end,1:2:end) .* rf; 
    If{1,iter} = InterpSonyImage(Ref_img); 
 
%     figure,
%     imshow(If{1,iter},[])
%     title(['Green interporation ', num2str(iter)])
%     axis equal;
%     colorbar;
%     caxis([0 5])
    
%     bandwidth_bak = 18;    %16 <--before basler sensor// 90<--pi sensor// bandwidth for lowpassfilter background // you need to change
%     dark_bias = min(If{1,iter}(:));
%     ori_level = mean2(If{1,iter});
%     [If{1,iter},myimg_filtered] = intensitycompensation(If{1,iter},bandwidth_bak,bandwidth_bak, dark_bias);
%     avg_holo1{1,iter} = If{1,iter}/mean2(If{1,iter}) * ori_level; %apply law pass filter to average hologram image

%     hFig = figure(100+iter);
%     set(hFig, 'Position', [1 1 1700 700])
%     imagesc(r0{1,iter});
%     axis equal;
%     caxis([0 1023]);
%     colorbar;
%     set(gcf,'Color',[1 1 1]);
%     title(filename,'interpreter','none');

    image=img_G{1,iter}(1160:1280, 1880:2000);

    normalized=image-min(image(:));
    normalized=normalized/max(normalized(:));

    [y,x]=imhist(normalized);

    result = entropy(y);

    figure,
    imshow(img_G{1,iter},[]);
%     title(filename,'interpreter','none');
%     noise=std(img_G{1,iter}(:)) / mean(img_G{1,iter}(:));
%     xlabel(['Noise Ratio is ' ,num2str(result)]); 
%     set(gca, 'XLim', [1000,1400], 'YLim', [700,1000])
    
%     figure,
%     imshow(avg_holo1{1,iter},[])
%     title(['interporation filtered image',num2str(iter)])
%     axis equal;
%     colorbar;
%     caxis([0 1000])
end

% image_pro=raw_stack{1,1};
% 
% rr=image_pro(1:2:end,1:2:end);
% mean(rr(:))
% 
% rg1=image_pro(1:2:end,2:2:end);
% mean(rg1(:))
% rg2=image_pro(2:2:end,1:2:end);
% mean(rg2(:))
% rgt=[rg1; rg2];
% mean(rgt(:))
% 
% rb=image_pro(2:2:end,2:2:end);
% mean(rb(:))
% 
% rrch=rr*( mean(rgt(:))/ mean(rr(:)) );
% rbch=rb*( mean(rgt(:))/ mean(rb(:)) );
% 
% [ht,wd]=size(raw_stack{1,1});
% 
% image_ch=zeros(size(raw_stack{1,1}));
% for i=1:100
%     for j=1:100
%         if mod(i,2) == 1 && mod(j,2) == 1
%             image_ch(i,j)=raw_stack{1,1}(i,j) *( 142.0/ 91.1 ) ;
%         elseif mod(i,2) == 1 && mod(j,2) == 0
%             image_ch(i,j)=raw_stack{1,1}(i,j);
%         elseif mod(i,2) == 0 && mod(j,2) == 1
%             image_ch(i,j)=raw_stack{1,1}(i,j);
%         elseif mod(i,2) == 0 && mod(j,2) == 0     
%             image_ch(i,j)=raw_stack{1,1}(i,j) *( 76.9/ 91.1 );
%         end 
%     end
% end