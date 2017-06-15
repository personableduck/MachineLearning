

selectFiles_tiff;

length_frame_stack = length(file_1);

for iter=1:length_frame_stack

filename=char(file_1(:,iter));

[Ref_img] = read_img_wrapper_dk(path_files,filename); %for color other sensor

[img_width,img_height]=size(Ref_img);

r0 = Ref_img(1:2:end,2:2:end);
r = Ref_img(1:2:end,2:2:end) ./ Ref_img(2:2:end,1:2:end);
r_ = [r,flipud(r),r;fliplr(r),r,fliplr(r);r,flipud(r),r];
h = fspecial('gaussian', [60 60], 30);
rf = imfilter(r_, h);
rf = rf(img_width/2+1:img_width, img_height/2+1:img_height); %this figures for half of image size// need to check

Ref_img(2:2:end,1:2:end) =  Ref_img(2:2:end,1:2:end) .* rf;

If = InterpSonyImage(Ref_img);

hFig = figure(100+iter);
set(hFig, 'Position', [1 1 1700 700])
imagesc(r0(300:525,725:1025))
axis equal
caxis([0 1023])
colorbar
set(gcf,'Color',[1 1 1]);
title(filename,'interpreter','none')

end