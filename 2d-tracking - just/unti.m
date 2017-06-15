
holo_re=holo_stack(:,:,1);

h_rec = Prop_SSA(holo_re,pixelsize/interp_factor,pixelsize/interp_factor,...
                    -z2_list(k),lambda_V/1000/n_glass);
h_rec=abs(h_rec);

figure, imshow(h_rec,[])

%-----------------------------

[rows, columns, numberofColorChannels] = size(h_rec);
if numberofColorChannels > 1
  grayImage = rgb2gray(h_rec);
else
  grayImage  = h_rec;
end

binaryImage = grayImage >= 200;
numberOfWhitePixels = sum(binaryImage(:));



%----------------------------

threshold_h = mean(h_rec);
objMask_h = (h_rec>threshold_h);
intenMask_h = objMask_h.*abs(h_rec);
objects_h = regionprops(objMask_h, intenMask_h, 'Area', 'Centroid', 'BoundingBox', 'MinIntensity');
centroidOrig_h = cell2mat({objects_h.Centroid}');

figure
imshow(intenMask_h,[])

hold on
plot((centroidOrig_h(:,1)),(centroidOrig_h(:,2)),'.b')

%--------------------------------

sum_holo = holo_stack(:,:,1) - holo_stack(:,:,2);
sum_rec = Prop_SSA(sum_holo,pixelsize/interp_factor,pixelsize/interp_factor,...
                    -z2_list(1),lambda_V/1000/n_glass);
sum_rec = abs(sum_rec);
recon_stack_sum(:,:,1) = sum_rec;

threshold_s1=max(holo_stack(:,:,1));
threshold_s= mean2(sum_rec)*3;
objMask_s = (sum_rec>threshold_s);

intenMask_s = objMask_s.*abs(recon_stack_sum(:,:,1));
objects_s = regionprops(objMask_s, intenMask_s, 'Area', 'Centroid', 'BoundingBox', 'MinIntensity');

centroidOrig_s = cell2mat({objects_s.Centroid}');
minIntensities_s = cell2mat({objects_s.MinIntensity}');

figure
imshow(recon_stack_sum(:,:,1),[])

hold on
plot((centroidOrig_s(:,1)),(centroidOrig_s(:,2)),'.b')