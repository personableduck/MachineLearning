%clear;%%%%%%%%%%%%%%
%depth
z2_list = [825];
num_z2_slices = length(z2_list);

rotation = 0;
interp_factor = 1;

%Justin's reconst parameters
pixelsize = 2.2;    % (um)
lambda_V = 645;
n_glass = 1.5; %justin 1.5
% n_glass = 1.332; %current setup

% % Select frames to process
selectFiles;
%%%%%%%%%%%%%%%%%%%%%%%%the following two lines occured in the end of
%%%%%%%%%%%%%%%%%%%%%%%%selectFIles
read_dir = path_files;
%single_image = ~iscell(file_1);

avg_holo = sperm_stack_avg(file_1,read_dir); %form average hologram for background subtraction

height = size(avg_holo,1)*interp_factor;
width = size(avg_holo,2)*interp_factor;

recon_avg_holo = zeros(height,width);%reconstructed average frame

%code from Wei to smooth the intensity using low pass filtering
bandwidth_bak = 16;    
dark_bias = min(avg_holo(:));
ori_level = mean2(avg_holo);
avg_holo = intensitycompensation(avg_holo,bandwidth_bak,bandwidth_bak, dark_bias);
avg_holo = avg_holo/mean2(avg_holo) * ori_level;

%do the reconstruction for avg_holo
rec = Prop_SSA(avg_holo,pixelsize/interp_factor,pixelsize/interp_factor,...
                    -z2_list(1),lambda_V/1000/n_glass);
recon_avg_holo = rec; 



%BW
%recon_avg_holo_BW = zeros(height,width);%reconstructed average frame in BW
% amplitude = abs(recon_avg_holo);
% 
% thresholdValue = mean2(amplitude)*0.75;
% recon_avg_holo_BW_2 = amplitude > thresholdValue;
%learning_sperm_recon_avg_1to9;

%--------

% intenMask = objMask.*abs(recon_stack(:,:,iter));
% objects = regionprops(objMask, intenMask, 'Area', 'Centroid', 'BoundingBox', 'MinIntensity');
% %_________

result = zeros(height, width);
for i = 3:height-2
   for j = 3:width-2
       y = j;
        x = i;
        inner_avg = (abs(recon_avg_holo(x,y)) + abs(recon_avg_holo(x+1,y)) + abs(recon_avg_holo(x-1,y))+ abs(recon_avg_holo(x,y+1)) + abs(recon_avg_holo(x,y-1)) + abs(recon_avg_holo(x+1,y+1))+ abs(recon_avg_holo(x+1,y-1)) + abs(recon_avg_holo(x-1,y+1)) + abs(recon_avg_holo(x-1,y-1)))/9;
         outer_avg = (abs(recon_avg_holo(x+2,y+2)) + abs(recon_avg_holo(x+1,y+2)) + abs(recon_avg_holo(x,y+2))+ abs(recon_avg_holo(x-1,y+2)) + abs(recon_avg_holo(x-2,y+2)) + abs(recon_avg_holo(x-2,y+1))+ abs(recon_avg_holo(x-2,y)) + abs(recon_avg_holo(x-2,y-1)) + abs(recon_avg_holo(x-2,y-2)) + abs(recon_avg_holo(x-1,y-2)) + abs(recon_avg_holo(x,y-2)) + abs(recon_avg_holo(x+1,y-2))+ abs(recon_avg_holo(x+2,y-2)) + abs(recon_avg_holo(x+2,y-1)) + abs(recon_avg_holo(x+2,y))+ abs(recon_avg_holo(x+2,y+1)) )/16;
        result(x,y)=outer_avg-inner_avg;
   end
end
recon_avg_holo_BW = result > mean2(abs(result))*6.1;%0.055;
recon_avg_holo_BW = recon_avg_holo_BW & objMask_sum;
%recon_avg_holo_BW = recon_avg_holo_BW_2;
figure
imshow(recon_avg_holo_BW,[]);

filter_size_2 = bwareafilt(recon_avg_holo_BW,[1 8]);
figure
imshow(filter_size_2,[]);

%get size and centroid of each region
prop = regionprops(recon_avg_holo_BW, 'Area', 'Centroid');

region_size = cell2mat({prop.Area}');
region_centroid = cell2mat({prop.Centroid}');

sperm_position = zeros(size(prop,1),2);
sperm_count = 0;

%threshold for sperm size
max_sperm_size = 8;
for i = 1:size(prop,1)
   if (region_size(i,1) <= max_sperm_size)
       sperm_count = sperm_count + 1;
       sperm_position(sperm_count, :) = region_centroid(i,:);
   end
end

%show reconstructed avg_holo
figure;
imshow(abs(recon_avg_holo),[]);
title('recon_avg_holo');
hold on
plot(sperm_position(1:sperm_count,1),sperm_position(1:sperm_count,2),'g.');

%show reconstructed avg_holo in BW
figure;
imshow(recon_avg_holo_BW,[]);
title('recon_avg_holo_BW');
hold on
plot(sperm_position(1:sperm_count,1),sperm_position(1:sperm_count,2),'g.');




