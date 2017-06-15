% %Reconstruction parameters
%%% DEPTHS
%justin's data:
%     IUI-2466x5-1_1000umG - 1765
%     IUI-2466x5-1_150umG  - 1000
%Data presented 10.29 - 2400
%2fps, 1 and 2 - 2700
%5fps - 3000
%10fps, 1 and 2 - 2700
z2_list = [2061.328125]; %802.8846
num_z2_slices = length(z2_list);

rotation = 0;
interp_factor = 1;

%Justin's reconst parameters
pixelsize = 2.2;    % (um)
lambda_V = 645;
n_glass = 1.332; %justin 1.5
% n_glass = 1.332; %current setup

%Localization parameters (in pixels)
distanceThreshold = 30; %in microns <--------------------- need to control dpend on frame rate
r_cencal = 10;       % radius of the area for calculating centroids (px)
tolerance_z = 50;    % for enclosing searching within 15 z2 micron

%convert to pixels
gridSize = 1/pixelsize;
distanceThreshold = distanceThreshold*gridSize; 

%Minimum number of frames tracked for a particular sperm head to be of
%interest.
num_frames_min_thresh = 3;

%Point of interest
POI_specified = 0;

% % Select frames to process
selectFiles;
read_dir = path_files;

single_image = ~iscell(file_1);

% SINGLE IMAGE
if (single_image)
 length_frame_stack = 1; % for single image
 [Ref_img, pixelsize] = read_img_wrapper(read_dir,filename,rotation);
 height = size(Ref_img,1)*interp_factor;
 width = size(Ref_img,2)*interp_factor;
 recon_stack = zeros(height,width,1);
 holo_stack = zeros(size(Ref_img,1), size(Ref_img,2), 1);
%MULTIPLE IMAGES
else
 length_frame_stack = length(file_1); % for multiple images
 avg_holo = sperm_stack_avg(file_1,read_dir); %form average hologram for background subtraction
 height = size(avg_holo,1)*interp_factor;
 width = size(avg_holo,2)*interp_factor;
 recon_stack = zeros(height,width,length_frame_stack);
 recon_stack_diff = zeros(height,width,length_frame_stack-1);
 holo_stack = zeros(size(avg_holo,1),size(avg_holo,2),length_frame_stack);
end

recon_avg_holo = zeros(height,width);
holo_stack_avg = zeros(size(avg_holo,1), size(avg_holo,2));

%----- law filter pass function ( make equal intensity on image)
% bandwidth_bak = 16;    
% dark_bias = min(avg_holo(:));
% ori_level = mean2(avg_holo);
% avg_holo = intensitycompensation(avg_holo,bandwidth_bak,bandwidth_bak, dark_bias);
% avg_holo = avg_holo/mean2(avg_holo) * ori_level; %apply law pass filter to average hologram image
%-----------------------------------------------------------------------------------
% 
% rec_avg = Prop_SSA(avg_holo,pixelsize/interp_factor,pixelsize/interp_factor,...
%                     -z2_list(1),lambda_V/1000/n_glass);
% recon_avg_holo = abs(rec_avg); % make an average hologram image(leave only stationary objects)
% 
% threshold_sum = mean2(recon_avg_holo)*0.80; %Because low pass filter, we need to control the (*80) figure.
% objMask_sum = (recon_avg_holo<threshold_sum); % intensity filter function, get rid of background.
% 
% filter_size = bwareafilt(objMask_sum,[1 8]); % filter pixel size, delte big pixel size, uper version of matlab 2014b 


%----jeff's size filter
% 
% max_sperm_size = 8;
% for i = 1:size(prop,1)
%    if (objMask_sum(i,1) <= max_sperm_size)
%        sperm_count = sperm_count + 1;
%        sperm_position(sperm_count, :) = region_centroid(i,:);
%    end
% end

%-------Jeff's 9 pixel compare filter 
% 
% result = zeros(height, width);
% for i = 3:height-2
%    for j = 3:width-2
%        y = j;
%         x = i;
%         inner_avg = (abs(recon_avg_holo(x,y)) + abs(recon_avg_holo(x+1,y)) + abs(recon_avg_holo(x-1,y))+ abs(recon_avg_holo(x,y+1)) + abs(recon_avg_holo(x,y-1)) + abs(recon_avg_holo(x+1,y+1))+ abs(recon_avg_holo(x+1,y-1)) + abs(recon_avg_holo(x-1,y+1)) + abs(recon_avg_holo(x-1,y-1)))/9;
%          outer_avg = (abs(recon_avg_holo(x+2,y+2)) + abs(recon_avg_holo(x+1,y+2)) + abs(recon_avg_holo(x,y+2))+ abs(recon_avg_holo(x-1,y+2)) + abs(recon_avg_holo(x-2,y+2)) + abs(recon_avg_holo(x-2,y+1))+ abs(recon_avg_holo(x-2,y)) + abs(recon_avg_holo(x-2,y-1)) + abs(recon_avg_holo(x-2,y-2)) + abs(recon_avg_holo(x-1,y-2)) + abs(recon_avg_holo(x,y-2)) + abs(recon_avg_holo(x+1,y-2))+ abs(recon_avg_holo(x+2,y-2)) + abs(recon_avg_holo(x+2,y-1)) + abs(recon_avg_holo(x+2,y))+ abs(recon_avg_holo(x+2,y+1)) )/16;
%         result(x,y)=outer_avg-inner_avg;
%    end
% end
% recon_avg_holo_BW = result > mean2(abs(result))*4;%0.055;
% filter_size_2 = recon_avg_holo_BW & filter_size;
% 
% intenMask_avg = filter_size_2.*recon_avg_holo;
% objects_cont = regionprops(filter_size_2, intenMask_avg, 'Area', 'Centroid', 'BoundingBox', 'MinIntensity');
% 
% centroidOrig_avg=cell2mat({objects_cont.Centroid}'); %find stationary sperm's coordinate, number

%-----finish finding stationary sperms' function 
%-----start finding moving objects 

coor_from = cell(1,length_frame_stack-1);
coor_to = cell(1,length_frame_stack-1);

for iter = 1:length_frame_stack
    
    %single file
    if (single_image)
        filename=char(file_1);
    %multiple files
    else
        filename=char(file_1(:,iter));
    end
    
    [Ref_img, pixelsize] = read_img_wrapper2(read_dir,filename,rotation);
    
    if( iter == 1)
        %For each frame
        I = Ref_img; 
        %background
        [ ~, bg ] = correct_illum( I, 30, 40, 0.05, 3 ); 
    end
    
    holo_stack(:,:,iter) = Ref_img;
    
    I= Ref_img;
    Ic = I ./ bg ./ 2; 
    Ref_img = Ic;
    
    %%%Avg Frame Substraction - Needs Fixing
%     Ref_img = sperm_bck_sub(Ref_img,avg_holo);   
%     recon_stack(:,:,iter) = rec;


%%% MODULE 1 - RECONSTRUCTION
    height = size(Ref_img,1)*interp_factor; 
    width = size(Ref_img,2)*interp_factor; 
    for k = 1:num_z2_slices
        %%Only 1 Angle
        %%Backpropagation  -- at predetermined z2 slices
        rec = Prop_SSA(Ref_img,pixelsize/interp_factor,pixelsize/interp_factor,...
                    -z2_list(k),lambda_V/1000/n_glass);
        
%         rec = abs(rec);
%         rec = angle(rec);
        %%Place image into stack
        recon_stack(:,:,iter) = rec;

    end

%%% MODULE 2 - DIFFERENCE HOLOGRAM SUBSTRACTION AND OBJECT DETECTION
    %If previous frame exists
    if(iter-1>0)
        diff_holo = holo_stack(:,:,iter) - holo_stack(:,:,iter-1);
        rec = Prop_SSA(diff_holo,pixelsize/interp_factor,pixelsize/interp_factor,...
                    -z2_list(k),lambda_V/1000/n_glass);
        rec = abs(rec);
        recon_stack_diff(:,:,iter-1) = rec;
        
        %Get points on image -- Currently, using simple thresholding
        threshold = mean2(rec)*4; %<--------------------- need to control dpend on resolution rate
        objMask = (rec>threshold);
%         objMask = bwareaopen(objMask,2);  %Added for noisy sample (ours, not justins) to remove 1 pixel values
        intenMask = objMask.*abs(recon_stack(:,:,iter));
        objects = regionprops(objMask, intenMask, 'Area', 'Centroid', 'BoundingBox', 'MinIntensity');
        
        %TODO: FILTER OUT BIG OBJECTS BASED ON SIZE
        
        %If mask returns empty, skip diff holo tracking for these two
        %frames
        if(isempty(objects))
            continue
        else
            centroidOrig = cell2mat({objects.Centroid}');
            minIntensities = cell2mat({objects.MinIntensity}');
            
            min_intensities_xyz(1,iter-1) = {minIntensities};
            xyz(1,iter-1)={centroidOrig}; % for IDL track, indicate whole moving objects coordinate before slice
            
            %For each point, find closest point, check distance
            %Take lower pix val to be from traveled, take higher pix val to be to traveled
            [centroidFrom,centroidTo] = MatchPoints(centroidOrig,distanceThreshold,minIntensities);
            
            coor_mid(iter-1) = {centroidOrig};
            coor_from(iter-1) = {centroidFrom};
            coor_to(iter-1) = {centroidTo};
        end    
        
    end

    
    %%OLD LOCALIZATION
    %%Cluster points in 2D     -- find and match points and determine optimal z2.
%     r_nearby = 2*r_cencal;
%     [peak_coor,peak_value,node_cluster] = getClusters(recon_stack,z2_list,r_nearby,tolerance_z);
%     
%     coor_fin{iter} = peak_coor';
end
disp('Reconstruction Finished!');

figure, imshow(abs(recon_stack(:,:,1)),[]) % first reconstructed image
%-------plot moving objects
hold on
for jkjk=1:iter-1
    hold on
    plot((xyz{1,jkjk}(:,1)),(xyz{1,jkjk}(:,2)),'g.')
end

sz_xyz=size(xyz);
sum_szxyz=0;

for loop_szxyz=1:sz_xyz(1,2)
    num_szxyz=size(xyz{1,loop_szxyz});
    sum_szxyz=sum_szxyz+num_szxyz(1,1);
end

avg_szxyz=sum_szxyz/sz_xyz(1,2)
%--------------------plot stationary sperms

% hold on
% plot(centroidOrig_avg(:,1),centroidOrig_avg(:,2),'y.')


% %% MODULE 3 - 2D TRACKING
% %%OLD 3D TRACKING
% % %%3D Head Tracking (TIME for 2D tracking case) -- Takes coordinates for detected head of sperm from each of
% % %%frame and performs 3d tracking
% tracksSimple = run_3d_head_tracking(length_frame_stack,coor_fin);
% 
% numFrames = [tracksSimple.end_frame] - [tracksSimple.start_frame];
% idxTracks = find(numFrames > num_frames_min_thresh);
% tst = tracksSimple(idxTracks);
% 
% numTracks = size(tracksSimple,1);
% goodTracks = zeros(numTracks,1);

%------------
% axis on
% set(gca, 'XLim', [1100,1400], 'YLim', [700, 1000])
% plot([1100; 1400], [700; 700], '--w',  [1100; 1400], [1000; 1000], '--w', 'LineWidth', 2)
% plot([1100; 1100], [700; 1000], '--w',  [1400; 1400], [700; 1000], '--w', 'LineWidth', 2)
% plot([1100; 1400], [700; 700], '--w',  [1100; 1400], [1000; 1000], '--w', 'LineWidth', 2)
% plot([1100; 1100], [700; 1000], '--w',  [1400; 1400], [700; 1000], '--w', 'LineWidth', 2)