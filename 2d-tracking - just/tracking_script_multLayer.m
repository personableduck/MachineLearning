% %Reconstruction parameters
%justin - 1765
%Data presented 10.29 - 2400
%2fps, 1 and 2 - 2700
%5fps - 3000
%10fps, 1 and 2 - 2700
z2_list = [1765];
num_z2_slices = length(z2_list);

rotation = 0;
interp_factor = 1;

%Justin's reconst parameters
pixelsize = 2.2;    % (um)
lambda_V = 645;
n_glass = 1.5; %justin 1.5
% n_glass = 1.332; %current setup

%Localization parameters (in pixels)
distanceThreshold = 30; %in microns
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

coor_from = cell(1,length_frame_stack-1);
coor_to = cell(1,length_frame_stack-1);
%For each frame (temporal)
for iter = 1:length_frame_stack
    
    %single file
    if (single_image)
        filename=char(file_1);
    %multiple files
    else
        filename=char(file_1(:,iter));
    end
    
    [Ref_img, pixelsize] = read_img_wrapper(read_dir,filename,rotation);
    holo_stack(:,:,iter) = Ref_img;
%     Ref_img = sperm_bck_sub(Ref_img,avg_holo);   %Still needs fixing
%     recon_stack(:,:,iter) = rec;

    height = size(Ref_img,1)*interp_factor; 
    width = size(Ref_img,2)*interp_factor; 
    recon_stack_ang1 = zeros(height,width,num_z2_slices);
    for k = 1:num_z2_slices
level = 0.065;
        %%Only 1 Angle
        %%Backpropagation  -- at predetermined z2 slices
        rec = Prop_SSA(Ref_img,pixelsize/interp_factor,pixelsize/interp_factor,...
                    -z2_list(k),lambda_V/1000/n_glass);
%         rec = abs(rec);
%         rec = angle(rec);
        %%Place image into stack
        recon_stack_ang1(:,:,k) = rec;
        recon_stack(:,:,iter) = rec;

    end
    
    %If previous frame exists
    if(iter-1>0)
        diff_holo = holo_stack(:,:,iter) - holo_stack(:,:,iter-1);
        rec = Prop_SSA(diff_holo,pixelsize/interp_factor,pixelsize/interp_factor,...
                    -z2_list(k),lambda_V/1000/n_glass);
        rec = abs(rec);
        recon_stack_diff(:,:,iter-1) = rec;
        
        %Get points on image -- Currently, using simple thresholding
        threshold = mean2(rec)*10;
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
            
            %For each point, find closest point, check distance
            %Take lower pix val to be from traveled, take higher pix val to be to traveled
            [centroidFrom,centroidTo] = MatchPoints(centroidOrig,distanceThreshold,minIntensities);
            
            coor_from(iter-1) = {centroidFrom};
            coor_to(iter-1) = {centroidTo};
        end    
        
    end

    
    %%OLD LOCALIZATION
    %%Cluster points in 2D     -- find and match points and determine optimal z2.
%     r_nearby = 2*r_cencal;
%     [peak_coor,peak_value,node_cluster] = getClusters(recon_stack_ang1,z2_list,r_nearby,tolerance_z);
% %     
    %%MAYBE COMBINE PEAK COORS HERE TO ACCOUNT FOR MULTIPLE SAME NODE
    %%CLUSTERS??
    
    %%Pick specific sperm to process -- if debugging and looking at only sperm
%     if(POI_specified == 1)
%         [currIdx, currDist] = find_min_dist(POI_x,POI_y,peak_coor(1,:),peak_coor(2,:));
%     end
%     
%     coor_fin{iter} = peak_coor';
end
disp('Reconstruction Finished!');



%%%OLD 3D TRACKING
% % %%3D Head Tracking (TIME for 2D tracking case) -- Takes coordinates for detected head of sperm from each of
% % %%frame and performs 3d tracking
% tracksSimple = run_3d_head_tracking(length_frame_stack,coor_fin);

% numFrames = [tracksSimple.end_frame] - [tracksSimple.start_frame];
% idxTracks = find(numFrames > num_frames_min_thresh);
% tst = tracksSimple(idxTracks);
% 
% numTracks = size(tracksSimple,1);
% goodTracks = zeros(numTracks,1);

