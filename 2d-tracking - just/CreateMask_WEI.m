% Create mask
% Created by Wei LUO @ EE, UCLA
% Date (MM-DD-YYYY): 05-19-2013
% version: 1.1

% this code creates a stack of images that will be used for tail reconstruction
% clear all
% close all
% clc


function [ mask ] = CreateMask(input_filed, th_method, thresholds, fill_holes, num_erode, num_dilate, dialatefactor, num_blur)

% input_filed = m_object;

% threshold method
if (nargin < 2)
	th_method = 'amplitude'; % can be 'phase' or 'amplitude' 
end

% step one thresholding
switch th_method
	case 'phase'
		m_img = angle(input_filed); % fprintf('Phase! \n')
	otherwise
		m_img = abs(input_filed); % fprintf('Amplitude! \n')
end

% thresholds: lower boudary and upper boundaries
if (nargin<3)
	thresholds = [(mean(m_img(:))), max(m_img(:)) ];
end

th_min = min(thresholds); 0.52;
th_max = max(thresholds); 0.6;

% flag for filling holes
if(nargin<4)
    fill_holes = 1;
end
% number of erosion
if(nargin<5)
    num_erode = 1;
end

% number of dialations
if(nargin<6)
    num_dilate = 1;
end

if(nargin<7)
    dialatefactor = 1;
end

% number of guassian blurs
if(nargin<8)
    num_blur = 1;
end

%% now start creating the mask

% thresholding
mask = (m_img>th_min).* (m_img<th_max);


% fill holes
if (fill_holes>0)
    mask = imfill(mask,'holes');
end



% erosion
se = strel('disk',1);
for i = 1:num_erode
   mask = imerode(mask,se);
end




% dialate the objects
for i = 1:num_dilate
	mask = imdilate(mask, strel('disk',dialatefactor));
end

% figure
% subplot(1,2,1)
% imagesc(m_img); axis image; colormap gray; title('image')
% subplot(1,2,2)
% imagesc(mask); axis image; colormap gray; title('dialate mask')
% pause

% smooth the boudaries
for i = 1:num_blur
    mask = SmoothPatternBoundary(double(mask), 3);
end
% label_obj = bwconncomp(mask,8);
% imagesc(mask); axis image; colormap gray; title('Mask')

end