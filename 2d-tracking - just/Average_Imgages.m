% Create Image Stack
% Created by Wei LUO @ EE, UCLA
% Date (MM-DD-YYYY): 05-19-2013
% version: 1.1


% this function gets the average of a stack of images

function [averaged_img, num_avg, valid_idx] =  Average_Imgages(img_stack, avg_idx)

% INPUTS:
% img_stack: the stack of the images
% avg_idx: indices of the images that is used for getting average

% OUTPUTS:
% avg_img: averaged image
% num_avg: number of images that are used for average
% valid_idx: the valid indices of images that are used for average

% code starts here

% validate the images
% get the actual number of images
[num_row, num_col, num_img] = size(img_stack);

if (nargin<2)
	avg_idx = 1:num_img;
end

temp_flag = (avg_idx>0) .* (avg_idx<=num_img);
valid_idx = find(temp_flag);

num_avg = length(valid_idx);

averaged_img = zeros(num_row, num_col);

for i = 1:num_avg
	averaged_img = averaged_img + img_stack(:,:,valid_idx(i));
end

if (num_avg>0)
	averaged_img = averaged_img/num_avg;
end


if(length(num_avg)==0)
    averaged_img = zeros(size(img_stack,1), size(img_stack,2));
end

end