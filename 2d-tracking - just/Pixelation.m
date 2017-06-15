% Created by Wei Luo @ EE, UCLA
% Date(MM-DD-YYYY): 02-23-2014

% ver: 1.1: commented out: put the last row and last column to be the same as the second last 
% this function creates a pixel lated image using a specified pixel function

function [output_img, down_num_row, down_num_col] = Pixelation(img_in, SR_factor_x, SR_factor_y, PF)

% variables
% img_in: input images
% SR_factor_x, SR_factor_y: super resolution factors
% PF: pixel function,

% if the pixel function is not specified, create a flat pixel function
if (nargin<4)
	PF = ones(SR_factor_y, SR_factor_x);
end

% normalize the pixel function, total summation should be one
PF = PF / sum(PF(:));

[row_old, col_old] = size(img_in);

row_new = row_old - mod(row_old, SR_factor_y);
col_new = col_old - mod(col_old, SR_factor_x);

% get cropped image with rounded size
img_in_new = img_in(1:row_new, 1:col_new);

% consider the pixel function as a delta function
row_bias = floor(SR_factor_y/2) + 1;
col_bias = floor(SR_factor_x/2) + 1;

row_range = row_bias:SR_factor_y:row_new;
col_range = col_bias:SR_factor_x:col_new;

% convolve with the pixel function
new_HiRes_I = conv2(img_in_new, PF, 'same');

output_img = new_HiRes_I(row_range, col_range);

% put the last row and last column to be the same as the second last
% output_img(end,:) = output_img(end-1,:);
% output_img(:,end) = output_img(:,end-1);


[down_num_row, down_num_col] = size(output_img);

end
