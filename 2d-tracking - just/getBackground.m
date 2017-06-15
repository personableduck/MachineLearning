% get background
%% Created by Wei LUO @ EE, UCLA
% Date (MM-DD-YYYY): 29-08-2014
% version: 1.0

function [background, est_center] = getBackground(img_in, center, radius)

% variables: 
% img_in: input image center: in the form of [row, col], the
% center of the region that is used for background estimation
% radius: size of the region that is used for background estimation

if (nargin<2)
    msg = 'Choose the center of DC (background) estimation';
    pixeltype = 'NONE';
	[row_center col_center] = ChooseCenter_Msg(img_in, pixeltype, msg);
else
	row_center = center(1);
	col_center = center(2);
end
	row_range = (row_center-radius):(row_center+radius);
	col_range = (col_center-radius):(col_center+radius);

	DC_img = img_in(row_range, col_range);
	background = mean2(DC_img);
	est_center = [row_center col_center];


end