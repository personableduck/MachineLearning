%% Created by Wei LUO @ EE, UCLA
% Date (MM-DD-YYYY): 10-09-2015
% Version: 1.1
%%

close all
clc
clear

% this code subtract the background (i.e. averaged image from a specified window)

readFolder = 'H:\Data\151202_300FPS\200X\ROI7-1\';
% readFolder = '/Users/weiluo/Documents/Data/151202_300FPS/ROI7-1_tif/ROI7-1/';

writeFolder = 'H:\Data\151202_300FPS\200X\ROI7-1_avg\';
% writeFolder = '/Users/weiluo/Documents/Data/151202_300FPS/ROI7-1_tif/ROI7-1_avg/';


keyword = '*.tif';
% list all files
files = dir([readFolder keyword]);
cc = cellstr([{files(:).name}]);
[cs,index] = sort_nat(cc,'ascend');
readin_files = files(index);   % now the files are sorted by name
num_of_files = length(readin_files);	%get number of the files

window_size_pos = 300; % number of later frames that are used for avg 
window_size_neg = 300; % number of previous frames that are used for avg 

window_distance = 1; % distance of the window from the image

num_row = 301;
num_col = 800;

num_bits = 16;
%% start processing the data
if(exist(writeFolder, 'dir') == 0)
	mkdir(writeFolder);
end

figure
for i = 1:num_of_files
	m_distance = ((1:num_of_files) - i);
    m_avg_mask_pos = ((m_distance - window_distance)<window_size_pos) .* (m_distance>window_distance);
    m_avg_mask_neg = (abs(-m_distance - window_distance)<window_size_neg) .* (m_distance<(-window_distance));
    
	% m_avg_mask = (m_distance >= window_distance);
    m_avg_mask = ((m_avg_mask_pos + m_avg_mask_neg)>0);
	avg_indices = find(m_avg_mask);

	% average the images
	m_num_avg = length(avg_indices);
	avg_image = zeros(num_row, num_col);
	for j = 1:m_num_avg
		m_filename = [readFolder num2str(avg_indices(j)) '.tif'];
		m_temp = double(imread(m_filename));
		avg_image = avg_image + m_temp / mean2(m_temp); % remember to normalize the image before summing them up!
    end
    
    if(m_num_avg>1)
        avg_image = avg_image/mean2(avg_image);
    end
	% background subtraction
	m_readfilename = [readFolder num2str(i) '.tif'];
	m_temp = double(imread(m_readfilename));
    
    if(m_num_avg>1)
    	m_temp = m_temp/mean2(m_temp) - avg_image + 1; % remember to put the background image back        
    else
        m_temp = m_temp/mean2(m_temp);
    end    

	% set the average level to half of the dynamic range
	m_temp = m_temp * 2^(num_bits-1);
	m_temp = m_temp .* (m_temp >=0);
	m_temp(m_temp>(2^num_bits - 1)) = (2^num_bits - 1);
	m_temp = uint16(m_temp);

    
    plot(m_avg_mask, 'o-')
    hold on
    stem(i, 1, 'r')
    hold off
% 	imagesc(m_temp)
% 	axis image;
% 	colormap gray
	title([num2str(i)])
	pause(0.1)

	% save the image
	m_saveDir = [writeFolder num2str(i) '.tif'];
	imwrite(m_temp, m_saveDir, 'tif');
end