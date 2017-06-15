%% Created by Wei LUO @ EE, UCLA
% Date (MM-DD-YYYY): 10-10-2015
% Version: 1.1
%%

% this code combines two reconstructed images together
% the images are created by script test_subtractionPhaseRetrieval.m

close all
clc
clear

% specify the folders for combination
folder_Ang1 = '/Users/weiluo/Documents/Data/151202_300FPS/ROI7-1_tif/ROI7-1_avg/PR_reconst_525_Ang1_nm/';
folder_Ang2 = '/Users/weiluo/Documents/Data/151202_300FPS/ROI7-1_tif/ROI7-1_avg/PR_reconst_530_Ang2_nm/';
Folder_Write = '/Users/weiluo/Documents/Data/151202_300FPS/ROI7-1_tif/ROI7-1_avg/PR_reconst_combined/';

% folder_Ang1 = 'H:\Data\151005_300FPS\200X_sample2\5-8\5-8-tif5-8_BakSubtracted\PR_reconst_527.1nm\';
% folder_Ang2 = 'H:\Data\151005_300FPS\200X_sample2\5-8\5-8-tif5-8_BakSubtracted\PR_reconst_530_Ang2_nm\';
% Folder_Write = 'H:\Data\151005_300FPS\200X_sample2\5-8\5-8-tif5-8_BakSubtracted\Ang2_nm_comp\';

crop_row_Ang1 = 1:800;
crop_col_Ang1 = 2631:3320; % 1061:1750; 

crop_row_Ang2 = 1:800;
crop_col_Ang2 = 2631:3320; % 1061:1750; 


% start combining
file_type = [folder_Ang1 '*.tif'];
files = dir(file_type);
cc = cellstr([{files(:).name}]);
[cs,index] = sort_nat(cc,'ascend');
readin_files = files(index);   % now the files are sorted by name

num_of_files = length(readin_files);	%get number of the files

split_col_idx = 400 * 6; % index of the column that is used for image combination

%% start combining the images
if(exist(Folder_Write, 'dir')==0)
	mkdir(Folder_Write)
end

for i = 1:num_of_files

	fprintf('Combining %d of %d\n', i, num_of_files);

	% get image 1
	m_readDir1 = [folder_Ang1 files(i).name];
	m_img1 = imread(m_readDir1);
	
	% m_img_combine = zeros(size([m_img1(crop_row_Ang1, crop_col_Ang1), m_img1(crop_row_Ang2, crop_col_Ang2)]));

	m_readDir2 = [folder_Ang2  strrep(files(i).name, '.tiff','.tif') ];
	if(exist(m_readDir2, 'file') == 0)
		fprintf('File No. %d does not exist in foler 2! skipping...\n', i)
	else
		% get image 2
		m_img2 = imread(m_readDir2);

		% combine the two images
        comb_img1 = m_img1(crop_row_Ang1, crop_col_Ang1);
        comb_img1 = floor(comb_img1/256);
        comb_img2 = m_img2(crop_row_Ang2, crop_col_Ang2);
        
        m_img_combine = [comb_img1 comb_img2];
% 		m_img_combine(1:end, 1:split_col_idx) = m_img1(1:end, 1:split_col_idx);
% 		m_img_combine(1:end, (split_col_idx+1):end) = m_img2(1:end, (split_col_idx+1):end);

        imagesc(m_img_combine)
        axis image; colormap gray
        title(['No. ' num2str(i) ' of ' num2str(num_of_files)])
        pause(0.1)
        
		% save the combined image
        m_img_combine = uint8(m_img_combine);
		m_saveDir = [Folder_Write strrep(files(i).name, '.tif', '_comb.jpeg')];
		imwrite(m_img_combine, m_saveDir, 'jpeg');

	end

end