%% Created by Wei LUO @ EE, UCLA
% Date (MM-DD-YYYY): 10-09-2015
% Version: 1.1
%%
close all
clc

readFolder = '/Users/Wei/Downloads/SpermImages10FPS/';
write_dir = [readFolder 'reconst2/'];
if(exist(write_dir, 'dir')==0)
    mkdir(write_dir);
end

keyword = '*.tiff';
% list all files
files = dir([readFolder keyword]);
cc = cellstr([{files(:).name}]);
[cs,index] = sort_nat(cc,'ascend');
readin_files = files(index);   % now the files are sorted by name
num_of_files = length(readin_files);	%get number of the files

z = 2400; % unit: um
n_ref = 1.332;
crop_col = (273:757);
crop_row = (139:678);
lambda = 0.645; % wavelength, unit: um
pixel_size = 2.2; % unit:um

interp_factor = 4;


smooth_edge = 30;

% illumination angles
theta_deg = 0;
phi_deg = 0;

% phase retrieval parameters
level = 0.065;


% coversion parameters
show_max = 1; 0.12; 0.1;
show_min = 0; 0.05; 0.03;
num_bit = 16;
convert_type = 'tif';

%% USER SHOULD NOT CHANGE FROM HERE

grid_size = pixel_size /interp_factor;

for img_count = 1:23; 1:num_of_files;
    tic
    m_filename = [readFolder readin_files(img_count).name];
    A = imread(m_filename);
    A = double(A(crop_row, crop_col)); % crop the FOV
    A = A/mean2(A);
    DC = 1;    
    [ A ] = SmoothFOVBoundary(A, smooth_edge, DC);


    % start reconstruction
    grid_size = pixel_size / interp_factor;


    % interpolate the image
    interp_extendflag = 1;
    [ m_image_interp ] = Interpo_img(A, interp_factor, interp_factor, 'spline', interp_extendflag, floor(interp_factor/2) );

   % back propagation
    [field_propagated, fx, fy, x, y] = BackPropagation(m_image_interp, grid_size, grid_size, lambda, n_ref, theta_deg, phi_deg, n_ref, -z);
    field_bak = BackPropagation(m_image_interp, grid_size, grid_size, lambda, n_ref, theta_deg, phi_deg, n_ref, z);

    %% define the mask
    m_mask = (angle(field_bak)>level);
    
    % erode
    SE=strel('disk',3 );
    m_mask=imerode(m_mask, SE);
    
    % dilate
    m_mask=imdilate(m_mask, SE);
    
    % smoothen the mask
    edge_num = 4;
    [ m_mask, PF ] = SmoothPatternBoundary(double(m_mask), edge_num );
    % [ m_mask, PF ] = SmoothPatternBoundary(m_mask, edge_num );

    field_real = real(field_propagated);
    field_imag = imag(field_propagated);

    field_real = field_real.*(1-m_mask) + mean(field_real(:)).* m_mask ;
    field_imag = field_imag.*(1-m_mask) + mean(field_imag(:)).* m_mask ;

    % propagate to mirror plane
    field_propagated = field_real + 1j * field_imag;

    [row_num, col_num] = size(field_propagated);

    [ illumination_field, theta_deg, phi_deg] = ...
        Planewave_field( lambda, n_ref, grid_size, grid_size, row_num, col_num, theta_deg, phi_deg, 'nochange');

    field_propagated = field_propagated .* illumination_field;

    [field_propagated_mirror, fx, fy, x, y] = ...
        Wave_Propagation(field_propagated, grid_size, grid_size, lambda, n_ref, -2*z);

    field_propagated_mirror = field_propagated_mirror .* conj(illumination_field);

    % shift back the image
    shift_v_pix = round(-2*z * tand(theta_deg) * sind(phi_deg) / grid_size);
    shift_h_pix = round(-2*z * tand(theta_deg) * cosd(phi_deg) / grid_size);
    field_propagated_mirror = circshift(field_propagated_mirror, [-shift_v_pix, -shift_h_pix]);
    field_propagated_mirror = field_propagated_mirror / exp(1j*angle(mean2(field_propagated_mirror)));

    
    % convert the image to a displayable dynamic range
    m_show_image = angle(field_propagated_mirror);
    m_show_image = m_mask;
    m_show_image(m_show_image>show_max) = show_max;
    m_show_image(m_show_image<show_min) = show_min;
    m_show_image = m_show_image - show_min;
	m_show_image = m_show_image/(show_max - show_min) * (2^num_bit-1);
	m_show_image = uint16(m_show_image);

    % save the image
	m_save_file_name = strrep(files(img_count).name, '.tiff', ['_PR_reconst.tif']);
	m_write_full_dir = [write_dir m_save_file_name];
	imwrite(m_show_image, m_write_full_dir, 'tif')
    
    fprintf('%d of %d reconstructed in %5.3f sec\n', img_count, num_of_files, toc)
   
    imagesc(m_mask); axis image; colormap gray
end