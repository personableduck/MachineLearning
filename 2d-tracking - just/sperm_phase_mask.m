%% Created by Wei LUO @ EE, UCLA
% Date (MM-DD-YYYY): 12-14-2015
% Version: 1.0

% this function generates the thresholded mask
function [mask_out, new_map] = sperm_phase_mask(field_in, level)

% mask_out: mask defined by the function
% new_map: map that is used to generate the mask

    if(nargin<2)
        level = 0.3; % default value for sperm head
        % Bovine: 0.3
    end

%     figure; imagesc(angle(field_in)); axis image; colormap gray
    
    a = field_in - mean2(field_in);
    b = abs(a);
    b = b./max(b(:));
    % figure; imagesc(b); axis image; colormap gray
    % [level, em] = graythresh(b); 
    mask_out = double(b>level);
    new_map = b;
end