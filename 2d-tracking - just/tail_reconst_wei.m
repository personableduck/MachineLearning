% tail_reconst
% Created by Wei LUO @ EE, UCLA
% Date (MM-DD-YYYY): 09-02-2014
% Version: 1.0

% this code generates 3D structure of the tail using

function [reconst_vector] = tail_reconst(projection1, projection2, ...
    theta_deg_1, phi_deg_1, theta_deg_2, phi_deg_2, x_range, y_range, z_range)

% variables
% projection1, projection2: a M-2 matrix. First column is x coordinates and
% the second column is y coordinates
% theta_deg and phi_deg are the azumith of the line
% r_tolarance: if the distance between two lines 

    % define the line using parameter approach, t being the parameter
    % x = t * sind(theta_deg) * cosd(phi_deg) - x0
    % y = t * sind(theta_deg) * sind(phi_deg) - y0
    % z = t * cosd(theta_deg) - z0
% first angle
% x1 = x + z * tand(theta_deg_1) * cosd(phi_deg_1);
% y1 = y + z * tand(theta_deg_1) * sind(phi_deg_1);
x1 = projection1(:,1);
y1 = projection1(:,2);

% second angle
% x2 = x + z * tand(theta_deg_2) * cosd(phi_deg_2);
% y2 = y + z * tand(theta_deg_2) * sind(phi_deg_2);
x2 = projection2(:,1);
y2 = projection2(:,2);

% figure
% plot(x1, y1, 'b.-', 'linewidth', 2);
% hold on
% plot(x2, y2, 'r.-', 'linewidth', 2);
% legend({'3D Tail', 'Angle 1', 'Angle 2'})
% axis image
% set(gca, 'fontsize', 20)
%% now start generate the 2D projections

num_row = length(x_range);
num_col = length(y_range);

dx = x_range(2) - x_range(1);
dy = y_range(2) - y_range(1);


% frist angle
x1_idx = ceil((x1 - min(x_range))/dx);
y1_idx = ceil((y1 - min(y_range))/dy);

proj_1 = zeros(num_row, num_col);

for i = 1:length(x1_idx)
    proj_1(y1_idx(i), x1_idx(i)) = 1;
end

% second angle
x2_idx = ceil((x2 - min(x_range))/dx);
y2_idx = ceil((y2 - min(y_range))/dy);

proj_2 = zeros(num_row, num_col);

for i = 1:length(x1_idx)
    proj_2(y2_idx(i), x2_idx(i)) = 1;
end

% figure
% subplot(1,2,1)
% imagesc(x_range, y_range, proj_1); axis image; colormap gray
% title('Angle 1')
% subplot(1,2,2)
% imagesc(x_range, y_range, proj_2); axis image; colormap gray
% title('Angle 2')

% create the mask
scan_z = z_range;


reconst_x = zeros(1, 10000);
reconst_y = zeros(1, 10000);
reconst_z = zeros(1, 10000);

[X_RANGE, Y_RANGE] = meshgrid(x_range, y_range);

add_idx = 1;
imge_count = 1;
% figure
% pause
for i = 1:length(scan_z)
    % shift the first projection
    temp1 = proj_1;
    temp_shift_row = -round(scan_z(i) * tand(theta_deg_1) * sind(phi_deg_1)/dy);
    temp1 = circshift(temp1, [temp_shift_row, 0]);
    if(temp_shift_row>0)
        temp1(1:temp_shift_row, :) = temp1(1:temp_shift_row, :) * 0;
    else
        temp1((end-temp_shift_row):end, :) = temp1((end-temp_shift_row):end, :) * 0;
    end
    
    temp_shift_col = -round(scan_z(i) * tand(theta_deg_1) * cosd(phi_deg_1)/dx);
    temp1 = circshift(temp1, [0, temp_shift_col]);
    if(temp_shift_col>0)
        temp1(:, 1:temp_shift_col) = temp1(:, 1:temp_shift_col) * 0;
    else
        temp1(:, (end-temp_shift_col):end) = temp1(:, (end-temp_shift_col):end) * 0;
    end
    
    %shift the second projection
    temp2 = proj_2;
    temp_shift_row = -round(scan_z(i) * tand(theta_deg_2) * sind(phi_deg_2)/dy);
    temp2 = circshift(temp2, [temp_shift_row, 0]);
    if(temp_shift_row>0)
        temp2(1:temp_shift_row, :) = temp2(1:temp_shift_row, :) * 0;
    else
        temp2((end-temp_shift_row):end, :) = temp2((end-temp_shift_row):end, :) * 0;
    end
    
    temp_shift_col = -round(scan_z(i) * tand(theta_deg_2) * cosd(phi_deg_2)/dx);
    temp2 = circshift(temp2, [0, temp_shift_col]);
    if(temp_shift_col>0)
        size(temp2)
        temp_shift_col
        temp2(:, 1:temp_shift_col) = temp2(:, 1:temp_shift_col) * 0;
    else
        temp2(:, (end-temp_shift_col):end) = temp2(:, (end-temp_shift_col):end) * 0;
    end    
    
    % get the total mask
    totalPmask = temp1 .* temp2;
    temp_idx = find(totalPmask>0);
    
    string_length = length(temp_idx);
    temp_X = X_RANGE(temp_idx);
    temp_Y = Y_RANGE(temp_idx);
    temp_Z = ones(1, length(temp_Y)) * scan_z(i);
    
    if(string_length>0)
        reconst_x(add_idx: (add_idx + string_length - 1)) = temp_X.';
        reconst_y(add_idx: (add_idx + string_length - 1)) = temp_Y.';
        reconst_z(add_idx: (add_idx + string_length - 1)) = temp_Z.';    
        add_idx = add_idx + string_length;
    end
    
%     subplot(1,3,1)
%     [ temp1_show ] = SmoothPatternBoundary(temp1, 4 );
%     imagesc(x_range, y_range, temp1_show); axis image; colormap gray; title('Angle 1')
%     axis([-10 50 -45 5])
%     subplot(1,3,2)
%     [ temp2_show ] = SmoothPatternBoundary(temp2, 4 );
%     imagesc(x_range, y_range, temp2_show); axis image; colormap gray; title('Angle 1')
%     axis([-10 50 -45 5])
%     subplot(1,3,3)
%     [ totalPmask_show ] = SmoothPatternBoundary(totalPmask, 6 );
%     totalPmask_show = totalPmask_show/(max(totalPmask_show(:))+1e-6);
%     imagesc(x_range, y_range, totalPmask_show); axis image; colormap gray; title('Combine')
%     axis([-10 50 -45 5])
%     suptitle(['z = ' num2str(scan_z(i)) ' um'])
%     caxis([0 1.85])
%     colormap jet
%     pause(1)
    
    % savefilename = [savefolder num2str(imge_count) '.jpg'];
    % saveas(gca, savefilename)
    imge_count = imge_count + 1;
end

    reconst_x = reconst_x(1:(add_idx-1));
    reconst_y = reconst_y(1:(add_idx-1));
    reconst_z = reconst_z(1:(add_idx-1));    
    
%     figure
%     plot3(reconst_x, reconst_y, reconst_z, 'o')
%     grid on
%     hold on
%     plot3(x,y,z, 'r.-', 'linewidth', 3)
%     set(gca, 'fontsize', 20)    
        
    reconst_vector = [reconst_x(:), reconst_y(:), reconst_z(:)];
end