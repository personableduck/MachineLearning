%% this function is used for compensation of the intensity distribution
function [compensated_img,myimg_filtered] = intensitycompensation(inputimage,bandwidth_x,bandwidth_y, dark_bias)

testimage = inputimage;

%green channels (G1, G2) are used
[num_row, num_column] = size(testimage);

%background level is from experiment
if(nargin>3)
    bacground_level = dark_bias; %64.1440;
else
    bacground_level = 0;
end

myimg = inputimage - bacground_level;

FT_myimg = fftshift(fft2(ifftshift(myimg)));

FWHM_x = bandwidth_x;
FWHM_y = bandwidth_y;

Fx = 1:1:num_column;
center_x = floor(num_column/2) + 1;
c_x = FWHM_x/2/sqrt(2*log(2));

Fy = 1:1:num_row;
center_y = floor(num_row/2) + 1;
c_y = FWHM_y/2/sqrt(2*log(2));

[X, Y] = meshgrid(Fx,Fy);
Gaussian_filter = exp(-((X-center_x)/c_x).^2/2).*...
exp(-((Y-center_y)/c_y).^2/2);

FT_myimg_filtered = Gaussian_filter .* FT_myimg;
myimg_filtered = fftshift(ifft2(ifftshift(FT_myimg_filtered)));
myimg_filtered = abs(myimg_filtered);

compensated_img = myimg./myimg_filtered;

end