function u_out = Prop_SSA(u_in, dx, dy, z, lambda)
% Simplified SSA-based Propogation Diffraction Calculator
% Based on spatial spectral analysis, calculate the distribution of a scalar field
% after propagation in free space (or a uniform medium). 
% Constant and perturbation part will be treated separately.
% u_in: input profile, dx,dy: the spacing between grids
% width_out: width of the output plane (square), z: propagation distance, 
% lambda: wavelength (in medium)
%
% arbitrary spatial unit can be used as long as they are consistent between
% dx, dy, z, and lambda.
% 
% x is rightward, y is downward
%

[n_y_in n_x_in] = size(u_in);       % number of grids in source
TF = TFGenerator(n_x_in,n_y_in,dx,dy,z,lambda);   % generate
prod = fft2(u_in).*ifftshift(TF);
clear TF;                       % the product for inverse Fourier transform
u_out = ifft2(prod);            % use inverse Fourier transform to calculate the output 
end