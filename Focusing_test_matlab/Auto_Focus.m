% function z = Auto_Focus_v5(I,l,u,delta_z,pixel_size,refidx,lambda,tf_phase,ROI,tf_rf,tf_plot)
% Auto Focus using Tamura Coefficient to find z2 of lens-free imaging setup
%
% 05/15/2015, created by Chris, chris.yichen.wu@gmail.com
% 07/22/2015, modified by Chris, changed the way to define propagation 
%   to address complaints.:(... Also added phase focus capability :)
% 11/14/2015, modified by Chris, changes include:
%   - fixed a search bug, thanks to Max
%   - use a handle to access to rough focus indicator figure
%   - in rough focus, find the most significant peak instead of finding maximum
% 03/27/2016, v5.1 by Chris.
%   - input ROI can be a scalar, which defines the r_ROI from the center
% 07/30/2016, v6 by Chris.
%   - *speed improvement* by ~30%-60%, including:
%   - FFT first, then each propagate only requires one IFFT
%   - fftshift the angular spectrum, thus saving the time to fftshift
%   - persistent variable of the grid generated
%   - Note: input images of odd number of pixels may have negligible (1e-16) numerical artifacts
% 08/10/2016, v6.1 by Chris.
%   - bug fix, propagate by z2 instead of -z2 to avoid ambiguity
% 08/13/2016, v6.2 by Chris & Yicheng Li
%   - bug fix, the rough focus propagated -z2_list
%
% Input Variables:
%   I           - Ny-by-Nx double, input hologram for auto-focus
%   [l,u]       - scalar, search range of z2, in m
%   delta_z     - scalar, tolerance of focus distance, in m
%   pixel_size  - scalar, pixel_size, in m
%   refidx      - scalar, refractive index, in 1
%   lambda      - scalar, wavelength, in m
%   tf_phase    - boolean, 1 for focusing phase, 0 for amplitude, default: 0 (amplitude)
%   ROI         - 1-by-4 double, ROI for focus 1x4 matrix [lx,ly,ux,uy] in pixels, (NOTE we don't propagate a small region,
%               instead, we propagate whole 'I', then only calculate the focus criteria of this ROI. If you want
%               to propagate small region, you may alternatively crop the 'I' first before feed into this code. DEFAULT is whole ROI)   
%   tf_rf       - boolean, rough search to check concavity or not, default: 1 (do rough focus)
%   tf_plot     - boolean, plot the rough search curve or not, default: 1 (plot the rough search curve)
%
% Output Variables:
%   z           - scalar, auto-focused propagation distance, in m
%   Tamuras     - (optional) 2-by-4M double, the first auto-focus curve
%

function [z,Tamuras] = Auto_Focus(I,l,u,delta_z,pixel_size,refidx,lambda,tf_phase,ROI,tf_rf,tf_plot,tf_verbose)

% clear the persistent variables in the propagation function
clear Auto_Focus

% deal with default cases
[Ny,Nx] = size(I);
% cut one pixel at the end if the pixel count is odd
if mod(Ny,2),   Ny = Ny-1;  I = I(1:Ny,:);     end
if mod(Nx,2),   Nx = Nx-1;  I = I(:,1:Nx);     end

if nargin <= 7 || isempty(tf_phase),    tf_phase = 0;                           end
if nargin <= 8 || isempty(ROI)
%     [Ny,Nx] = size(I)  
    ROI = [1,1,Nx,Ny];  
elseif isscalar(ROI)    % this means we are giving a radius of the ROI
%     [Ny,Nx] = size(I);
    yc = floor(Ny/2)+1;
    xc = floor(Nx/2)+1;
    ROI = [xc-ROI,yc-ROI,xc+ROI-1,yc+ROI-1];
end
if nargin <= 9 || isempty(tf_rf),       tf_rf = 1;                              end
if nargin <= 10 || isempty(tf_plot),    tf_plot = 1;                            end
if nargin <= 11 || isempty(tf_verbose), tf_verbose = true;                      end
% ROI = ROI-zeros(1,4);

% rough focus to make sure of concavity
if tf_rf
    [l,u,Tamuras] = rough_focus(I,l,u,ROI,tf_phase,pixel_size,refidx,lambda,tf_plot,tf_verbose);
end

if isnan(l)||isnan(u)
    z = NaN;    % there's no peak in the rough searched region
else
    % fine search using golden ratio in concave region
    alpha = (sqrt(5)-1)/2;
    p = u-alpha*(u-l);
    fp = focus_error(I,p,ROI,tf_phase,pixel_size,refidx,lambda);
    q = l+alpha*(u-l);
    fq = focus_error(I,q,ROI,tf_phase,pixel_size,refidx,lambda);
    while u-l > delta_z
        if fp < fq
            l = p;
            p = q;
            fp = fq;
            q = l+alpha*(u-l);
            fq = focus_error(I,q,ROI,tf_phase,pixel_size,refidx,lambda);
        else
            u = q;
            q = p;
            fq = fp;
            p = u-alpha*(u-l);
            fp = focus_error(I,p,ROI,tf_phase,pixel_size,refidx,lambda);
        end 
    end
    z = (l+u)/2;
end

end

function [l,u,Tamuras] = rough_focus(I,l0,u0,ROI,tf_phase,pixel_size,refidx,lambda,tf_plot,tf_verbose)
% recursively shrink [zl,zu] to make search region convex

% flexible parameters, increasing this may increase your computational time
M = 3;    % totally 4M+1 rough search points

% kernal starts here
N = 2*M;  % totally 2N+1 search Points for concavity

z_list = linspace(l0,u0,2*N+1);
g_list = zeros(1,2*N+1);
for ii = 1:2*N+1
    g_list(ii) = focus_error(I,z_list(ii),ROI,tf_phase,pixel_size,refidx,lambda);
end
Tamuras = [z_list;g_list];

if tf_plot, h = figure(999); plot(z_list,g_list); drawnow;      end

flag_concave = 0;
while ~flag_concave
    if tf_plot, figure(h); plot(z_list,g_list); drawnow;   end
    flag_concave = 1;
    % check concavity
    for ii = 2:2*N
%         if (g_list(ii-1)+g_list(ii+1)) > (2*g_list(ii))   % concave criteria
        if g_list(ii) <= min([g_list(ii-1),g_list(ii+1)])   % quasi-concave criteria
            flag_concave = 0;
            if tf_verbose    
%                 disp('z range is not concave,half-shrinking...'); 
            end
            break;
        end
    end
    [~,lsor] = findpeaks(g_list,'SortStr','descend');
    if isempty(lsor)
        if tf_verbose,        warning('No Peak Found: Adjust your search range'); end
        z_list(:) = NaN;
        break;
    end
    max_ind = lsor(1);
    if max_ind == 1
        if tf_verbose,        warning('Lowerbound Too Large: Decrease Lowerbound'); end
        z_list(:) = NaN;
        break;
    end
    if max_ind == 2*N+1
        if tf_verbose,         warning('Upperbound Too Small: Increase Upperbound'); end
        z_list(:) = NaN;
        break;
    end
    if max_ind-M < 1
        start_ind = 1;
    elseif max_ind+M>2*N+1
        start_ind = 2*N+1-2*M;
    else
        start_ind = max_ind-M;
    end
    temp_z = z_list;
    temp_g = g_list;
    % keep N+1 points unchanged to simplify calculation complexity
    for ii = 1:N+1
        z_list(2*ii-1) = temp_z(start_ind+(ii-1));
        g_list(2*ii-1) = temp_g(start_ind+(ii-1));
    end
    % calculate focus error of additional N points, total 2N+1 points
    for ii = 1:N
        z_list(2*ii) = (z_list(2*ii-1)+z_list(2*ii+1))/2;
        g_list(2*ii) = focus_error(I,z_list(2*ii),ROI,tf_phase,pixel_size,refidx,lambda);
    end

end
l = z_list(1);
u = z_list(2*N+1);

end

function err = focus_error(I,z,ROI,tf_phase,pixel_size,refidx,lambda)

% define the propagation function
PROP = @(x,y) Propagate_multi(x,pixel_size,refidx,lambda,y);

U = PROP(I,z);
U = U(ROI(2):ROI(4),ROI(1):ROI(3));
if tf_phase == 0
    err = Tamura_coeff(abs(U));
elseif tf_phase == 1
    phase_hold = angle(U)+pi;
    err = Tamura_coeff(abs(phase_hold));
else
    err_amp = Tamura_coeff(abs(U));
    phase_hold = angle(U)+pi;
    err_phase = Tamura_coeff(abs(phase_hold));
    err = (err_amp.^(1-tf_phase))*(err_phase.^tf_phase);
end

end

function C = Tamura_coeff(I)    % Tamura coefficient criteria

% N_crop = 10;     % crop the edge to avoid edge effect
% [Ny,Nx] = size(I);
% I = I(N_crop:Ny-N_crop,N_crop:Nx-N_crop);

std_I = std(I(:));
mean_I = mean(I(:));
C = log(sqrt(std_I/mean_I)); %%%%% negative is minimun positive is max

end

function Uz = Propagate_multi(U,pixelsize,refidx,lambda,z)
% optimized for propagation of the same image multiple times
% - use persistent variables to eliminate repeated computation

persistent U_hat
persistent Fhh Fvv

FT = @(x) fft2(x);
IFT = @(x) ifft2(x);

if isempty(U_hat)
    U_hat = FT(U);
    [NFv, NFh] = size(U_hat);
    Fs = 1/pixelsize;
    Fh = Fs/NFh .* (-ceil((NFh-1)/2) : floor((NFh-1)/2));
    Fv = Fs/NFv .* (-ceil((NFv-1)/2) : floor((NFv-1)/2)); 
    [Fhh, Fvv] = meshgrid(Fh, Fv);
end

lamdaeff = lambda/refidx;
H = exp(1j.*2.*pi.*z./lamdaeff.*(sqrt(1-(lamdaeff.*Fvv).^2-(lamdaeff*Fhh).^2)-1)); % added -1 here to avoid DC term
H((Fhh.^2+Fvv.^2) >= 1/lamdaeff^2) = 0;
H = ifftshift(H);

Uz = IFT(U_hat.*H);

end
