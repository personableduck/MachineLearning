% function z = Auto_Focus_v4(I,l,u,delta_z,pixel_size,refidx,lambda,focusChannel,ROI,focusMeas)

% Auto Focus using Tamura Coefficient

% 05/15/2015, created by Chris, chris.yichen.wu@gmail.com

% 07/22/2015, modified by Chris, changed the way to define propagation 

% to address complaints.:(... Also added phase focus capability :)

% Modified by Max zybmax09@gmail.com

% 

% Input Variables:

%   I       - input hologram for auto-focus

%   [l,u]   - search range of z2, in m

%   delta_z - tolerance of focus distance

%   pixel_size - pixel_size, in m

%   refidx  - refractive index, in 1

%   lambda  - wavelength, in m

%   focusChannel  - 1 for focusing phase, 0 for amplitude. 2 for

%           Background subtraction

%   focusMeas  - 'Tamura', 'Sobel', 'Binariness', 'Concentration', 'Gini',

%             'Energy', 'Entropy', 'AutoCorr', 'Range'; Default is 'Tamura'

%   rect  - [lh, lv, wh, wv];

%

% Output Variables:

%   z       - auto-focused propagation distance 

%



function z = Auto_Focus_v5(I,l,u,delta_z,pixel_size,refidx,lambda,focusChannel,focusMeas,focusDrv,polarity,rect,tf_verb,tf_fig)



if isempty(focusChannel), 

    focusChannel = 0;  %Default focus channel is amplitude

end

    

if isempty(focusMeas),

    focusMeas = 'Tamura';

end



if isempty(focusDrv),

    focusDrv = 1;

end



if isempty(polarity), 

    polarity = 1;

end



if isempty(rect),

    [Ny,Nx] = size(I);

    ROI = [1,1,Nx,Ny];

else

    rect = round(rect); ROI = [rect(1), rect(2), rect(1)+rect(3), rect(2)+rect(4)];

end



if isempty(tf_verb), 

    tf_verb = true;

end



if isempty(tf_fig), 

    tf_fig = true;

end



if strcmpi(focusMeas, 'Tamura'),

    foc = @(a) Tamura_coeff(a);

elseif strcmpi(focusMeas, 'Sobel'),

    foc = @(a) gradient_Sobel(a);

elseif strcmpi(focusMeas, 'Binariness'),

    foc = @(a) Binariness(a);

elseif strcmpi(focusMeas, 'Concentration'),

    foc = @(a) Concentration(a);

elseif strcmpi(focusMeas, 'Gini'),

    foc = @(a) Gini(a);

elseif strcmpi(focusMeas, 'Energy'),

    foc = @(a) Energy(a);

elseif strcmpi(focusMeas, 'Entropy'),

    foc = @(a) Entropy(a);

elseif strcmpi(focusMeas, 'AutoCorr'),

    foc = @(a) AutoCorr(double(a));

elseif strcmpi(focusMeas, 'Range'),

    foc = @(a) Range(a);

else

    foc = @(a) Tamura_coeff(a);

end



% rough focus to make sure of concavity

[l,u] = rough_focus(I,l,u,ROI,focusChannel,pixel_size,refidx,lambda,foc,focusDrv,polarity,tf_verb,tf_fig);



if isnan(l)||isnan(u)

    z = NaN;

else

    % fine search using golden ratio in concave region

    alpha = (sqrt(5)-1)/2;

    p = u-alpha*(u-l);

    fp = focus_error(I,p,ROI,focusChannel,pixel_size,refidx,lambda,foc,focusDrv,polarity);

    q = l+alpha*(u-l);

    fq = focus_error(I,q,ROI,focusChannel,pixel_size,refidx,lambda,foc,focusDrv,polarity);



    while u-l > delta_z

        if fp < fq

            l = p;

            p = q;

            fp = fq;

            q = l+alpha*(u-l);

            fq = focus_error(I,q,ROI,focusChannel,pixel_size,refidx,lambda,foc,focusDrv,polarity);

        else

            u = q;

            q = p;

            fq = fp;

            p = u-alpha*(u-l);

            fp = focus_error(I,p,ROI,focusChannel,pixel_size,refidx,lambda,foc,focusDrv,polarity);

        end 

            disp(l)

            disp(u);

    end



    z = (l+u)/2;

end



% toc





end



% recursively shrink [zl,zu] to make search region convex

function [l,u] = rough_focus(I,l0,u0,ROI,focusChannel,pixel_size,refidx,lambda,foc,focusDrv,polarity,tf_verb,tf_fig)



% flexible parameters

M = 3;    % totally 4M+1 rough search points



if 4*M + 1 < round((u0-l0)/1e-5), 

    M = round(round((u0-l0)/1e-5) / 4); 

end



% kernal starts here

N = 2*M;  % totally 2N+1 search Points for concavity



z_list = linspace(l0,u0,2*N+1);

g_list = zeros(1,2*N+1);

for ii = 1:2*N+1

    g_list(ii) = focus_error(I,z_list(ii),ROI,focusChannel,pixel_size,refidx,lambda,foc,focusDrv,polarity);

end

if tf_fig, figure(20); plot(z_list,g_list); drawnow; end

flag_concave = 0;

while ~flag_concave

    flag_concave = 1;

    % check concavity

    for ii = 2:N*2

%         if (g_list(ii-1)+g_list(ii+1)) > (2*g_list(ii))   % concave criteria

        if g_list(ii) <= min([g_list(ii-1),g_list(ii+1)])   % quasi-concave criteria

            flag_concave = 0;

            if tf_verb, disp('z range is not concave,half-shrinking...'); end

            break;

        end

    end

    [~,max_ind] = max(g_list);

    if max_ind == 1

        warning('Lowerbound Too Large: Decrease Lowerbound');

        z_list(:) = NaN;

        break;

    end

    if max_ind == 2*N+1

        warning('Upperbound Too Small: Increase Upperbound');

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

        g_list(2*ii) = focus_error(I,z_list(2*ii),ROI,focusChannel,pixel_size,refidx,lambda,foc,focusDrv,polarity);

    end

if tf_fig, figure(20); plot(z_list,g_list); drawnow; end

%     disp(diff(z_list));

%     disp(g_list);

end

l = z_list(1);

u = z_list(2*N+1);



end



function err = focus_error(I,z,ROI,focusChannel,pixel_size,refidx,lambda,foc,focusDrv,polarity)



dz = 0.1e-6; %the delta_z for creating the derivatives

if focusDrv == 0, 

    zs = z; 

elseif focusDrv == 1, 

    zs = [z+dz, z-dz];

elseif focusDrv == 2, 

    zs = [z+dz, z, z-dz];

end

errs = zeros(size(zs));

% define the propagation function

PROP = @(x,y) Propagate_ver3(x,pixel_size,refidx,lambda,y,false);

Us = zeros([size(I), numel(zs)]);

for ii = 1: numel(zs), 

    Us(:,:,ii) = PROP(I,-zs(ii));

end

Us = Us(ROI(2):ROI(4),ROI(1):ROI(3),:);

% find focus measure for each image

for ii = 1: numel(zs), 

    U = Us(:,:,ii);

    if focusChannel == 0

        errs(ii) = foc(abs(U));

    elseif focusChannel == 1

        phase_hold = angle(U*exp(1i*2*pi*refidx*zs(ii)/lambda));

        for kk = 1: 5,

            phase_hold = mod(phase_hold-mean2(phase_hold)+pi, 2*pi) - pi;

        end

        errs(ii) = foc(abs(phase_hold + pi));

    elseif focusChannel == 2

        errs(ii) = foc(abs(U-mean(U(:))));

    else

        err_amp = foc(abs(U));

        phase_hold = angle(U*exp(1i*2*pi*refidx*zs(ii)/lambda));

        err_phase = foc(abs(phase_hold));

        errs(ii) = (err_amp.^(1-focusChannel))*(err_phase.^focusChannel);

    end

end

% find derivative

if focusDrv == 0, 

    err = errs(1); 

elseif focusDrv == 1, 

    err = (errs(1)-errs(2))/2/dz;

elseif focusDrv == 2, 

    err = (errs(1)+errs(3)-2*errs(2))/dz^2;

end

err = err*polarity;

end





function f = Tamura_coeff(I)



N_crop = 10;     % crop the edge to avoid edge effect

[Ny,Nx] = size(I);

I = I(N_crop:Ny-N_crop,N_crop:Nx-N_crop);



std_I = std(I(:));

mean_I = mean(I(:));

f = sqrt(std_I/mean_I);





end



function f = gradient_Sobel(I)



N_crop = 10;     % crop the edge to avoid edge effect

[Ny,Nx] = size(I);

I = abs(I);



Sobel_y = -fspecial('sobel');

Sobel_x = Sobel_y';



gx = filter2(Sobel_x,I);

gy = filter2(Sobel_y,I);



gx_crop = gx(N_crop:Ny-N_crop,N_crop:Nx-N_crop);

gy_crop = gy(N_crop:Ny-N_crop,N_crop:Nx-N_crop);



f = sum(sqrt(gx_crop(:).^2+gy_crop(:).^2));



end



function f = Binariness(I)



I = I(:) - min(I(:));

I = I/max(I(:));

binariness = I.^2 .* (I-1).^2;

f = sum(binariness(:));



end



function f = Concentration(I)



p = I.^4;  p = sum(p(:));

q = I.^2;  q = sum(q(:));

f = p/(q.^2);





end



function f = Gini(I)



I = sort(I(:));

c = sum(abs(I(:)))*ones(length(I),1);

N = length(I);

g = (I./c).* (N + 0.5 - (1:N)')./N;

f = 1 - 2*(sum(g));



end



function f = Energy(I)

    f = sum(abs(I(:)));

end



function f = Entropy(I)



    binrange = linspace(min(I(:)),max(I(:)),1000);

    N = histc(I(:),binrange);

    N(N==0) = 0.0000001;

    N = N/max(N);

    f = -1*sum(N.*log2(N));

    

end



function f = AutoCorr(I)



    [v,h] = size(I);

    m = mean(I(:));

    x1 = [I(:,2:end), m*ones(v,1)];

    x2 = [I(:,3:end), m*ones(v,2)];

    y1 = [I(2:end,:); m*ones(1,h)];

    y2 = [I(3:end,:); m*ones(2,h)];

    

    corr1 = I.*x1 + I.*y1;

    corr2 = I.*x2 + I.*y2;

    f = sum(corr1(:))-sum(corr2(:));



end



function f = Range(I)



    binrange = linspace(min(I(:)),max(I(:)),1000);

    N = histc(I(:),binrange);

    f = max(N(N>0))-min(N(N>0));

    

end





% function Uz = Propagate_ver3(U,pixelsize,refidx,lambda,z,tf_freqmask)

% Propagation in frequency by division of angular spectrum

% Note: (0,0) assumed to be centered naturally

% Inputs:   spectrum    -   input image in frequency domain

%           pixelsize   -   pixel size of sensor, in m

%           refidx      -   refractive index

%           lambda      -   wavelength, in m

%           z           -   propagation distance

%           tf_freqmask -   0 or 1, add frequency mask to avoid ringing

% 

% 11/20/2014, Created by (Chris) Yichen WU, chris.yichen.wu@gmail.com

% 01/05/2015, ver2 by (Chris) Yichen WU, chris.yichen.wu@gmail.com

%       -   deleted so-called freq_mask

%       -   combined PropGeneral to one function

% 01/15/2015, ver3 by (Chris) Yichen WU, chris.yichen.wu@gmail.com

%       -   added freq_mask



function Uz = Propagate_ver3(U,pixelsize,refidx,lambda,z,tf_freqmask)



FT = @(x) ifftshift(fft2(fftshift(x)));

IFT = @(x) ifftshift(ifft2(fftshift(x)));



spectrum = FT(U);

[NFv, NFh] = size(spectrum);

Fs = 1/pixelsize;

Fh = Fs/NFh .* (-ceil((NFh-1)/2) : floor((NFh-1)/2));

Fv = Fs/NFv .* (-ceil((NFv-1)/2) : floor((NFv-1)/2)); 

[Fhh, Fvv] = meshgrid(Fh, Fv);



DiffLimMat = ones(NFv,NFh);

lamdaeff = lambda/refidx;

DiffLimMat((Fhh.^2+Fvv.^2) >= 1/lamdaeff^2) = 0;



H = exp(1j.*2.*pi.*z./lamdaeff.*(1-(lamdaeff.*Fvv).^2-(lamdaeff*Fhh).^2).^0.5);

H(~DiffLimMat) = 0;



spectrum_z = spectrum.*H;



% if nargout >1

%     DiffLim = DiffLimMat;

% end



if tf_freqmask

    freqmask = BandLimitTransferFunction(pixelsize, z, lambda, Fvv, Fhh);

    spectrum_z = spectrum_z.*freqmask;

end



Uz = IFT(spectrum_z);



end



function [freqmask] = BandLimitTransferFunction(pixelsize, z, lamda, Fvv, Fhh)



[hSize, vSize] = size(Fvv);

dU = (hSize*pixelsize)^-1;

dV = (vSize*pixelsize)^-1;

Ulimit = ((2*dU*z)^2+1)^-0.5/lamda;

Vlimit = ((2*dV*z)^2+1)^-0.5/lamda;



freqmask = (Fvv.^2./(Ulimit^2)+Fhh.^2.*(lamda^2))<=1 & (Fvv.^2.*(lamda^2)+Fhh.^2./(Vlimit^2))<=1;



end







