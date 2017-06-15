function [img_z, H] = Propagate( img, pixelsize, refidx, lamda, z, convunits, tf_zeropad, tf_freqmask )
%Written by Yibo Zhang, zybmax09@gmail.com, Oct 2013
%BACKPROPFUNC propagates img with pixelsize, refidx, wavelength by a
%distance of z
%ConvUnits=1 means using nm for wavelength, um for pixelsize and z
%Default: convunits=true, tf_zeropad=true, tf_freqmask=true

if nargin <= 7,
    tf_freqmask = true;
end
if nargin <= 6,
    tf_zeropad = true;
end
if nargin <= 5,
    convunits = true;
end

if convunits,
    lamda = lamda*1e-9;
    pixelsize = pixelsize*1e-6;
    z = z*1e-6;
end

Nh = size(img, 2); %horizontal sample points
Nv = size(img, 1); %vertical sample points
meanvalue = mean2(img);

if tf_zeropad,
    img = padarray(img, [Nv, Nh], meanvalue, 'post');
end

spectrum = fft2(img);
spectrum = fftshift(spectrum);

[NFv, NFh] = size(spectrum);
Fs = 1/pixelsize;
Fh = Fs/NFh .* (-ceil((NFh-1)/2) : floor((NFh-1)/2));
Fv = Fs/NFv .* (-ceil((NFv-1)/2) : floor((NFv-1)/2)); 
[Fhh, Fvv] = meshgrid(Fh, Fv);

H = PropGeneral(Fhh, Fvv, lamda, refidx, z);
freqmask = BandLimitTransferFunction(pixelsize, z, lamda, Fvv, Fhh);

spectrum_z = spectrum.*H;
if tf_freqmask,
    spectrum_z = spectrum_z.*freqmask;
end

spectrum_z = ifftshift(spectrum_z);

img_z = ifft2(spectrum_z);
img_z = img_z(1:Nv, 1:Nh);

end


function H = PropGeneral(Fhh, Fvv, lamda, refidx, z)

DiffLimMat = ones(size(Fhh));
lamdaeff = lamda/refidx;
DiffLimMat((Fhh.^2+Fvv.^2) >= 1/lamdaeff^2) = 0;
H = exp(1j.*2.*pi.*z./lamdaeff.*(1-(lamdaeff.*Fvv).^2-(lamdaeff*Fhh).^2).^0.5);
H(~DiffLimMat) = 0;

end

function [freqmask] = BandLimitTransferFunction(pixelsize, z, lamda, Fvv, Fhh)

[hSize, vSize] = size(Fvv);
dU = (hSize*pixelsize)^-1;
dV = (vSize*pixelsize)^-1;
Ulimit = ((2*dU*z)^2+1)^-0.5/lamda;
Vlimit = ((2*dV*z)^2+1)^-0.5/lamda;

freqmask = (Fvv.^2./(Ulimit^2)+Fhh.^2.*(lamda^2))<=1 & (Fvv.^2.*(lamda^2)+Fhh.^2./(Vlimit^2))<=1;

end