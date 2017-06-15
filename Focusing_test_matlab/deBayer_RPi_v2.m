% function to debayer images on RPi Camera v2
% uses bilinear interpolation
%
% Input Variables:
%   Iin         - Ny-by-Nx double, input image (Bayer pattern)
%   channel     - string, 'R' or 'G' or 'B'
%
% Output Variables:
%   Iout        - Ny-by-Nx double, output (debayered) image
%
%

function Iout = deBayer_RPi_v2(Iin,channel)

% the grids for different channels
[Ny,Nx] = size(Iin);
y = 1:Ny;
x = 1:Nx;
[x_mat,y_mat] = meshgrid(x,y);
xR = 1:2:Nx;
yR = 1:2:Ny;
[xR_mat,yR_mat] = meshgrid(xR,yR);
xB = 2:2:Nx;
yB = 2:2:Ny;
[xB_mat,yB_mat] = meshgrid(xB,yB);
xG1 = 1:2:Nx;
yG1 = 2:2:Ny;
xG2 = 2:2:Nx;
yG2 = 1:2:Ny;

% debayer for given channel = 'R','G','B'
switch(channel)
    case 'R'
        I = Iin(yR,xR);
        Iout = interp2(xR_mat,yR_mat,I,x_mat,y_mat,'linear',0);
    case 'G'
        Iout = Iin;
        Iout(yR,xR) = (Iout(yG1,xG1)+Iout(mod(yG1+1,Ny)+1,xG1)+...
            Iout(yG2,xG2)+Iout(yG2,mod(xG2-3,Nx)+1))/4;
        Iout(yB,xB) = (Iout(yG1,xG1)+Iout(yG1,mod(xG1+1,Nx)+1)+...
            Iout(yG2,xG2)+Iout(mod(yG2-3,Ny)+1,xG2))/4;
    case 'B'
        I = Iin(yB,xB);
        Iout = interp2(xB_mat,yB_mat,I,x_mat,y_mat,'linear',0);
    otherwise
        error('Give input of channel as R,G or B');
end

end