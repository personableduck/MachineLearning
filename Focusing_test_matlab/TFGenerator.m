function TF = TFGenerator(n_x_in,n_y_in,dx,dy,z,lambda)
% transfer function generator for SSA-based Propogation

if mod(n_x_in,2)==1     
    fx = 1/dx/n_x_in*(-(n_x_in-1)/2:(n_x_in-1)/2);     % odd
else
    fx = 1/dx/n_x_in*(-(n_x_in)/2:(n_x_in/2-1));       % even
end
if mod(n_y_in,2)==1     
    fy = 1/dy/n_y_in*(-(n_y_in-1)/2:(n_y_in-1)/2);     % odd
else
    fy = 1/dy/n_y_in*(-(n_y_in)/2:(n_y_in/2-1));       % even
end

[fxx fyy] = meshgrid(fx,fy);                  
fzz_2 = (1/lambda)^2-fxx.^2-fyy.^2;     % the square of fzz
circ = fzz_2>0;                         % a filter for the non-propagating components
fzz = sqrt(fzz_2.*circ);                % spatial frequency in z, only the purely real components
clear fxx fyy fzz_2;
exp_z = exp(-1i*2*pi*(fzz.*z));                      % the phase term with fz*z
clear fzz;
TF = exp_z.*circ;                   % transfer function
end