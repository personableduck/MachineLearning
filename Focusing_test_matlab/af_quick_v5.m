
function [ opt_x, opt_val ] = af_quick_v5( A, z_range, ps, RI, wav, ...
    foc_ch, foc_ms, foc_dr, polrty, rect, init_res, delta, quick_zoom, tf_verb, tf_fig )

if isempty(init_res),   init_res = ( z_range(2)-z_range(1) ) / 16;  end
if isempty(delta),      delta = ( z_range(2)-z_range(1) ) * 1e-6;   end
if isempty(quick_zoom), tf_qz = false; else tf_qz = true;           end

if isempty(tf_verb),    tf_verb = true;                             end
if isempty(tf_fig),     tf_fig = true;                              end

if strcmpi(foc_ms, 'Tamura'),
    foc_f = @(a) Tamura(a); % focus function
elseif strcmpi(foc_ms, 'Sobel'),
    foc_f = @(a) gradient_Sobel(a);
elseif strcmpi(foc_ms, 'Binariness'),
    foc_f = @(a) Binariness(a);
elseif strcmpi(foc_ms, 'Concentration'),
    foc_f = @(a) Concentration(a);
elseif strcmpi(foc_ms, 'Gini'),
    foc_f = @(a) Gini(a);
elseif strcmpi(foc_ms, 'Energy'),
    foc_f = @(a) Energy(a);
elseif strcmpi(foc_ms, 'Entropy'),
    foc_f = @(a) Entropy(a);
elseif strcmpi(foc_ms, 'AutoCorr'),
    foc_f = @(a) AutoCorr(double(a));
elseif strcmpi(foc_ms, 'Range'),
    foc_f = @(a) Range(a);
elseif strcmpi(foc_ms, 'GiniOfGrad'),
    foc_f = @(a) GiniOfGrad(a);
else
    foc_f = @(a) Tamura(a);
end

switch foc_ch
    case 0
        foc_t = @(x) abs(x);     % focus transformation
    case 1
        foc_t = @(x) get_phase(x);
    case 2
        foc_t = @(x) abs( x - mean2(x) );
end

foc_v = @(z) polrty .* foc_f(foc_t(im_rect(Propagate(A, ps, RI, wav, z, true, false, false), rect)));

switch foc_dr
    case 0
        eval_func = @(z) foc_v(z);
    case 1
        eval_func = @(z) (foc_v(z+0.1) - foc_v(z-0.1))/0.2;
    case 2
        eval_func = @(z) (foc_v(z+0.1) + foc_v(z-0.1) - 2.*foc(z)) / 0.01;
end

% rough search

l = z_range(1); u = z_range(2); % lower and upper bounds
M = ceil( (u-l) / init_res / 4);% totally 4M+1 points, or 2N+1 points
x_list = [l: (u-l)/4/M: u];      % search points of x
f_list = zeros(size(x_list));   % search points of the function

for ii = 1: 4*M+1,
    f_list(ii) = eval_func( x_list(ii) );
end

if tf_fig, figure(20); clf; plot(x_list,f_list, '.b'); drawnow; title('focus curve'); xlabel('z(\mum)'); ylabel('focus measure'); end % plot result

[opt_val, idx] = min( f_list(:) );
opt_x = x_list(idx);

if idx == 1,
    warning('Lowerbound Too Large: Decrease Lowerbound');
    opt_x = NaN;
    opt_val = NaN;
    return;
end

if idx == numel(x_list),
    warning('Upperbound Too Small: Increase Upperbound');
    opt_x = NaN;
    opt_val = NaN;
    return;
end

% quick zoom
if tf_qz,
    M = 4;
    l = max( opt_x - quick_zoom/2, z_range(1) );
    u = min( opt_x + quick_zoom/2, z_range(2) );
    x_list = [l: (u-l)/4/M: u];       % search points of x
    f_list = zeros(size(x_list));   % search points of the function
    for ii = 1: numel(x_list),
        f_list(ii) = eval_func( x_list(ii) );
    end
end

% recursively shrink [zl,zu] to make search region unimodal (more strict
% than quasi-convex)
N = 2*M;
flag_qcvx = 0;

while ~flag_qcvx
    flag_qcvx = 1;
    
    % check unimodal
    for ii = 2: N*2,
        if f_list(ii) > max( f_list(ii-1),f_list(ii+1) ),    % unimodal criteria
            flag_qcvx = 0;
            if tf_verb, disp('z range is not concave,half-shrinking...'); end
            break;
        end
    end
    
    [~,min_ind] = min(f_list);
    if min_ind == 1
        warning('Lowerbound Too Large: Decrease Lowerbound');
        x_list(:) = NaN;
        break;
    end
    if min_ind == 4*M+1
        warning('Upperbound Too Small: Increase Upperbound');
        x_list(:) = NaN;
        break;
    end
    
    if min_ind-M < 1
        start_ind = 1;
    elseif min_ind+M>2*N+1
        start_ind = 2*N+1-2*M;
    else
        start_ind = min_ind-M;
    end
    
    temp_x = x_list;
    temp_f = f_list;
    
    % keep N+1 points unchanged to simplify calculation complexity
    for ii = 1:N+1
        x_list(2*ii-1) = temp_x(start_ind+(ii-1));
        f_list(2*ii-1) = temp_f(start_ind+(ii-1));
    end
    
    % calculate function value of additional N points, total 2N+1 points
    for ii = 1:N
        x_list(2*ii) = (x_list(2*ii-1)+x_list(2*ii+1))/2;
        f_list(2*ii) = eval_func( x_list(2*ii) );
    end
    
    if tf_fig, figure(20); hold on; plot(x_list,f_list,'.b'); drawnow; end
    
end

l = x_list(1);
u = x_list(2*N+1);

% gold-ratio bisection

if isnan(l)||isnan(u)
    opt_x = NaN;
    
else
    
    % fine search using golden ratio in unimodal region
    alpha = (sqrt(5)-1)/2;
    p = u-alpha*(u-l);
    fp = eval_func(p);
    q = l+alpha*(u-l);
    fq = eval_func(q);
    
    if tf_fig, figure(20); hold on; plot(p, fp, '.b'); hold on; plot(q, fq, '.b'); drawnow; end
    
    while u-l > delta,
        if fp > fq
            l = p;
            p = q;
            fp = fq;
            q = l+alpha*(u-l);
            fq = eval_func(q);
            
            if tf_fig, figure(20); hold on; plot(q, fq, '.b'); drawnow; end
            
        else
            u = q;
            q = p;
            fq = fp;
            p = u-alpha*(u-l);
            fp = eval_func(p);
            
            if tf_fig, figure(20); hold on; plot(p, fp, '.b'); drawnow; end
            
        end
        
        
        
    end
    
    opt_x = (l+u)/2;
end

end

function f = Tamura(I)

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

function f = GiniOfGrad(I)

[gx, gy] = gradient(I);
g = sqrt(abs(gx).^2 + abs(gy).^2 );
f = Gini(g);

end

function B = im_rect(A, rect)

if isempty(rect), B = A;
else rect = round(rect); B = A(rect(2):rect(2)+rect(4), rect(1):rect(1)+rect(3));
end

end

function phi = get_phase(A)
phi = angle(A);
for kk = 1: 5,
    phi = mod(phi-mean2(phi)+pi, 2*pi) - pi;
end
end