% bas_dir     = 'G:\ColorPathology2016_Opt\three_wav_exp\20161201_lung_infection_HE_#1_3wavs\3x3_450_540_600';
% mid_dir     = @(ii) ['angle', num2str(ii+2)];
% las_dir     = 'wavelength1';
% 
wavelength = 633;
BP = @(x,z) Propagate(x, 1.12, 1, wavelength, z, true, false, true);

zs = zeros(8,1);
index_list = [1:7 13];
for ff = 1
    filename = sprintf('%d_1_1.bin', 1+(ff+2-1)*36);
    I = avg_holo;
    I = I(1308:2180, 1596:2756);
    z = af_quick_v5(I,[300 600], 1.12, 1,wavelength,0, 'GiniOfGrad', 1, -1, [], 10, 0.1, 40, true, true );
    If = BP(I, z);
    figure(1); clf; imshow(abs(If),[]); title(ff); drawnow;
    zs(ff) = z;
    disp(z);
end