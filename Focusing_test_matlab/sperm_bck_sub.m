function [out] = sperm_bck_sub(raw,holo_avg)
%% Sperm Stack Background Subtract

avg = holo_avg;

out = raw-avg+mean(avg(:));
end