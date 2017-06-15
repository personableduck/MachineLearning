


format_spec = '%f\n';
fid = fopen('./af_output.txt', 'r');
focus_curve = fscanf(fid, format_spec);
fclose(fid);
figure;
plot(focus_curve);