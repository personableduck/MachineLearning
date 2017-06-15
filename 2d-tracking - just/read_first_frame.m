filename=char(file_1(:,1));

[Ref_img, pixelsize] = read_img_wrapper(read_dir,filename,rotation);
    
bandwidth_bak = 16;
    
dark_bias = min(Ref_img(:));
ori_level = mean2(Ref_img);
Ref_img = intensitycompensation(Ref_img,bandwidth_bak,bandwidth_bak, dark_bias);
Ref_img = Ref_img/mean2(Ref_img) * ori_level;

rec = Prop_SSA(Ref_img,pixelsize/interp_factor,pixelsize/interp_factor,...
                    -z2_list(k),lambda_V/1000/n_glass);

%------------------

recon_stack(:,:,1) = rec;

first_frame= abs(recon_stack(:,:,1));

threshold_first = mean2(first_frame)*0.80; %Because low pass filter, we need to control the (*80) figure.
objMask_first = (first_frame<threshold_first); % intensity filter function, get rid of background.

compare_1=recon_stack_diff(:,:,1);
threshold_1 = mean2(compare_1)*6; %<--------------------- need to control dpend on resolution rate
objMask_1 = (compare_1>threshold_1);

object_first_one=objMask_first&objMask_1;
figure, imshow(object_first_one)