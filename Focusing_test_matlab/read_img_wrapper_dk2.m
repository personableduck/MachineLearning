function [Ref_img] = read_img_wrapper_dk(read_dir,filename)
%DK:// changing color image to gray image
%DK:// changed { double -> im2double }
%Ref_img = im2double(rgb2gray(imread([read_dir,filename])));
Ref_img = double(imread([read_dir,filename]));
Ref_img = Ref_img(:,:,2);
end