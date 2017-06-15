imgrr1=double(imread('ICI-11642v2_x3-1_0.0ms.bmp'));
imgrr2=double(imread('ICI-11642v2_x3-1_986.0ms.bmp'));
imgrr3=double(imread('ICI-11642v2_x3-1_496.0ms.bmp'));
imgrr4=double(imread('ICI-11642v2_x3-1_1441.0ms.bmp'));
imgrr5=double(imread('ICI-11642v2_x3-1_1896.0ms.bmp'));
imgrr6=double(imread('ICI-11642v2_x3-1_2352.0ms.bmp'));
imgrr7=double(imread('ICI-11642v2_x3-1_2814.0ms.bmp'));
imgrr8=double(imread('ICI-11642v2_x3-1_3260.0ms.bmp'));
imgrr9=double(imread('ICI-11642v2_x3-1_3720.0ms.bmp'));
imgrr10=double(imread('ICI-11642v2_x3-1_4164.0ms.bmp'));
imgrr11=double(imread('ICI-11642v2_x3-1_4658.0ms.bmp'));
imgrr12=double(imread('ICI-11642v2_x3-1_5434.0ms.bmp'));
imgrr13=double(imread('ICI-11642v2_x3-1_5900.0ms.bmp'));
imgrr14=double(imread('ICI-11642v2_x3-1_6352.0ms.bmp'));
imgrr15=double(imread('ICI-11642v2_x3-1_6827.0ms.bmp'));
imgrr16=double(imread('ICI-11642v2_x3-1_7296.0ms.bmp'));
imgrr17=double(imread('ICI-11642v2_x3-1_7748.0ms.bmp'));
imgrr18=double(imread('ICI-11642v2_x3-1_8217.0ms.bmp'));
imgg=imread('test_11.bmp');
hist=imread('Hist_detection.bmp');

%%%%%%
figure
imshow((abs(recon_stack(:,:,1))+abs(recon_stack(:,:,2)))/2,[])
hold on

dddeee= size(obj_mat_cor);
ymbc=dddeee(1,1);

for pictureN=1:ymbc
hold on
plot((obj_mat_cor{pictureN,1}(:,1)),(obj_mat_cor{pictureN,1}(:,2)),'g')
hold on
plot((obj_mat_cor{pictureN,1}(1,1)),(obj_mat_cor{pictureN,1}(1,2)),'b.')

end

%%%%%%%%

hold on
for jseed=1:233
hold on
plot((x_obj(:,jseed)),(y_obj(:,jseed)),'g')
end
hold on
plot((x_obj(1,:)),(y_obj(1,:)),'b.')

%%%%%%%%%

figure
imshow((abs(recon_stack(:,:,1))+abs(recon_stack(:,:,2)))/2,[])
for jkjk=1:iter-1
hold on
plot((coor_from{1,jkjk}(:,1)),(coor_from{1,jkjk}(:,2)),'g.')
plot((coor_to{1,jkjk}(:,1)),(coor_to{1,jkjk}(:,2)),'g.')
end
hold on
plot((coor_from{1,1}(:,1)),(coor_from{1,1}(:,2)),'b.')

%%%%%%

figure
imshow((abs(recon_stack(:,:,1))+abs(recon_stack(:,:,2)))/2,[])
hold on
for pictureN=1:534
hold on
for
    plot((obj_mat_cor{pictureN,1}(1,1):(obj_mat_cor{pictureN,1}(2,1))

plot((obj_mat_cor{pictureN,1}(:,1)),(obj_mat_cor{pictureN,1}(:,2)),'g')
hold on
plot((obj_mat_cor{pictureN,1}((iter-1)*2,1)),(obj_mat_cor{pictureN,1}((iter-1)*2,2)),'b.','MarkerSize',10))

end

%%%%%%%%

figure
subplot(1,2,1),imshow((abs(imgrr1)+abs(recon_stack(imgrr2)))/2,[])
subplot(1,2,2),imshow((abs(recon_stack(:,:,1))+abs(recon_stack(:,:,2)))/2,[])
hold on
for pictureN=1:dkjj
hold on
subplot(1,2,2),plot((obj_mat_cor{pictureN,1}(:,1)),(obj_mat_cor{pictureN,1}(:,2)),'m')
hold on
subplot(1,2,2),plot((obj_mat_cor{pictureN,1}(1,1)),(obj_mat_cor{pictureN,1}(1,2)),'g.')

end




figure
subplot(1,2,1),imshow((abs(imgrr1)+abs(recon_stack(imgrr2)))/2,[])
subplot(1,2,2),imshow((abs(recon_stack(:,:,1))+abs(recon_stack(:,:,2)))/2,[])
hold on
subplot(1,2,2),plot((coor_from{1,1}(:,1)),(coor_from{1,1}(:,2)),'m.')
hold on
subplot(1,2,2),plot((coor_to{1,1}(:,1)),(coor_to{1,1}(:,2)),'g.')

%%%%%%

figure
plot((coor_from{1,1}(:,1)),(coor_from{1,1}(:,2)),'bo')
hold on
plot((coor_to{1,1}(:,1)),(coor_to{1,1}(:,2)),'g.')
plot((coor_to{1,2}(:,1)),(coor_to{1,2}(:,2)),'g.')
plot((coor_from{1,2}(:,1)),(coor_from{1,2}(:,2)),'g.')
for jkjk=1:19
hold on
plot((coor_from{1,jkjk}(:,1)),(coor_from{1,jkjk}(:,2)),'g.')
plot((coor_to{1,jkjk}(:,1)),(coor_to{1,jkjk}(:,2)),'g.')
end
figure
imshow((abs(recon_stack(:,:,1)))
 imshow((abs(recon_stack(:,:,1)))
                                 |
Error: Expression or statement is incorrect--possibly unbalanced (, {, or [.
 
Did you mean:
imshow((abs(recon_stack(:,:,1))))
Warning: Image is too big to fit on screen; displaying at 33% 
> In images.internal.initSize (line 71)
  In imshow (line 305) 
imshow(abs(recon_stack(:,:,1)),[])
Warning: The initial magnification of the image is set to 'fit' in a docked figure. 
> In imshow (line 285) 
figure
imshow(abs(recon_stack(:,:,1)),[])
Warning: Image is too big to fit on screen; displaying at 33% 
> In images.internal.initSize (line 71)
  In imshow (line 305) 
hold on


for jkjk=1:iter-1
hold on
plot((coor_from{1,jkjk}(:,1)),(coor_from{1,jkjk}(:,2)),'g.')
plot((coor_to{1,jkjk}(:,1)),(coor_to{1,jkjk}(:,2)),'g.')
end
hold on
plot((coor_from{1,1}(:,1)),(coor_from{1,1}(:,2)),'bo')

figure
subplot(1,2,1),imshow((abs(imgrr1)+abs(recon_stack(imgrr2)))/2,[])
subplot(1,2,2),imshow((abs(recon_stack(:,:,5))+abs(recon_stack(:,:,6)))/2,[])
hold on
figure
subplot(1,2,1),imshow((abs(imgrr1)+abs(recon_stack(imgrr2)))/2,[])
subplot(1,2,2),imshow((abs(recon_stack(:,:,1))+abs(recon_stack(:,:,2)))/2,[])
hold on
subplot(1,2,2),plot((coor_from{1,1}(:,1)),(coor_from{1,1}(:,2)),'mo')
subplot(1,2,2),plot((coor_from{1,1}(:,1)),(coor_from{1,1}(:,2)),'m*')
hold on
subplot(1,2,2),plot((coor_to{1,1}(:,1)),(coor_to{1,1}(:,2)),'go')
subplot(1,2,2),plot((coor_to{1,1}(:,1)),(coor_to{1,1}(:,2)),'g*')
figure
figure
subplot(1,2,1),imshow((abs(imgrr1)+abs(recon_stack(imgrr2)))/2,[])
subplot(1,2,2),imshow((abs(recon_stack(:,:,1))+abs(recon_stack(:,:,2)))/2,[])
hold on
subplot(1,2,2),plot((coor_from{1,1}(:,1)),(coor_from{1,1}(:,2)),'m.')
hold on
subplot(1,2,2),plot((coor_to{1,1}(:,1)),(coor_to{1,1}(:,2)),'g.')
help plot

%------------differ image

min_recon_avg = min(min(rec_avg));
max_recon_avg = max(max(rec_avg));
rec_avg_flip = max_recon_avg-(rec_avg-min_recon_avg);
imshow(abs(rec_avg_flip),[])