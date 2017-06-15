imgrr3=imread('ICI-11642v2_x3-1_0.0ms.bmp');
imshow(imgrr3)
set(gca, 'XLim', [170, 330], 'YLim', [190, 360])
axis on

imgrr3=double(imread('ICI-11642v2_x3-1_986.0ms.bmp'));

imgrr1=imread('ICI-11642v2_x3-1_0.0ms.bmp');
imgrr2=imread('ICI-11642v2_x3-1_496.0ms.bmp');
imgrr3=imread('ICI-11642v2_x3-1_986.0ms.bmp');
imgrr4=imread('ICI-11642v2_x3-1_1441.0ms.bmp');
imgrr5=imread('ICI-11642v2_x3-1_1896.0ms.bmp');
imgrr6=imread('ICI-11642v2_x3-1_2352.0ms.bmp');
imgrr7=imread('ICI-11642v2_x3-1_2814.0ms.bmp');
imgrr8=imread('ICI-11642v2_x3-1_3260.0ms.bmp');
imgrr9=imread('ICI-11642v2_x3-1_3720.0ms.bmp');

figure
imshow(imgrr1)
imshow(imgrr1+imgrr2+imgrr3+imgrr4+imgrr5+imgrr6+imgrr7+imgrr8+imgrr9)
imshow((imgrr1+imgrr2+imgrr3+imgrr4+imgrr5+imgrr6+imgrr7+imgrr8+imgrr9)/2,[])


%amplitude
set(gca, 'XLim', [170, 330], 'YLim', [190, 360])
axis on

figure
imshow(abs(recon_stack(:,:,1)),[])
imshow(angle(recon_stack(:,:,1)),[]) %phase image

imshow(abs(recon_stack(:,:,3))-imgrr3,[])


I = imread('rice.png');
Iq = imsubtract(tryy,50);
figure, imshow(I), figure, imshow(Iq)

a123=abs(recon_stack(:,:,1));
a234=abs(recon_stack(:,:,1)); %different z

figure
imshow(a123-a234)


figure, imshow(abs(recon_stack(:,:,1))+abs(recon_stack(:,:,2))+abs(recon_stack(:,:,3))+abs(recon_stack(:,:,4))+abs(recon_stack(:,:,5))+abs(recon_stack(:,:,6))+abs(recon_stack(:,:,7))+abs(recon_stack(:,:,8))+abs(recon_stack(:,:,9)),[])


figure, imshow(a45a5-a45aaa)
a45a5=abs(recon_stack(:,:,1))+abs(recon_stack(:,:,2))+abs(recon_stack(:,:,3))+abs(recon_stack(:,:,4))+abs(recon_stack(:,:,5))+abs(recon_stack(:,:,6))+abs(recon_stack(:,:,7))+abs(recon_stack(:,:,8))+abs(recon_stack(:,:,9));
a45aaa=abs(recon_stack(:,:,2))+abs(recon_stack(:,:,3))+abs(recon_stack(:,:,4))+abs(recon_stack(:,:,5))+abs(recon_stack(:,:,6))+abs(recon_stack(:,:,7))+abs(recon_stack(:,:,8))+abs(recon_stack(:,:,9));

figure, imshow(((abs(recon_stack(:,:,1))+abs(recon_stack(:,:,2))+abs(recon_stack(:,:,3))+abs(recon_stack(:,:,4))+abs(recon_stack(:,:,5))+abs(recon_stack(:,:,6))+abs(recon_stack(:,:,7))+abs(recon_stack(:,:,8))+abs(recon_stack(:,:,9)))-(abs(recon_stack(:,:,1))+abs(recon_stack(:,:,2))+abs(recon_stack(:,:,3))+abs(recon_stack(:,:,4))+abs(recon_stack(:,:,5))+abs(recon_stack(:,:,6))+abs(recon_stack(:,:,7))+abs(recon_stack(:,:,8))+abs(recon_stack(:,:,9)))/1.3),[])




%%%%%%%%%
a34346=uint8(abs(recon_stack(:,:,3)));

%%%%%%%%%%%%%%%%%%%%

hold on
plot((centroidOrig(:,1)),(centroidOrig(:,2)),'.b')



