%% background image

% back_img=double(imread('C:\Users\1234\Downloads\d0-3.jpg'));
% 
% mask_value=mean(back_img(:))+2*std(back_img(:));
% object_img=back_img>mask_value;
% cover_img=object_img.*back_img;
% 
% % figure,imshow(object_img)
% 
% object_result = regionprops(object_img, 'Area', 'Centroid', 'BoundingBox');
% centroidOrig = cell2mat({object_result.Centroid}');
% box = cell2mat({object_result.BoundingBox}');
% 
% [width_1]=find(box(:,3)>180 & box(:,3)<250); %circle size filter 
% [height_2]=find(box(:,4)>180 & box(:,4)<250);
% result_row = intersect(width_1,height_2);
% 
% if (length(result_row) == 5)
%     display ('Success')
% else
%     display ('Fail: try to change circle size')
% end
% 
% figure,imshow(cover_img,[])
% hold on
% plot(centroidOrig(result_row(:),1),centroidOrig(result_row(:),2),'.')
% title('Background Entire image');
% 
% crop_size=140; 
% 
% crop_1=cover_img(centroidOrig(result_row(1),2)-crop_size:centroidOrig(result_row(1),2)+crop_size ,centroidOrig(result_row(1),1)-crop_size: centroidOrig(result_row(1),1)+crop_size);
% 
% for i=1:length(result_row)
%     crop_img{1,i}=cover_img( (centroidOrig(result_row(i),2)-crop_size):(centroidOrig(result_row(i),2)+crop_size) ,(centroidOrig(result_row(i),1)-crop_size): (centroidOrig(result_row(i),1)+crop_size));  
%     result_mean{1,i}= crop_img{1,i}(crop_img{1,i}~=0);
%     
%     figure,imshow(crop_img{1,i})
%     title(['Croped image ', num2str(i)]);
%     xlabel(['Mean intensity: ', num2str(mean(result_mean{1,i}(:)))]);
% 
% end    


%% strong response image
back_img=double(imread('C:\Users\1234\Downloads\d0.5-3.jpg'));
back_img=double(imread('C:\Users\1234\Downloads\d0-3.jpg'));
comp_row=0

mask_value=mean(back_img(:))+2*std(back_img(:));
object_img=back_img>mask_value;
cover_img=object_img.*back_img;

figure,imshow(back_img,[])

object_result = regionprops(object_img, 'Centroid','BoundingBox');
centroidOrig = cell2mat({object_result.Centroid}');
box = cell2mat({object_result.BoundingBox}');

[width_1]=find(box(:,3)>130 & box(:,3)<180); %minor circle size filter 
[height_1]=find(box(:,4)>130 & box(:,4)<180);
result_row = intersect(width_1,height_1);


if (length(result_row) < 5)
    [width_2]=find(box(:,3)>180 & box(:,3)<250); %major circle size filter 
    [height_2]=find(box(:,4)>180 & box(:,4)<250);
    comp_row = intersect(width_2,height_2);
end    

check_row=comp_row;

if (comp_row > 0)
    for j=1:length(comp_row)
        for i=1:length(result_row)
            if (norm( [ centroidOrig(result_row(i),1) - centroidOrig(comp_row(j),1), ...
                    centroidOrig(result_row(i),2) - centroidOrig(comp_row(j),2) ] ) <= 60) 
                check_row(j)=nan;
            end
        end
    end
    
    add_row=[];
    increase_num=0;
    for i=1:length(check_row)
        if ~isnan(check_row(i))
            increase_num=increase_num+1;
            add_row(increase_num)= check_row(i);
        end    
    end

    result_row=[result_row;add_row];
end

if (length(result_row) == 5)
    display ('Success')
else
    display ('Fail: try to change circle size')
end

figure,imshow(cover_img,[])
hold on
plot(centroidOrig(result_row(:),1),centroidOrig(result_row(:),2),'.')
title('Background Entire image');

imageSize = size(back_img);

for i=1:length(result_row)
    
    ci = [centroidOrig(result_row(i),2), centroidOrig(result_row(i),1), 73];     % center and radius of circle ([c_row, c_col, r])
    [xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
    mask = ((xx.^2 + yy.^2)<ci(3)^2);
    cover_img=mask.*back_img;  
    average_cal= cover_img(cover_img~=0);
    
    figure,imshow(cover_img,[])
    title(['Croped image ', num2str(i)]);
    xlabel(['Mean intensity: ', num2str(mean(average_cal(:)))]);

end


