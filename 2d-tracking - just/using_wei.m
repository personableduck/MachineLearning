
% figure, imshow(abs(recon_stack(:,:,1)),[])
% 
% for jkjk=1:8
% hold on
% plot((coor_mid{1,jkjk}(:,1)),(coor_mid{1,jkjk}(:,2)),'g.')
% end
% hold on
% plot((coor_mid{1,1}(:,1)),(coor_mid{1,1}(:,2)),'b.')


%----------------------------------
xyzs=0;
bb=0;


for k=1:iter-1 
    
a=size(xyz{1,k});
b=a(1,1);

for c=1:b
xyzs(c+bb,3)=k*2-1;
end

for jp=1:b-1

xyzs(jp+bb,1)=xyz{1,k}(jp,1);
xyzs(jp+bb,2)=xyz{1,k}(jp,2);

xyzs(b+bb,1)=xyz{1,k}(b,1);
xyzs(b+bb,2)=xyz{1,k}(b,2);
%------

threshold3 = 15; %threshold object tracking. need to control
threshold4=-(threshold3);

bk=b;

for jp2=1:bk-jp
    jp3=jp2+jp;
xx = xyz{1,k}(jp3,1) - xyz{1,k}(jp,1);
yy = xyz{1,k}(jp3,2) - xyz{1,k}(jp,2);
    
if((threshold4 <= xx & xx <= threshold3) & (threshold4 <= yy & yy <= threshold3))
xyzs(jp3+bb,3)=(k*2);       
break

end
end

end

%------


aa=size(xyzs);
bb=aa(1,1);  

end

[temp,ord] = sort(xyzs(:,3));
xyzs_ord = xyzs(ord,:);

trackk=track(xyzs_ord,20); %we have to find max distance 22 is max value

% figure, imshow(abs(recon_stack(:,:,1)),[]) 
% hold on
% plot((trackk(:,1)),(trackk(:,2)),'y.')

%----------------------------------------------------

jpy=max(trackk(:,4));

a1=size(trackk);
b1=a1(1,1);



for k2=1:jpy
    
    count_inc=0;
    
for j1=1:b1
    
    if trackk(j1,4) == k2
        count_inc=count_inc+1;
        object_trajectories{1,k2}(count_inc,:) = trackk(j1,:);
    end
    
end   
end   

figure, imshow(abs(recon_stack(:,:,1)),[]) 

for jljl=1:jpy
hold on
plot((object_trajectories{1,jljl}(:,1)),(object_trajectories{1,jljl}(:,2)),'g.')
hold on
plot((object_trajectories{1,jljl}(1,1)),(object_trajectories{1,jljl}(1,2)),'b.')
end

% set(gca, 'XLim', [80, 200], 'YLim', [70, 150])
% axis on
% 
% sqrt((122-126)^2+(114-115)^2)
% 
% %------
% jpy=max(xyzs(:,3));
% 
% a1=size(xyzs);
% b1=a1(1,1);
% 
% 
% 
% for k2=1:jpy
%     
%     count_inc=0;
%     
% for j1=1:b1
%     
%     if xyzs(j1,3) == k2
%         count_inc=count_inc+1;
%         object_trajectories{1,k2}(count_inc,:) = xyzs(j1,:);
%     end
%     
% end   
% end   
% 
% 
% 
% figure, imshow(recon_stack_diff(:,:,1),[])
% hold on
% plot((object_trajectories{1,1}(:,1)),(object_trajectories{1,1}(:,2)),'g.')
% hold on
% plot((object_trajectories{1,2}(:,1)),(object_trajectories{1,2}(:,2)),'b.')

% A=[1 2];
% B=[3 4];
% C=[A; B];