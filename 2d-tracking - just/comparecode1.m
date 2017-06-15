figure
imshow(abs(recon_stack(:,:,9)),[])%amplitude
set(gca, 'XLim', [170, 330], 'YLim', [190, 360])
axis on

for jkjk=1:iter-1
hold on
plot((coor_from{1,jkjk}(:,1)),(coor_from{1,jkjk}(:,2)),'g.')
plot((coor_to{1,jkjk}(:,1)),(coor_to{1,jkjk}(:,2)),'g.')
end
hold on
plot((coor_from{1,1}(:,1)),(coor_from{1,1}(:,2)),'bo')





hold on
for jseed=1:233
hold on
plot((x_obj(:,jseed)),(y_obj(:,jseed)),'g')
end
hold on
plot((x_obj(5,:)),(y_obj(5,:)),'b.')
