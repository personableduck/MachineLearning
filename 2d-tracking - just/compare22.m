set(gca, 'XLim', [2350, 2560], 'YLim', [1600, 1780])

figure
imshow((abs(recon_stack(:,:,1))+abs(recon_stack(:,:,2)))/2,[])

hold on
for jseed=1:688
hold on
plot((x_obj(:,jseed)),(y_obj(:,jseed)),'g')
end
hold on
plot((x_obj(5,:)),(y_obj(5,:)),'b.')

figure
hist(2.2*speed_avg,35)
xlabel('Speed(um/sec)')
ylabel('Counts')

median(2.2*speed_avg)
mean(2.2*speed_avg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
imshow((abs(recon_stack(:,:,1))+abs(recon_stack(:,:,2)))/2,[])
hold on
for pictureN=1:dkjj
hold on
plot((obj_mat_cor{pictureN,1}(:,1)),(obj_mat_cor{pictureN,1}(:,2)),'g.')
hold on
plot((obj_mat_cor{pictureN,1}(16,1)),(obj_mat_cor{pictureN,1}(16,2)),'b.')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hist(speed_all,25)
