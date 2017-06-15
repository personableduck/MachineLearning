% min_intensities_xyz;
% 
% if(minIntensities(1)<minIntensities(currIdx))
%         centroidToAppend = [currPointX,currPointY];
%         centroidFromAppend = [centroidOrig(currIdx,1),centroidOrig(currIdx,2)];
%     else
%         centroidToAppend = [centroidOrig(currIdx,1),centroidOrig(currIdx,2)];
%         centroidFromAppend = [currPointX,currPointY];
% end
% 
%     centroidTo = [centroidTo; centroidToAppend];
%     centroidFrom = [centroidFrom; centroidFromAppend];
%     
%     centroidOrig(currIdx,:) = [];
%     minIntensities(currIdx,:) = [];
%     centroidOrig(1,:) = [];
%     minIntensities(1,:) = [];
%     
% end

%-------------------------------------------
    

divided_frame=1; %dividing fram number <- k
    
maxsize_num_xyz=size(xyz{1,divided_frame}); %read size
read_size_xyz=maxsize_num_xyz(1,1);

xyzs(:,1)=xyz{1,1}(:,1);
xyzs(:,2)=xyz{1,1}(:,2);
xyzs(:,3)=0;

for loop_xyz=1:(read_size_xyz)-1

% xyzs(loop_xyz,1)=xyz{1,divided_frame}(loop_xyz,1); %copy elements of xyz
% xyzs(loop_xyz,2)=xyz{1,divided_frame}(loop_xyz,2);
% 
% xyzs(read_size_xyz,1)=xyz{1,divided_frame}(read_size_xyz,1); %last xyz coordinates
% xyzs(read_size_xyz,2)=xyz{1,divided_frame}(read_size_xyz,2); %last xyz coordinates  
    
    
distance_threshold_positive = 15; %*****threshold object tracking. need to control. maximum ditance to distinguish as a same object


for realative_loop=loop_xyz+1:read_size_xyz
    
    xx = xyz{1,divided_frame}(loop_xyz,1) - xyz{1,divided_frame}(realative_loop,1); % comparing x coordinates
    yy = xyz{1,divided_frame}(loop_xyz,2) - xyz{1,divided_frame}(realative_loop,2); % comparing y coordinates
    
    if(sqrt(xx^2+yy^2) <= distance_threshold_positive) % finding moving objects, match function
        
        if (min_intensities_xyz{1,divided_frame}(change_size_xyz) < min_intensities_xyz{1,divided_frame}(loop_xyz)) % compare intensity to decide before one. making oreder between 2 frames
        
            xyzs(loop_xyz,3)= (divided_frame*2);
            xyzs(realative_loop,3)= (divided_frame*2)-1; 
        else       
            xyzs(loop_xyz,3)= (divided_frame*2)-1;
            xyzs(realative_loop,3)= (divided_frame*2);
        end
break
    
    end

end
end


%------

finding_frame_num=max(xyzs(:,3)); %read devided frame number

size_func_xyzs=size(xyzs);
read_size_xyzs=size_func_xyzs(1,1);


for devided_frame_num_loop=1:finding_frame_num
    
    count_increase=0; %increasing number to count number
    
if devided_frame_num_loop == (finding_frame_num-1) | devided_frame_num_loop == (finding_frame_num) %start from before loop exmaple) excute only 3,4 or 5,6

    for put_xyzs_num=1:read_size_xyzs
        
        if xyzs(put_xyzs_num,3) == devided_frame_num_loop
        count_increase=count_increase+1;
        devided_coordinates{1,devided_frame_num_loop}(count_increase,:) = xyzs(put_xyzs_num,:);
        
        end
        
        
        if xyzs(put_xyzs_num,3) == 0
        count_increase=count_increase+1;
        devided_coordinates{1,devided_frame_num_loop}(count_increase,:) = xyzs(put_xyzs_num,:);
        devided_coordinates{1,devided_frame_num_loop}(count_increase,3) = devided_frame_num_loop;

        end
            
    end
    

end
end


figure, imshow(recon_stack_diff(:,:,1),[])
hold on
plot(devided_coordinates{1,1}(:,1),devided_coordinates{1,1}(:,2),'g.')
hold on
plot(devided_coordinates{1,2}(:,1),devided_coordinates{1,2}(:,2),'b.')

% [temp,ord] = sort(xyzs(:,3));
% xyzs_ord = xyzs(ord,:);