%next process increse for slow --> 400 % // incease rate 0.65 

%%%%%%SPEED

disfps= 10; %fps
timess= 2814/1000;  %fps*iter;

dkcormax= 1035;

spsp=0;

for ababs = 1: dkcormax
sumdis=0;

max12345=size(track_coor(:,:,ababs))

if(max12345(1,1)>1)

for objds = 1:max12345(1,1)-1
    
    obdistances= pixelsize*(sqrt((track_coor(objds+1,1,ababs)-track_coor(objds,1,ababs))^2 + (track_coor(objds+1,2,ababs)-track_coor(objds,2,ababs))^2)); % um
    
    sumdis = obdistances + sumdis;
    
end


speedss = sumdis/timess; %um/s
spsp=spsp+1;

if(obj_spd(spsp,1) < 2000)

    obj_spd(spsp,1) = speedss;
else
    obj_spd(spsp,1) = 0;

end

end

end

figure
hist(obj_spd,25)
xlabel('Speed(um/sec)')
ylabel('Counts')
median(obj_spd)
mean(obj_spd)
dkcormax %objectnumber