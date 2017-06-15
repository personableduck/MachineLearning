t = [1.5, 3, 5, 8]; %each exposure time
a = [133, 221, 362, 414]; %each mean value for green channel
figure; plot(t,a, '.--r')
figure; plot(t,a./t, '.--r')
figure; plot(t,(a-60)./t, '.--r') %dark current supression is 60
