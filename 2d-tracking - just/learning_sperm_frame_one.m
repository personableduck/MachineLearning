


test_point_one = [36 143; 113 138; 107 179; 126 197; 77 197; 142 165; 281 40; 305 37; 315 42; 261 23; 245 23];
inner_avg_one = zeros(1,size(test_point_one,1));
outer_avg_one = zeros(1,size(test_point_one,1));
for i = 1:size(test_point_one,1)
    y = test_point_one(i,1);
    x = test_point_one(i,2);
    inner_avg_one(1,i) = (abs(recon_stack(x,y,1)) + abs(recon_stack(x+1,y,1)) + abs(recon_stack(x-1,y,1))+ abs(recon_stack(x,y+1,1)) + abs(recon_stack(x,y-1,1)) + abs(recon_stack(x+1,y+1,1))+ abs(recon_stack(x+1,y-1,1)) + abs(recon_stack(x-1,y+1,1)) + abs(recon_stack(x-1,y-1,1)))/9;
    outer_avg_one(1,i) = (abs(recon_stack(x+2,y+2,1)) + abs(recon_stack(x+1,y+2,1)) + abs(recon_stack(x,y+2,1))+ abs(recon_stack(x-1,y+2,1)) + abs(recon_stack(x-2,y+2,1)) + abs(recon_stack(x-2,y+1,1))+ abs(recon_stack(x-2,y,1)) + abs(recon_stack(x-2,y-1,1)) + abs(recon_stack(x-2,y-2,1)) + abs(recon_stack(x-1,y-2,1)) + abs(recon_stack(x,y-2,1)) + abs(recon_stack(x+1,y-2,1))+ abs(recon_stack(x+2,y-2,1)) + abs(recon_stack(x+2,y-1,1)) + abs(recon_stack(x+2,y,1))+ abs(recon_stack(x+2,y+1,1)) )/16;
end
inner_avg_one_mean = mean(inner_avg_one);
outer_avg_one_mean = mean(outer_avg_one);
diff_avg_one = outer_avg_one_mean - inner_avg_one_mean;

test_point_two = [2346 1800; 2379 1803; 2392 1759; 2444 1759; 2327 1839];
inner_avg_two = zeros(1,size(test_point_two,1));
outer_avg_two = zeros(1,size(test_point_two,1));
for i = 1:size(test_point_two,1)
    y = test_point_two(i,1);
    x = test_point_two(i,2);
    inner_avg_two(1,i) = (abs(recon_stack(x,y,1)) + abs(recon_stack(x+1,y,1)) + abs(recon_stack(x-1,y,1))+ abs(recon_stack(x,y+1,1)) + abs(recon_stack(x,y-1,1)) + abs(recon_stack(x+1,y+1,1))+ abs(recon_stack(x+1,y-1,1)) + abs(recon_stack(x-1,y+1,1)) + abs(recon_stack(x-1,y-1,1)))/9;
    
    outer_avg_two(1,i) = (abs(recon_stack(x+2,y+2,1)) + abs(recon_stack(x+1,y+2,1)) + abs(recon_stack(x,y+2,1))+ abs(recon_stack(x-1,y+2,1)) + abs(recon_stack(x-2,y+2,1)) + abs(recon_stack(x-2,y+1,1))+ abs(recon_stack(x-2,y,1)) + abs(recon_stack(x-2,y-1,1)) + abs(recon_stack(x-2,y-2,1)) + abs(recon_stack(x-1,y-2,1)) + abs(recon_stack(x,y-2,1)) + abs(recon_stack(x+1,y-2,1))+ abs(recon_stack(x+2,y-2,1)) + abs(recon_stack(x+2,y-1,1)) + abs(recon_stack(x+2,y,1))+ abs(recon_stack(x+2,y+1,1)) )/16;
end
inner_avg_two_mean = mean(inner_avg_two);
outer_avg_two_mean = mean(outer_avg_two);
diff_avg_two = outer_avg_two_mean - inner_avg_two_mean;


test_point_three = [2279 180; 2316 209; 2285 138; 2229 170];
inner_avg_three = zeros(1,size(test_point_three,1));
outer_avg_three = zeros(1,size(test_point_three,1));
for i = 1:size(test_point_three,1)
    y = test_point_three(i,1);
    x = test_point_three(i,2);
    inner_avg_three(1,i) = (abs(recon_stack(x,y,1)) + abs(recon_stack(x+1,y,1)) + abs(recon_stack(x-1,y,1))+ abs(recon_stack(x,y+1,1)) + abs(recon_stack(x,y-1,1)) + abs(recon_stack(x+1,y+1,1))+ abs(recon_stack(x+1,y-1,1)) + abs(recon_stack(x-1,y+1,1)) + abs(recon_stack(x-1,y-1,1)))/9;
    
    outer_avg_three(1,i) = (abs(recon_stack(x+2,y+2,1)) + abs(recon_stack(x+1,y+2,1)) + abs(recon_stack(x,y+2,1))+ abs(recon_stack(x-1,y+2,1)) + abs(recon_stack(x-2,y+2,1)) + abs(recon_stack(x-2,y+1,1))+ abs(recon_stack(x-2,y,1)) + abs(recon_stack(x-2,y-1,1)) + abs(recon_stack(x-2,y-2,1)) + abs(recon_stack(x-1,y-2,1)) + abs(recon_stack(x,y-2,1)) + abs(recon_stack(x+1,y-2,1))+ abs(recon_stack(x+2,y-2,1)) + abs(recon_stack(x+2,y-1,1)) + abs(recon_stack(x+2,y,1))+ abs(recon_stack(x+2,y+1,1)) )/16;
end
inner_avg_three_mean = mean(inner_avg_three);
outer_avg_three_mean = mean(outer_avg_three);
diff_avg_three = outer_avg_three_mean - inner_avg_three_mean;


result = zeros(height, width);

for i = 3:height-2
   for j = 3:width-2
       y = j;
        x = i;
        inner_avg = (abs(recon_stack(x,y,1)) + abs(recon_stack(x+1,y,1)) + abs(recon_stack(x-1,y,1))+ abs(recon_stack(x,y+1,1)) + abs(recon_stack(x,y-1,1)) + abs(recon_stack(x+1,y+1,1))+ abs(recon_stack(x+1,y-1,1)) + abs(recon_stack(x-1,y+1,1)) + abs(recon_stack(x-1,y-1,1)))/9;
    
         outer_avg = (abs(recon_stack(x+2,y+2,1)) + abs(recon_stack(x+1,y+2,1)) + abs(recon_stack(x,y+2,1))+ abs(recon_stack(x-1,y+2,1)) + abs(recon_stack(x-2,y+2,1)) + abs(recon_stack(x-2,y+1,1))+ abs(recon_stack(x-2,y,1)) + abs(recon_stack(x-2,y-1,1)) + abs(recon_stack(x-2,y-2,1)) + abs(recon_stack(x-1,y-2,1)) + abs(recon_stack(x,y-2,1)) + abs(recon_stack(x+1,y-2,1))+ abs(recon_stack(x+2,y-2,1)) + abs(recon_stack(x+2,y-1,1)) + abs(recon_stack(x+2,y,1))+ abs(recon_stack(x+2,y+1,1)) )/16;
        result(x,y)=outer_avg-inner_avg;
   end
end

result_BW = result < 9;
figure
imshow(result_BW,[]);

