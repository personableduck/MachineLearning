function new_factor = find_factor(matrix_xy,size_mt,threshold_mt,k)

b = size_mt;

xyzs = matrix_xy;

threshold3 = threshold_mt; %threshold object tracking
threshold4=-(threshold3);

for jp=1:b-1
    for jp2=1:b-1
        jp3=jp2+1;
        xx = xyzs(jp3,1) - xyzs(jp,1);
        yy = xyzs(jp3,2) - xyzs(jp,2);
    
        if((threshold4 <= xx & xx <= threshold3) & (threshold4 <= yy & yy <= threshold3))
        xyzs(jp3,3)=(k*2); 
        end
        
    end
    
end


new_factor = new_xy;