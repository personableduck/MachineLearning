function [ imageInterp ] = InterpSonyImage( image )
    
    imageInterp = image;
    cordH = 3:2:size(image,2)-2;
    cordV = 3:2:size(image,1)-2;
    imageInterp(cordV, cordH) = 0.25*(image(cordV + 1, cordH) + image(cordV, cordH + 1) + ...
        image(cordV - 1, cordH) + image(cordV, cordH - 1));
    cordH = 4:2:size(image,2)-3;
    cordV = 4:2:size(image,1)-3;
    imageInterp(cordV, cordH) = 0.25*(image(cordV + 1, cordH) + image(cordV, cordH + 1) + ...
        image(cordV - 1, cordH) + image(cordV, cordH - 1));
    
end

