
divide_ROI

mean(z)

mod(4,2)
z_new=NaN(1,Nb_ROI);

if mod(divide_ROI,2) == 0
    avrg_z=( z(((divide_ROI^2)/2) - (divide_ROI/2)) + z( ((divide_ROI^2)/2) - (divide_ROI/2) +1) + z( ((divide_ROI^2)/2) + (divide_ROI/2) ) + z( ((divide_ROI^2)/2) + (divide_ROI/2) + 1) )/4; 
    for i1=(divide_ROI-2)/2-1:-1:0
        
        jump1=(divide_ROI+1)*i1+1;
        record=z(jump1)-avrg_z;
        if i1== (divide_ROI-2)/2-1
            z_new(jump1)=z(jump1);
        else
            z_new(jump1)=z( (divide_ROI+1)*(i1+1)+1 ) + record;
        end
        
        jump2=(divide_ROI-1)*i1+divide_ROI;
        record=z(jump2)-avrg_z;
        if i1== (divide_ROI-2)/2-1
            z_new(jump2)=z(jump2);
        else
            z_new(jump2)=z( (divide_ROI-1)*(i1+1)+divide_ROI ) + record;
        end
        
        len=jump2-jump1-1;
        diff_1= (z_new(jump2)-z_new(jump1)) / (len+1);
        for jj=1:len
            z_new(jump1+jj)=
        end
    end
    
    for i2=(divide_ROI-2)/2-1:-1:0
        jump=(divide_ROI-1)*i2+divide_ROI;
        record=z(jump)-avrg_z;
        if i2== (divide_ROI-2)/2-1
            z_new(jump)=z(jump);
        else
            z_new(jump)=z( (divide_ROI-1)*(i2+1)+divide_ROI ) + record;
        end
    end
    
    
else
    avrg_z=ceil( (divide_ROI^2)/2 );
end

z_new=NaN(1,Nb_ROI);