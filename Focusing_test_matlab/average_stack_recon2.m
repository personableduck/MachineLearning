z2_list=374;
pixelsize=1.12;

av_recon = angle(Prop_SSA(avg_holo,1.12,1.12,-z2_list,0.625));
figure,imshow(av_recon,[])

for i=-10:10
    if i == -10
        av_recon = angle(BP(avg_holo, avrg_z+0.1*i));
        object=mean(av_recon(:)) - std(av_recon(:)) > av_recon;
%         figure,imshow(object,[])
        av_recon_stack=object .* av_recon;
%         figure,imshow(abs(av_recon_stack),[])
    else
        av_recon = angle(BP(avg_holo, avrg_z+i));
        object=mean(av_recon(:)) - std(av_recon(:)) > av_recon;
        av_recon_stack=object .* av_recon;
        
        av_recon_stack=abs(av_recon_stack) .* abs(av_recon);
    end
end

figure,imshow(av_recon_stack,[])
