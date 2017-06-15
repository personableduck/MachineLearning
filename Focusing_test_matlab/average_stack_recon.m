figure,imshow(avg_holo,[])
avrg_z = 382.7267;

BP = @(x,z) Propagate(x, 1.12, 1, 633, z, true, false, true);
av_recon = angle( BP(avg_holo, avrg_z));
figure,imshow(av_recon,[])

for i=-10:10
    if i == -10
        av_recon = angle(BP(avg_holo, avrg_z+i));
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

