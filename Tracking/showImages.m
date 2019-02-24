function showImages(meas)
    figure(1);
    a_h=gca;
    axis square;

    for i = 1:37
        imagesc(a_h,squeeze(meas.runData.images(end,i,:,:)));
        pause(1);
    end
end