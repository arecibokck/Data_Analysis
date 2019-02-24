function showImages(meas)
    figure(1);
    a_h=gca;
    axis square;

    for i = 1:5
        image = squeeze(meas.runData.images(i,1,:,:));
        image_size = size(image,1);
        sector_size = floor(image_size/3);
        l = [1:1:image_size];
        l_1 = repmat(sector_size, 1, image_size);
        l_2 = repmat(2*sector_size, 1, image_size);
        imagesc(a_h, image);
        xlim([1 image_size]);
        ylim([1 image_size]);
        plot(l, l_1, 'LineStyle','-','Color',[1.,1,1,0.6],'LineWidth',2);
        hold on;
        plot(l, l_2, 'LineStyle','-','Color',[1.,1,1,0.6],'LineWidth',2);
        hold on;
        plot(l_1, l, 'LineStyle','-','Color',[1.,1,1,0.6],'LineWidth',2);
        hold on;
        plot(l_2, l, 'LineStyle','-','Color',[1.,1,1,0.6],'LineWidth',2);
        pause(1);
    end
end