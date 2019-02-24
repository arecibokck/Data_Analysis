function secMWSBS(measData)
    peak_find_thres = 200;
    for r = 1:size(measData.runData.images,1)
        
        first_image = squeeze(measData.runData.images(r,1,:,:));
        first_image_peaks = FastPeakFind(first_image, peak_find_thres);
        first_im_count = sectorCount(first_image_peaks, size(first_image,1));
        
        second_image = squeeze(measData.runData.images(r,2,:,:));
        second_image_peaks = FastPeakFind(second_image, peak_find_thres);
        second_im_count = sectorCount(second_image_peaks, size(second_image,1));
        
        diff_image = squeeze(measData.runData.images(r,3,:,:));
        diff_image_peaks = FastPeakFind(diff_image, peak_find_thres);
        diff_im_count = sectorCount(diff_image_peaks, size(diff_image,1));
        
        result = first_im_count - diff_im_count ./ second_im_count - diff_im_count;
        
    end    
end


function atomCounts = sectorCount(peakData, image_size)
    sector_size = floor(image_size/3);
    atomCounts = zeros(1,9);
    
    x = peakData(1:2:end);
    y = peakData(2:2:end);
    
    for i = 1:length(x)
        for x_i = [0, sector_size, 2*sector_size]
            for y_i = [0, sector_size, 2*sector_size]
                if x_i == 0 && y_i == 0
                    t = 1;
                elseif x_i == 0 && y_i == sector_size
                    t = 2;
                elseif x_i == 0 && y_i == 2*sector_size
                    t = 3;
                elseif x_i == sector_size && y_i == 0
                    t = 4;
                elseif x_i == sector_size && y_i == sector_size
                    t = 5;
                elseif x_i == sector_size && y_i == 2*sector_size
                    t = 6;
                elseif x_i == 2*sector_size && y_i == 0
                    t = 7;
                elseif x_i == 2*sector_size && y_i == sector_size
                    t = 8;
                elseif x_i == 2*sector_size && y_i == 2*sector_size
                    t = 9;
                end
                if x(i) >= (x_i+1) && x(i) < (x_i+sector_size) && y(i) >= (y_i+1) && y(i) < (y_i+sector_size)
                    atomCounts(t) = atomCounts(t) + 1;
                end
            end
        end
    end
end