function ret = MirroredFlatTopGaussianDistribution(offset, flatWidth, stddev, x)
    % ret = MirroredFlatTopGaussianDistribution(offset, flatWidth, stddev, x)
    
    ret = exp(-max((abs(x-offset)-flatWidth),0).^2./(2*stddev^2));
    % normalize
    ret = ret./sum(ret);
end
