function ret = HigherOrderGaussianDistribution(mean, var, exponent, x)
    ret = exp(-((x-mean).^2./(2*var)).^exponent);
    ret = ret./sum(ret);
end
