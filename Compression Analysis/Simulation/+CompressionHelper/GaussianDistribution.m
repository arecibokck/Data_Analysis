function ret = GaussianDistribution(mean, var, x)
    % return the gauss distribution with parameters mean, var.
    % p = gaussdis(mean, var, x);
    ret = 1/sqrt(2*pi*var)*exp(-((x-mean).^2)/(2*var));
    ret = ret./sum(ret);
end
