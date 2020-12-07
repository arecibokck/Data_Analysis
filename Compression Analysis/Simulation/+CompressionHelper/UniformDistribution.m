function ret = UniformDistribution(this, a, b, x)
    ret = zeros(length(x),1);
    ret(x>a & x<b) = 1/(b - a);
end