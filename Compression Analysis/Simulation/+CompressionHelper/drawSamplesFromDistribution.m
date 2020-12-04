function ret = drawSamplesFromDistribution(NumberOfSamples, Xvals, ProbabilityDensityFunction)
    % ret = drawSamplesFromDistribution(NumberOfAtoms, velocities, ProbabilityDensityFunction)
    % ret is a N*1 array of randomly sampled values from a distribution defined by P(X) 
    %
    % input: NumberOfSamples   
    %        Xvals   
    %        ProbabilityDensityFunction   
    %
    % output:  ret    N*1 array of randomly sampled values
    %
    % Note:  this implements "Inverse transform sampling"
    
    % - normalize in L1-norm
    ProbabilityDensityFunction     = ProbabilityDensityFunction./sum(ProbabilityDensityFunction); 
    % - compute cumulative distribution
    CumulativeDistributionFunction = cumsum(ProbabilityDensityFunction); 
    % - append a 0 at the beginning of the distribution if the input starts with a non-zero entry (e.g. when sampling a gaussian-distribution)
    if(CumulativeDistributionFunction(1)>0) 
        CumulativeDistributionFunction = [0 CumulativeDistributionFunction];
        dx = Xvals(2) - Xvals(1);
        x0 = Xvals(1) - dx;
        Xvals = [x0 Xvals];
    end 
    % - remove repeted values
    [CumulativeDistributionFunction, index] = unique(CumulativeDistributionFunction); 
    
    % - create uniformly distributed random numbers between 0 and 1
    rnd = rand(NumberOfSamples, 1); 
    % - calculate output
    ret = interp1(CumulativeDistributionFunction, Xvals(index), rnd, 'linear', 0);
end