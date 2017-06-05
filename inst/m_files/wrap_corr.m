function SigHat = wrap_corr(m,Ne,cummap,S,cutoff)
% USAGE: compute the shrinkage estimator of covariance matrix in Wen and Stephens (2010)
% INPUT:
%	m: the number of individuals in the reference panel, integer
%	Ne: the effective population size (diploid), integer
%	cummap: cumulative genetic map in cM, numSNP by 1
%	Hpanel: the estimated covariance matrix (cov(Hpanel))
%	cutoff: the hard threshold for small entries being zero, scalar 
% OUTPUT:
%	SigHat: estimated covariance matrix of haplotype, numSNP by
%	numSNP

nmsum = sum(1 ./ (1:(2*m-1)));
theta = (1/nmsum) / (2*m + 1/nmsum);

S = triu(S);
numSNP = size(S,1);
numSNP = size(S,1);
for i = 1:numSNP
    for j = (i+1):numSNP
        if cummap(j) < cummap(i)
            error('Negative recombination rate is produced.')
        end
        rho = 4 * Ne * (cummap(j) - cummap(i)) / 100;
        shrinkage = exp(-rho/(2*m));
        % hard thresholding to obtain sparse and banded estimator
        if shrinkage < cutoff
            shrinkage = 0;
        end
        S(i, j) = shrinkage * S(i, j);
    end
end

  % copy the upper half (NO diagonal) to the lower half
  S = S + triu(S, 1)';
	
  % SigHat is derived from Li and Stephens model (2003)
  SigHat = (1-theta)^2 * S + 0.5*theta * (1-0.5*theta) * eye(numSNP);
end