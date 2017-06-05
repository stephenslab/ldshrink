function S = compute_shrinkage(m, Ne,Hpanel, cummap, cutoff)
% USAGE: compute the shrinkage estimator of covariance matrix in Wen and Stephens (2010)
% INPUT:
%	m: the number of individuals in the reference panel, integer
%	Ne: the effective population size (diploid), integer
%	cummap: cumulative genetic map in cM, numSNP by 1
%	Hpanel: the (phased) haplotypes from a reference panel, numIND by numSNP
%	cutoff: the hard threshold for small entries being zero, scalar 
% OUTPUT:
%	SigHat: estimated covariance matrix of haplotype, numSNP by numSNP

% theta is related to mutation
  nmsum = sum(1 ./ (1:(2*m-1)));
  theta = (1/nmsum) / (2*m + 1/nmsum);
	
  % S is obtained from Sigma_panel by shrinking off-diagonal entries toward 0
  % NB: Hpanel is a phased haplotype matrix
  S = cov(Hpanel);
  S = triu(S);


  % compute the values on and above the 1st diagonal
  numSNP = size(S,1);
  shrinkage=zeros(numSNP,numSNP);
  for i = 1:numSNP
      for j = (i+1):numSNP
          if cummap(j) < cummap(i)
              error('Negative recombination rate is produced.')
          end
          rho = 4 * Ne * (cummap(j) - cummap(i)) / 100;
          shrinkage(i,j) = exp(-rho/(2*m));
          % hard thresholding to obtain sparse and banded estimator
          if shrinkage(i,j) < cutoff
              shrinkage(i,j) = 0;
          end
          S(i, j) = shrinkage(i,j) * S(i, j);
      end
  end
end