function [R] = get_corr_d(m, Ne, cummap, Hpanel, cutoff)
% USAGE: compute LD matrix using the shrinkage estimator in Wen and Stephens (2010)
% INPUT:
%	m: the number of individuals in the reference panel, integer
%	Ne: the effective population size (diploid), integer
%	cummap: cumulative genetic map in cM, numSNP by 1
%	Hpanel: the (phased) haplotypes from a reference panel, numIND by numSNP
%	cutoff: the hard threshold for small entries being zero, scalar 
% OUTPUT:
%	R: the estimated LD matrix, numSNP by numSNP, sparse matrix
%	BR: the banded storage of R

  % compute the shrinkage estimator of covariance matrix 
  SigHat = shrink_cov(m, Ne, cummap, Hpanel, cutoff);

  % convert covariance to correlation
  R = corrcov(SigHat);
  clear SigHat;

  % get the bandwidth of R
end



