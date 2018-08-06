function FDR = FDRAnalysis(Corrs,PermCorrs)

% function FDR = FDRAnalysis(Corrs,PermCorrs) performs a false discovery
% rate analysis of m-by-m correlation matrix Corrs, in comparison with the
% permutation-based m-by-m-by-Repeats correlation matrix PermCorrs. It
% computes how many above-diagonal entries of Corrs are above different
% possible correlation thresholds (namely -1:0.01:1). It looks at the
% empirical fraction of permuted correlations above each of those
% thresholds, on the basis of which it computes expected numbers of false
% positives. And from that, it estimates false discovery rate. These
% information provide guidance to the user for selecting a correlation
% cutoff to form a relevance network.

[m,dummy1] = size(Corrs);
[dummy2,dummy3,Repeats] = size(PermCorrs);

% What is the set of correlation thresholds that we will test?
FDR.CorrThresh = -1:0.01:1;

% For each threshold...
for i=1:length(FDR.CorrThresh)
    % Threshold
    CT = FDR.CorrThresh(i);
    % How many above-diagonal entries are above that threshold?
    FDR.NCorrAbove(i) = (sum(sum(Corrs>=CT))-m)/2;
    % What fraction of permutations are above that threshold?
    FDR.FPermAbove(i) = (sum(sum(sum(PermCorrs>=CT)))-m*Repeats)/(Repeats*m*(m-1));
    % Expected false positives
    FDR.EFalsePos(i) = FDR.FPermAbove(i)*m*(m-1)/2;
    % Estimated false discovery rate
    FDR.EFDR(i) = FDR.EFalsePos(i)/FDR.NCorrAbove(i);
    if FDR.EFDR(i)>1
        FDR.EFDR(i) = 1;
    end
end



