function Corrs = BayesianCorrelation_Grouped(ReadCounts,Groups,SkipCovU)

% function Corrs_mmr = BayesianPermutation_Grouped(ReadCounts,Groups,Repeats,SkipCovU)
% computes the grouped Bayesian correlations between all pairs of m
% entities across n conditions. The m-by-n ReadCounts input matrix
% specifies the numbers of reads for each entity (rows) and condition
% (columns). The second input, Groups, is a 1-by-n vector of group numbers,
% pecifying to which group each condition belongs. For instance, if the
% first two conditions are group 1, second three conditions are group 3,
% and third three conditions are group 2, we would have Groups = [1 1 3 3
% 3 2 2 2]. The optional 3rd argument, if true, tells the function not to
% compute the second half of the covariance term (covariance of
% uncertainties in levels, averaged across conditions & groups). The term
% tends to be very small compared to everything else, yet slow to compute.
% Further, it is unchanged in permutation computations, so in most cases,
% it can be safely ignored. The answer is returned in the m-by-m matrix
% Corrs.

if nargin<3
    SkipCovU = 0;
end

% How big are things?
[m,s] = size(ReadCounts);
UniqGroups = unique(Groups);
g = length(UniqGroups);

% Priors
PriorAlphas_ms = ones(m,s)/(m-1);

% Compute posteriors and concentration parameters (sums of alphas across
% entities)
PosteriorAlphas_ms = ReadCounts + PriorAlphas_ms;
TotalAlphas_s = sum(PosteriorAlphas_ms,1);
TotalAlphas_ms = repmat(TotalAlphas_s,m,1);

% Posterior means and variances by sample
MeanU_ms = PosteriorAlphas_ms./TotalAlphas_ms;
VarU_ms = PosteriorAlphas_ms.*(TotalAlphas_ms-PosteriorAlphas_ms)./(TotalAlphas_ms.*TotalAlphas_ms.*(TotalAlphas_ms+1));

% Posterior means and variances by group
MeanU_mg = [];
VarU_mg = [];
for i=1:g
    I = find(Groups==UniqGroups(i));
    MeanU_mg = [MeanU_mg mean(MeanU_ms(:,I),2)];
    X = MeanU_ms(:,I)-repmat(MeanU_mg(:,end),1,length(I));
    VarU_mg = [VarU_mg (mean(VarU_ms(:,I),2)+mean(X.^2,2))/length(I)];
end

% Variance across groups
MeanGMeanU_m = mean(MeanU_mg,2);
X = MeanU_mg-repmat(MeanGMeanU_m,1,g);
VarGMeanU_m = mean(X.^2,2);

% Total variance
MeanGVarU_m = mean(VarU_mg,2);
VarGU_m = VarGMeanU_m + MeanGVarU_m;

% Mean across groups of covariance within each group
MeanGCovU_mm = zeros(m);
if ~SkipCovU
    for i=1:g
        TempCov = zeros(m);
        J = find(Groups==UniqGroups(i));
        for j=1:length(J)
            s = J(j);
            TempCov = TempCov - MeanU_ms(:,s)*MeanU_ms(:,s)'/(TotalAlphas_s(s)+1);
        end
        TempCov = TempCov/length(J);
        MeanGCovU_mm = MeanGCovU_mm + TempCov;
    end
    MeanGCovU_mm = MeanGCovU_mm/g;
end

% Covariance across groups of mean within each group
CovGMeanU_mm = X*X'/g;

% Total covariance
CovGU_mm = CovGMeanU_mm + MeanGCovU_mm;

% Correlation
S_m = 1./sqrt(VarGU_m);
SS_mm = S_m*S_m';
Corrs = CovGU_mm .* SS_mm;
Corrs(1:(m+1):m^2) = 1;


