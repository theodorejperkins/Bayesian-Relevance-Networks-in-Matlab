function Corrs = PearsonCorrelation_Grouped(ReadCounts,Groups)

% function Corrs = PearsonCorrelation_Grouped(ReadCounts,Groups) computes
% the grouped Pearson correlations between all pairs of m entities across n
% conditions. The m-by-n ReadCounts input matrix specifies the numbers of
% reads for each entity (rows) and condition (columns). The second input,
% Groups, is a 1-by-n vector of group numbers, specifying to which group
% each condition belongs. For instance, if the first two conditions are
% group 1, second three conditions are group 3, and third three conditions
% are group 2, we would have Groups = [1 1 3 3 3 2 2 2]. The correlation is
% computed by normalizing read counts by the total reads in each column,
% then averaging those within groups, and computing the correlation across
% groups. The answer is return int he m-by-m matrix Corrs.

% Sizes of things
[m,n] = size(ReadCounts);
UniqGroups = unique(Groups);
g = length(UniqGroups);

% Empirical frations
F_mn = zeros(m,n);
for j=1:n
    F_mn(:,j) = ReadCounts(:,j)/sum(ReadCounts(:,j));
end

% Group means
MeanS_mg = zeros(m,g);
for i=1:g
    I = find(Groups==UniqGroups(i));
    MeanS_mg(:,i) = mean(F_mn(:,I),2);
end

% Entity means across groups
MeanGS_m = mean(MeanS_mg,2);

% Entity variances across groups
X_mg = MeanS_mg - repmat(MeanGS_m,1,g);
VarG_m = mean(X_mg.^2,2);
OneByStd_m = 1./sqrt(VarG_m);
OneByStd_mm = OneByStd_m*(OneByStd_m');

% Entity covariances across groups
CovG_mm = X_mg*(X_mg')/g;

% Correlations
Corrs = CovG_mm .* OneByStd_mm;
Corrs(1:(m+1):(m^2)) = 1;
I = find(isnan(Corrs));
if ~isempty(I)
    Corrs(I) = 0;
end




