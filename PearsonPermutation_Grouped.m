function Corrs_mmr = PearsonPermutation_Grouped(ReadCounts,Groups,Repeats)

% function Corrs_mmr = PearsonPermutation_Grouped(ReadCounts,Groups,Repeats)
% estimates a null distribution for corrleations computed by the function 
% PearsonCorrelation_Grouped. The inputs are ReadCounts (m-by-n matrix),
% Groups (1-by-n vector), and repeats, a positive integer. The ReadCounts
% and Groups inputs have the same meaning as for the function 
% PearsonCorrelation_Grouped. The final input specifies the number of
% random permutations to test. The output is an m-by-m-by-Repeats matrix of
% correlations computed from permutations. These can be used as estimates
% of null distributions for each pair of entities, or can be combined to
% form a single, overall null distribution.

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

% Permuted covariances and correlations
Corrs_mmr = zeros(m,m,Repeats);
rng('default');
for r=1:Repeats
    disp(['Repeat ',num2str(r)]);
    % Permute X temporarily
    TempX_mg = X_mg;
    for i=1:m
        TempX_mg(i,:) = TempX_mg(i,randperm(g));
    end
    % Permuted entity covariances across groups
    TempCov = X_mg*(TempX_mg')/g;
    % Correlations
    TempCorrs = TempCov .* OneByStd_mm;
    TempCorrs(1:(m+1):(m^2)) = 1;
    Corrs_mmr(:,:,r) = TempCorrs;
end
I = find(isnan(Corrs_mmr));
if ~isempty(I)
    Corrs_mmr(I) = 0;
end



