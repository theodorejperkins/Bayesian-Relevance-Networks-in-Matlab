function Corrs_mmr = BayesianPermutation_Grouped(ReadCounts,Groups,Repeats,SkipCovU)

% function Corrs_mmr = BayesianPermutation_Grouped(ReadCounts,Groups,Repeats,SkipCovU)
% estimates a null distribution for corrleations computed by the function 
% BayesianCorrelation_Grouped. The inputs are ReadCounts (m-by-n matrix),
% Groups (1-by-n vector), Repeats, a positive integer, and SkipCovU. The 
% ReadCounts, Groups and SkipCovU inputs have the same meaning as for the
% function BayesianCorrelation_Grouped. The Repeats input specifies the
% number of random permutations to test. The output is an m-by-m-by-Repeats
% matrix of correlations computed from permutations. These can be used as
% estimates of null distributions for each pair of entities, or can be
% combined to form a single, overall null distribution.

if nargin<4
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

% Permuted covariance and correlation
S_m = 1./sqrt(VarGU_m);
SS_mm = S_m*S_m';
Corrs_mmr = zeros(m,m,Repeats);
for r=1:Repeats
    disp(['Repeat ',num2str(r)]);
    % Permute one copy of X
    TempX = X;
    for i=1:m
        TempX(i,:) = TempX(i,randperm(g));
    end
    % Permuted covariance
    CovGMeanU_mm = X*TempX'/g;
    CovGU_mm = CovGMeanU_mm + MeanGCovU_mm;
    % Permuted correlation
    Corrs = CovGU_mm .* SS_mm;
    Corrs(1:(m+1):m^2) = 1;
    Corrs_mmr(:,:,r) = Corrs;
end




