function [p,f_obs,f_perm] = rm_anova2_np(Y,S,F1,F2,nperm)
%
% function stats = rm_anova2(Y,S,F1,F2,FACTNAMES)
%
% Two-factor, within-subject repeated measures ANOVA.
% For designs with two within-subject factors.
%
% Parameters:
%    Y          dependent variable (numeric) in a column vector
%    S          grouping variable for SUBJECT
%    F1         grouping variable for factor #1
%    F2         grouping variable for factor #2
%    FACTNAMES  a cell array w/ two char arrays: {'factor1', 'factor2'}
%
%    Y should be a 1-d column vector with all of your data (numeric).
%    The grouping variables should also be 1-d numeric, each with same
%    length as Y. Each entry in each of the grouping vectors indicates the
%    level # (or subject #) of the corresponding entry in Y.
%
% Returns:
%    stats is a cell array with the usual ANOVA table:
%      Source / ss / df / ms / F / p
%
% Notes:
%    Program does not do any input validation, so it is up to you to make
%    sure that you have passed in the parameters in the correct form:
%
%       Y, S, F1, and F2 must be numeric vectors all of the same length.
%
%       There must be at least one value in Y for each possible combination
%       of S, F1, and F2 (i.e. there must be at least one measurement per
%       subject per condition).
%
%       If there is more than one measurement per subject X condition, then
%       the program will take the mean of those measurements.
%
% Aaron Schurger (2005.02.04)
%   Derived from Keppel & Wickens (2004) "Design and Analysis" ch. 18
%

%
% Revision history...
%
% 11 December 2009 (Aaron Schurger)
% 
% Fixed error under "bracket terms"
% was: expY = sum(Y.^2);
% now: expY = sum(sum(sum(MEANS.^2)));
%

% get observed stats
f_obs = run_anova(Y,S,F1,F2);

% get n. participants
npp = numel(unique(S));

% define perm F vector
permF = nan(nperm,1);

% get permutated interaction distribution
for perm = 1 : nperm
    
    % begin permuting Y
    Yp = nan(size(Y));
    
    % cycle through each participant
    for pp = 1 : npp
                    
        % get Y avlues for this particiapnt
        vals = Y(S==pp);
        
        % reorganise data
        Yp(S==pp) = vals(randperm(numel(vals)));
    end
    
    % get F
    permF(perm,1) = run_anova(Yp,S,F1,F2);    
end

% sort F distribution
f_perm = sort(permF);

% find position of observed F in sorted permutation distribution and
% calculate p
p = sum(f_perm>=f_obs) ./ nperm;

end


function fAB = run_anova(Y,S,F1,F2)

stats = cell(4,5);

F1_lvls = unique(F1);
F2_lvls = unique(F2);
Subjs = unique(S);

a = length(F1_lvls); % # of levels in factor 1
b = length(F2_lvls); % # of levels in factor 2
n = length(Subjs); % # of subjects

INDS = cell(a,b,n); % this will hold arrays of indices
CELLS = cell(a,b,n); % this will hold the data for each subject X condition
MEANS = zeros(a,b,n); % this will hold the means for each subj X condition

% Calculate means for each subject X condition.
% Keep data in CELLS, because in future we may want to allow options for
% how to compute the means (e.g. leaving out outliers > 3stdev, etc...).
for i=1:a % F1
    for j=1:b % F2
        for k=1:n % Subjs
            INDS{i,j,k} = find(F1==F1_lvls(i) & F2==F2_lvls(j) & S==Subjs(k));
            CELLS{i,j,k} = Y(INDS{i,j,k});
            MEANS(i,j,k) = mean(CELLS{i,j,k});
        end
    end
end

% make tables (see table 18.1, p. 402)
AB = reshape(sum(MEANS,3),a,b); % across subjects
AS = reshape(sum(MEANS,2),a,n); % across factor 2
BS = reshape(sum(MEANS,1),b,n); % across factor 1

A = sum(AB,2); % sum across columns, so result is ax1 column vector
B = sum(AB,1); % sum across rows, so result is 1xb row vector
S = sum(AS,1); % sum across columns, so result is 1xs row vector
T = sum(sum(A)); % could sum either A or B or S, choice is arbitrary

% degrees of freedom
dfAB = (a-1)*(b-1);
dfABS = (a-1)*(b-1)*(n-1);

% bracket terms (expected value)
expA = sum(A.^2)./(b*n);
expB = sum(B.^2)./(a*n);
expAB = sum(sum(AB.^2))./n;
expS = sum(S.^2)./(a*b);
expAS = sum(sum(AS.^2))./b;
expBS = sum(sum(BS.^2))./a;
expY = sum(sum(sum(MEANS.^2))); %sum(Y.^2);
expT = T^2 / (a*b*n);

% sums of squares
ssAB = expAB - expA - expB + expT;
ssABS = expY - expAB - expAS - expBS + expA + expB + expS - expT;

% mean squares
msAB = ssAB / dfAB;
msABS = ssABS / dfABS;

% f statistic
fAB = msAB / msABS;

end