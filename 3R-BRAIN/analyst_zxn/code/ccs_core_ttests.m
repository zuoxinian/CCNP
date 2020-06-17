function [h,p,ci,stats] = ccs_core_ttests(x,y,type,alpha)
%% Extension of the built-in function - ttest
% 
% INPUT 
%   x - sample vector 1
%   y - sample vector 2
%   type - 'paired' or 'independent'
%   alpha - significance level
% 
% OUTPUT
%   same as ttest but more statistics with stats structure such as
%   stats.eftsize (cohen's d) and its 95% CI
%      d = 0.01  --> very small effect size
%      d = 0.20  --> small effect size
%      d = 0.50  --> medium effect size
%      d = 0.80  --> large effect size
%      d = 1.20  --> very large effect size
%      d = 2.00  --> huge effect size
%   stats.power (see more details in sampsizepwr)
%
% AUTHOR
%   Xi-Nian Zuo (https://zuoxinian.github.io)
%   Created in Beijing Normal University, 06/07/2020.
% REFERENCE
%   Shinichi Nakagawa and Innes C Cuthill (2007) Effect Size, Confidence Interval 
%   and Statistical Significance: A Practical Guide for Biologists. Biol Rev Camb 
%   Philos Soc. 82(4): 591-605. 
%
%% Code working from here
if nargin < 4
    alpha = 0.05;
end
nx = numel(x); ny = numel(y);
nmax = max([nx ny]); nmin = min([nx ny]);
if strcmp(type,'paired')
    if nx==ny
        [h,p,ci,stats] = ttest(x,y,'Alpha',alpha);
        %effect size
        corrXY = corr(x,y);
        tmpeff = stats.tstat*sqrt(2*(1-corrXY)/(stats.df+1));
        tmpsed = sqrt(2*(1-corrXY)/(stats.df+1) + tmpeff^2/(2*stats.df));
        tmpeff_up = tmpeff + 1.96*tmpsed;
        tmpeff_low = tmpeff - 1.96*tmpsed;
        stats.eftsize = [tmpeff_up;tmpeff;tmpeff_low];
        %power
        tmppwr = sampsizepwr('t',[0 stats.sd],mean(x-y),[],stats.df+1,'Tail','both');
        stats.power = tmppwr;
    end
end
if strcmp(type,'independent')
    [h,p,ci,stats] = ttest2(x,y,'Alpha',alpha);
    %effect size
    tmpeff = stats.tstat*sqrt((stats.df+2)/(numel(x)+numel(y)));
    tmpsed = sqrt((stats.df+2)/(numel(x)+numel(y)) + tmpeff^2/2*stats.df);
    tmpeff_up = tmpeff + 1.96*tmpsed;
    tmpeff_low = tmpeff - 1.96*tmpsed;
    stats.eftsize = [tmpeff_up;tmpeff;tmpeff_low];
    %power
    tmppwr = sampsizepwr('t2', [mean(x) std(x)], mean(y), [], nmin, ...
        'Ratio', nmax/nmin, 'Tail','both');
    stats.power = tmppwr;
end

end

