function refractoryTestFit(cells)
%REFRACTORYTESTFIT   Test number of modes of refractory value distribution.
%   REFRACTORYTESTFIT takes the logarithm of refractory values and fits the
%   mixture of Gaussians up to 5 components. It tests the number of modes
%   using a parametric bootstrap approach. Goodness-of-fit is measured by 
%   the Kolmogorov-Smirnov test statistic for comparison between the
%   original sample and the bootstrap distribution.
%
%   Reference: Fisher NI (1993) Statistical analysis of circular data,
%   Cambridge University Press, Cambridge. pp. 100-102.
%
%   See also FITGMDIST and VMCOMPONENTS_EM.

%   Balazs Hangya, Tamas Laszlovszky
%   Laboratory of Systems Neuroscience
%   Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu

% Reset random seed
rng('shuffle');

% Get cellIDs
cb = whichcb;
global RESDIR;

% Load ACG for corresponding cells
acgFile = which(['ACG_matrices_' cells '.mat']);
load(acgFile);
allChAT = cellids;

% Load refractory periods
RefractoryL = log(Refractory);

% Model selection
max_modes = 5;
[p_compare, AIC, BIC] = deal(nan(1,max_modes));
for NumModes = 1:max_modes
    [p_compare(NumModes), AIC(NumModes), BIC(NumModes),...
        plotX(NumModes,:), plotY(NumModes,:)] = ...
        modelselect(RefractoryL,NumModes);
end

% FIGURE1 PANEL E
% Plot models
H1 = figure;
hold on;
for i = 1:NumModes
   plot(plotX(6-i,:), plotY(6-i,:), 'LineWidth', 3)
end
title('GM model fitting')
legend('5 Gauss', '4 Gauss', '3 Gauss', '2 Gauss', '1 Gauss')
saveas(H1,[RESDIR cells '\fitGMdist' cells '\GM_models_all_' cells '.fig'])
saveas(H1,[RESDIR cells '\fitGMdist' cells '\GM_models_all_' cells '.jpeg'])
close(H1);

% SUPPLEMENTARY FIGURE S1 H-I
% Plot AIC and BIC
H2 = figure;
hold on;
plot(1:5, AIC, 'Color', [0.6 0 0], 'LineWidth', 3)
plot(1:5, BIC, 'Color', [0 0 0.6], 'LineWidth', 3)
plot(1:5, AIC, 'o', 'Color', [0.6 0 0], 'MarkerSize', 12)
plot(1:5, BIC, 'o', 'Color', [0 0 0.6], 'MarkerSize', 12)
title('AIC and BIC')
xlabel('Number of Gaussian mixtures');
xticks(1:5)
legend('AIC', 'BIC')
saveas(H2,[RESDIR cells '\fitGMdist' cells '\AIC_BIC_' cells '.fig'])
saveas(H2,[RESDIR cells '\fitGMdist' cells '\AIC_BIC_' cells '.jpeg'])
close(H2);

% Plot P-values
H3 = figure;
hold on;
plot(1:5, p_compare, 'Color', [0 0.6 0], 'LineWidth', 3)
plot(1:5, p_compare, 'o', 'Color', [0 0.6 0], 'MarkerSize', 12)
title('P value for Gaussian mixture selection')
xlabel('Number of Gaussian mixtures');
xticks(1:5)
saveas(H3,[RESDIR cells '\fitGMdist' cells '\P_values_' cells '.fig'])
saveas(H3,[RESDIR cells '\fitGMdist' cells '\P_values_' cells '.jpeg'])
close(H3);

% Plotting hist
H4 = figure;
hold on;
yyaxis left % Left axis for log(Refractory) hist values
edges = [1 2.^(0:0.8:10)];
cnts = (edges(1:end-1) + edges(2:end)) / 2;
cnts = [0 cnts];
nm = histc(log(Refractory),log(cnts));
stairs(log(cnts),nm, 'LineWidth', 3, 'Color', [0.6 0.6 0.6])
ylabel('BFCN#')
y_lim = ylim;
y_lim(2) = y_lim(2)+5;
yticks([0 y_lim(2)])
ylim([0 y_lim(2)])
% Plotting GM fits
clr = [0.8 0 0; 0 0.8 0; 0 0 0.8];
yyaxis right
for i = 1:2
   plot(plotX(i,:), plotY(i,:), '-', 'LineWidth', 3, 'Color', clr(i,:))
end
yticks([])
xlabel('Refractory (ms)')
legend('Refractory', 'One mode', 'Two modes')
legend('boxoff')
xlim([0 log(400)])
xticks([0 log(9) log(90) log(400)])
xticklabels({'0', '9', '90', '400'});
ax1 = gca;
ax1.YAxis(2).Visible = 'off';   % remove right y-axis
ax1.YAxis(1).Color = [0 0 0]; % set left y-axis black
setmyplot_tamas;
axis square;
saveas(H4,[RESDIR cells '\fitGMdist' cells '\GM_models_' cells '.fig'])
saveas(H4,[RESDIR cells '\fitGMdist' cells '\GM_models_' cells '.jpeg'])
close(H4);

save([RESDIR cells '\fitGMdist' cells '\Refractory_GMfitNew_' cells '.mat'], 'AIC','BIC','p_compare',...
    'plotX', 'plotY', 'allChAT', 'Refractory', 'RefractoryL', '-v7.3');
% -------------------------------------------------------------------------
function [p_compare, AIC, BIC, plotX, plotY] = modelselect(RefractoryL,NumModes)

% Fit a mixture of 2 Gaussian distributions on log refractory 
gm = fitgmdist(RefractoryL, NumModes,'RegularizationValue',0.1);
L = gm.NegativeLogLikelihood;  % note: directly applying the formula resulted in a more reliable log-likelihood value
x = -15:0.2:15;
y = Pmodel(gm,x,NumModes);
% figure
% plot(x,y)
% close;
plotX = x;
plotY = y;

% AIC, BIC
Porig = Pmodel(gm,RefractoryL,NumModes);
Lborig = -log(prod(Porig));   % negative log-likelihood
Q = - Lborig;
sample_size = length(RefractoryL);   % sample size of original distribution
num_params = 3 * NumModes - 1;
AIC = 2 * num_params - 2 * Q;
BIC = - 2 * Q + num_params * log(sample_size);
% Q = - L;
% gm.AIC
% gm.BIC

% Parametric bootsrtap
bss = 1000;  % bootstrap sample size
[Lb, p, KS] = deal(nan(1,bss));  % log-likelihood as goodness-of-fit measure
saveDat = nan(bss,sample_size);

for iB = 1:bss   % bootstrap
    R = Rmodel(gm,sample_size,NumModes);   % random sample from from the fitted distribution
    warning('OFF','stats:gmdistribution:FailedToConverge')
    gmb = fitgmdist(R,NumModes,'RegularizationValue',0.1);   % fit model on bootstrap distribution
    warning('backtrace')
    saveDat(iB, :) = R;
    
    % Negative log-likelihood
    P = Pmodel(gmb,R,NumModes);
    Lb(iB) = -log(prod(P));
    
    % Goodness-of-fit
    [~, p(iB), KS(iB)] = kstest(R,'cdf',[x' cdf(gmb,x')]);   % KS statistic
end

% Test for number of modes = k vs. > k using KS GOF
[~, porig, KSorig] = kstest(RefractoryL,'cdf',[x' cdf(gm,x')]);
KS_compare = sum(KS>KSorig) / bss;
p_compare = sum(p<porig) / bss;

% -------------------------------------------------------------------------
function R = Rmodel(gm,sample_size,NumModes)

% Draw
r = rand(sample_size,1);   % which component
N = randn(sample_size,1);  % standard normal 

% Determine which Gaussian to sample from
r3 = repmat(r,1,NumModes);
cC = cumsum(gm.ComponentProportion);
cC2 = repmat(cC,sample_size,1);
drc = r3 - cC2;
drc2 = drc > 0;
wD = sum(drc2,2) + 1;   % rank order of selected component based on r

% Generate sample
gmS = squeeze(gm.Sigma(1,1,:));
R = N .* gmS(wD) + gm.mu(wD);   % random sample from from the fitted distribution

% -------------------------------------------------------------------------
function P = Pmodel(gm,R,NumModes)

% Likelihood values based on the GM PDF
P = 0;
for iM = 1:NumModes
    P = P + gm.ComponentProportion(iM)*normpdf(R,gm.mu(iM),gm.Sigma(1,1,iM));
end