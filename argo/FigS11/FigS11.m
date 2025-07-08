
%% Figure 2 - Correlation between POC stocks estimated using variable versus fixed Ez

% modified 06/13/2024  by Mari Bif, marianabif@gmail.com

% To run this code, make sure bbp_mhw_processed.mat and chla_mhw_processed.mat are in the same
% directory of this .m file.

% Load processed files
load('bbp_mhw_processed.mat')
load('chla_mhw_processed.mat')

% Conversion constants
c1 = 31200;
c2 = 3.04;

% Preallocate
ez_depth = NaN(1, size(chla_small,2));
poc_variableEz = NaN(1, size(chla_small,2));
poc_fixedEz = NaN(1, size(chla_small,2));

% Loop through profiles
for i = 1:size(chla_small, 2)
    chl_prof = chla_small(:,i);
    bbp_prof = bbp_small(:,i);
    z = press1(:,i);

    % Clean and align
    valid = ~isnan(chl_prof) & ~isnan(bbp_prof) & ~isnan(z);
    chl_prof = chl_prof(valid);
    bbp_prof = bbp_prof(valid);
    z = z(valid);

    if isempty(chl_prof), continue, end

    [chl_max, idx_max] = max(chl_prof);
    thresh = 0.1 * chl_max;
    deeper = find(z > z(idx_max) & chl_prof <= thresh, 1);

    if isempty(deeper), continue, end

    ez = z(deeper);
    ez_depth(i) = ez;

    % Integrate bbp to Ez
    idx_ez = z <= ez;
    idx_100 = z <= 100;

    if sum(idx_ez) < 2 || sum(idx_100) < 2, continue, end

    int_bbp_ez = trapz(z(idx_ez), bbp_prof(idx_ez));
    int_bbp_100 = trapz(z(idx_100), bbp_prof(idx_100));

    poc_variableEz(i) = c1 * int_bbp_ez + c2;
    poc_fixedEz(i) = c1 * int_bbp_100 + c2;
end

% Convert MATLAB datenum to datetime
dates = datetime(datet, 'ConvertFrom', 'datenum');
valid = ~isnan(poc_variableEz) & ~isnan(poc_fixedEz);

% Create table
T = table(dates(valid)', poc_variableEz(valid)', poc_fixedEz(valid)', ...
    'VariableNames', {'Date','POC_variableEz','POC_fixed100m'});

% Monthly means
T.Year = year(T.Date);
T.Month = month(T.Date);
G = groupsummary(T, {'Year','Month'}, 'mean');

% Linear regression
x = T.POC_fixed100m;
y = T.POC_variableEz;
mdl = fitlm(x, y);
slope = mdl.Coefficients.Estimate(2);
intercept = mdl.Coefficients.Estimate(1);
r2 = mdl.Rsquared.Ordinary;

% Plot
figure;
set(gcf, 'Position', [100, 100, 500, 400]);  % [left, bottom, width, height]
set(groot, 'defaultAxesFontSize', 14);
scatter(x, y, 40, 'k', 'filled', 'MarkerFaceAlpha', 0.2); hold on;
scatter(G.mean_POC_fixed100m, G.mean_POC_variableEz, 30, 'r', 'filled');
xvals = linspace(min(x), max(x), 100);
plot(xvals, slope*xvals + intercept, 'k--', 'LineWidth', 2);
text(min(xvals)*1.05, max(y)*0.95, sprintf('y = %.2fx + %.2f\nR^2 = %.2f', slope, intercept, r2), ...
    'FontSize', 14);
xlabel('POC (0–100 m) [mg m^{-2}]'); xlim([0 6000])
ylabel('POC (0–Ez) [mg m^{-2}]'); ylim([0 6000])
legend({'All data','Monthly means','Linear fit'}, 'Location','best');
grid on;

%saveas(gcf,'FigS11','epsc');saveas(gcf,'FigS11.tif');saveas(gcf,'FigS11.svg');saveas(gcf,'FigS11.fig');