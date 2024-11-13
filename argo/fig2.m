%% Figure 2 -  Map, mean monthly chlorophyll to carbon (Chl:C) ratios and POC stocks (small particles) over the float’s time series.
% modified 11/12/2024  by Mari Bif, marianabif@gmail.com

% To run this code, make sure bbp_mhw_processed.mat and chla_mhw_processed.mat are in the same
% directory of this .m file.

% This figure uses a dataset that was extracted from the Argo Global Data Assembly Center (USGODAE; usgodae.org/ftp/outgoing/argo/) and processed for small and large particles. 
% Only quality-controlled files (Sprof) and data flagged as “good” were considered. 
% We focused the analysis in the upper ocean between the surface and 400 m, with data every 5 to 10 days and acquired between June 2010 and October 2022, 
% with a gap between February 2016 and August 2018. The study region was constrained between 47˚-52.5˚N and 137.5-146˚W. For more details on how to estimate 
% small and large particles, chlorophyll concentrations and mixed layer depths please see: Bif et al. and references wherein.

% User may need to adjust figure position for their computer screen
% (F1.Position , line ~192)
% To reproduce the map, please install m_map, A mapping package for Matlab:
% https://www-old.eoas.ubc.ca/~rich/map.html

clc
clear all
close all

load("bbp_mhw_processed.mat"); % Processed bbp signals, separating small and large particles
load("chla_mhw_processed.mat"); % Processed chla concentrations, corrected using HPLC data.

%% Integrated bbp calculations: here we propagate bbp concentrations ~6-7m of depth to depth=0 to calculate integrals
 
% Extract month and year vectors
datematrix = datevec(datet);month = datematrix(:,2)';year = datematrix(:,1)';
Ez = 100; % Ez is around 100 m for the entire time series, little interannual variability

% Define date matrices
datematrix = datevec(datet); 
month = datematrix(:, 2)';
year = datematrix(:, 1)';

% Preallocate variables
bbp_big_ez = zeros(1, size(press1, 2));
bbp_small_ez = zeros(1, size(press1, 2));
bbp_big_300 = zeros(1, size(press1, 2));
bbp_small_300 = zeros(1, size(press1, 2));

% Integrated bbp above Ez and between mld and Ez
for j = 1:size(press1, 2)
    int5 = find(press1(:, j) >= Ez, 1, 'first');
    int6 = find(press1(:, j) >= 300, 1, 'first');
    
    bbp_big_ez(j) = trapz(rmmissing(bbp2(1:int5, j)));
    bbp_small_ez(j) = trapz(rmmissing(bbp3(1:int5, j)));
    
    bbp_big_300(j) = trapz(rmmissing(bbp2(int5:int6, j)));
    bbp_small_300(j) = trapz(rmmissing(bbp3(int5:int6, j)));
end

% Calculate monthly and yearly means
unique_years = unique(year);
unique_months = unique(month);

[month_mean, year_mean, mld_mean, bbp_ezbig_mean, bbp_ezsmall_mean, ...
 bbp_300big_mean, bbp_300small_mean] = deal([]);

for uniqueyear = unique_years
    for uniquemonth = unique_months
        val = find(year == uniqueyear & month == uniquemonth);
        
        if ~isempty(val)
            month_mean(end+1) = uniquemonth;
            year_mean(end+1) = uniqueyear;
            mld_mean(end+1) = nanmean(mld(val));
            
            % Mean and std calculations for each depth horizon
            bbp_ezbig_mean(end+1) = nanmean(bbp_big_ez(val));
            bbp_ezsmall_mean(end+1) = nanmean(bbp_small_ez(val));
            bbp_300big_mean(end+1) = nanmean(bbp_big_300(val));
            bbp_300small_mean(end+1) = nanmean(bbp_small_300(val));
        end
    end
end

% Remove zero columns
remove_idx = all(~month_mean, 1);
bbp_ezbig_mean(:, remove_idx) = [];
bbp_ezsmall_mean(:, remove_idx) = [];
bbp_300big_mean(:, remove_idx) = [];
bbp_300small_mean(:, remove_idx) = [];

year_mean(:, remove_idx) = []; 
mld_mean(:, remove_idx) = []; 
month_mean(:, remove_idx) = [];

% Date conversions and POC calculations
date_bands = datenum(year_mean, month_mean, 1);
tick5 = date_bands;

t = datetime(min(datet) - 6, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM') + calmonths(1:65);
t = datenum(t);

% Conversion factors for bbp to POC
c1 = 31200; c2 = 3.04;

bbp_300small_mean = c1 * bbp_300small_mean + c2;
bbp_ezsmall_mean = c1 * bbp_ezsmall_mean + c2;
bbp_300big_mean = c1 * bbp_300big_mean + c2;
bbp_ezbig_mean = c1 * bbp_ezbig_mean + c2;


%% Unique years and months
unyr = unique(year_mean); 
unmo = unique(month_mean);
num_months = numel(unmo);
num_years = numel(unyr);

% Initialize matrices with NaNs
matrix_300small = nan(num_months, num_years);
matrix_300big = nan(num_months, num_years);
matrix_ezsmall = nan(num_months, num_years);
matrix_ezbig = nan(num_months, num_years);
matrix_month = nan(num_months, num_years);

% Fill matrices by matching months and years
for i = 1:num_years
    n = find(year_mean == unyr(i));
    m = month_mean(n);

    matrix_ezsmall(m, i) = bbp_ezsmall_mean(n);
    matrix_ezbig(m, i) = bbp_ezbig_mean(n);
    matrix_300small(m, i) = bbp_300small_mean(n);
    matrix_300big(m, i) = bbp_300big_mean(n);
    matrix_month(m, i) = m;
end

%% Estimating chl:c ratios

% Integrated chl above Ez
int5 = arrayfun(@(j) find(press1(:,j) >= Ez, 1, 'first'), 1:size(press1, 2));
chla_small_ez = arrayfun(@(j) trapz(rmmissing(chla_small(1:int5(j), j))), 1:size(press1, 2));
% Initialize output variables
num_samples = length(chla_small_ez);
chla_small_ez_mean = nan(1, num_samples);
chla_small_ezstd = nan(1, num_samples);
month_mean = nan(1, num_samples);
year_mean = nan(1, num_samples);
mld_mean = nan(1, num_samples);

% Loop over years and months
for y = min(unique(year)):max(unique(year))
    for m = min(unique(month)):max(unique(month))
        val = find(year == y & month == m);
        if ~isempty(val)
            month_mean(val) = m;
            year_mean(val) = y;
            mld_mean(val) = nanmean(mld(val));
            chla_small_ez_mean(val) = nanmean(chla_small_ez(val));
            chla_small_ezstd(val) = nanstd(chla_small_ez(val));
        end
    end
end

% Remove zero columns
valid_cols = ~isnan(chla_small_ez_mean);
month_mean = month_mean(valid_cols);
chla_small_ezstd = chla_small_ezstd(valid_cols);
year_mean = year_mean(valid_cols);
mld_mean = mld_mean(valid_cols);
chla_small_ez_mean = chla_small_ez_mean(valid_cols);

% Reorganizing matrix by year
unyr = unique(year_mean); unmo = unique(month_mean);
matrix_chla_small = NaN(length(unmo), length(unyr));
matrix_month = NaN(length(unmo), length(unyr));

for i = 1:length(unyr)
    for j = 1:length(unmo)
        idx = find(year_mean == unyr(i) & month_mean == unmo(j));
        matrix_chla_small(j, i) = nanmean(chla_small_ez_mean(idx));
        matrix_month(j, i) = unmo(j);
    end
end
% Chla:C per year
chl_c_ratios = matrix_chla_small ./ matrix_ezsmall;

%% Figures - Map, BBP and chl:c mean stocks per month per year 

handles=unique(floatn);
lon2 = lon-360; % Degrees W
hexcolors = ['#F194B8';'#8BABF1';'#9B8BF4';'#0073E6';'#FAAF90';'#C44601';'#350092';...
    '#1BD6C8';'#C44601';'#FAAF90';'#B9E192';'#029356'];
%colors1 = hex2rgb(hexcolors);
colors1 = hexcolors;

F1 = figure;
F1.Position = [10,10,1000,500];
%F1.Position = [321.8000 111.4000 1276 657.6000];

t = tiledlayout(F1,3,3);
nexttile(1,[2 1]);
m_proj('lambert','lat',[45 60],'long',[-155 -125]);
set(0,'DefaultAxesFontSize', 16);
hold on
    m_coast('patch', [192 192 192]/256,'edgecolor','k');
    m_grid('box','fancy','fontsize',14);
for i=1:size(floatn,2)
    if floatn(i) == 5903274
      f1=  m_scatter(lon2(i),lat(i),28,[0 0.65 0.31],'filled');
    elseif floatn(i) == 5903714
         f2=    m_scatter(lon2(i),lat(i),28,[1 0.51 0.26],'filled');
    elseif floatn(i) == 5905988
            f3=    m_scatter(lon2(i),lat(i),28,[0 0 1],'filled');
    end
end

% Float approaximate location in 2015 and 2019
m_line(-143,49.7,'marker','*','color','k','linewi',2.5,...
          'linest','none','markersize',10); m_text(-151,50,{'2019'},'fontsize',18,'rotation',20);
m_line(-139,49.7,'marker','*','color','k','linewi',2.5,...
          'linest','none','markersize',10); m_text(-137.5,49,{'2015'},'fontsize',18);
% Line P stations
m_line(-138.4,49.34,'marker','^','color','m','linewi',2.5,...
          'linest','none','markersize',10); m_text(-151,48,{'P26'},'fontsize',18,'color','m','rotation',20);
m_line(-145,50,'marker','^','color','m','linewi',2.5,...
          'linest','none','markersize',10); m_text(-137.5,51,{'P20'},'fontsize',18,'color','m');

legend([f1(1), f2(1), f3(1)], '5903274', '5903714','5905988',Location='north')
xL=xlim;yL=ylim;
text(0.7*xL(2),0.7*yL(2),'a','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize', 18)

% POC and Chl:C ratios
nexttile(2);
for i=1:6
plot(matrix_month(:,i),chl_c_ratios(:,i),'Color',colors1(i,:),'LineWidth',3); 
hold on
end
xlim([1 12]);xticks(1:12);ylabel('Chl:C')
title('Small particles - 1^s^t MHW','FontSize', 15);
xL=xlim;yL=ylim;
text(0.99*xL(2),99*yL(1),'0-100 m','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize', 16)
text(0.99*xL(2),0.8*yL(2),'b','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize', 18)
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});xtickangle(45);
legend({'2010','2011','2012','2013','2014','2015'},'NumColumns',2);
legend('Position', [0.18, 0.12, 0.14, 0.12]);
grid on
hold off

nexttile(3);
for i=8:12
plot(matrix_month(:,i),chl_c_ratios(:,i),'Color',colors1(i,:),'LineWidth',3); 
hold on
end
xlim([1 12]);xticks(1:12);ylabel('Chl:C');ylim([0 0.04]);
xL=xlim;yL=ylim;
text(0.99*xL(2),0.5*yL(1),'0-100 m','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize', 16)
text(0.99*xL(2),0.8*yL(2),'c','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize', 18)
title('Small particles - 2^n^d MHW','FontSize', 15); 
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});xtickangle(45);
legend({'2018','2019','2020','2021','2022'});
legend('Position', [0.91, 0.45, 0.07, 0.15]);
grid on
hold off


nexttile(5);
for i=1:6
plot(matrix_month(:,i),matrix_ezsmall(:,i),'Color',colors1(i,:),'LineWidth',3); 
hold on
end
xlim([1 12]);xticks(1:12);ylabel('POC (mg m^-^2)'); ylim([0 1100]);
xL=xlim;yL=ylim;
text(0.99*xL(2),0.5*yL(1),'0-100 m','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize', 16)
text(0.99*xL(2),0.8*yL(2),'d','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize', 18)
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});xtickangle(45);
grid on
hold off

nexttile(6);
for i=8:12
plot(matrix_month(:,i),matrix_ezsmall(:,i),'Color',colors1(i,:),'LineWidth',3); 
hold on
end
xlim([1 12]);xticks(1:12);ylabel('POC (mg m^-^2)')
xL=xlim;yL=ylim;
text(0.99*xL(2),0.5*yL(1),'0-100 m','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize', 16)
text(0.99*xL(2),0.8*yL(2),'e','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize', 18)
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});xtickangle(45);
grid on
hold off

nexttile(8);
for i=1:6
plot(matrix_month(:,i),matrix_300small(:,i),'Color',colors1(i,:),'LineWidth',3); 
hold on
end
xlim([1 12]);xticks(1:12);ylabel('POC (mg m^-^2)')
ylim([0 300]);
xL=xlim;yL=ylim;
text(0.99*xL(2),0.5*yL(1),'100-300 m','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize', 16)
text(0.99*xL(2),0.8*yL(2),'f','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize', 18)
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});xtickangle(45);
grid on
hold off

nexttile(9);
for i=8:12
plot(matrix_month(:,i),matrix_300small(:,i),'Color',colors1(i,:),'LineWidth',3); 
hold on
end
xlim([1 12]);xticks(1:12);ylabel('POC (mg m^-^2)');ylim([0 1600]);
xL=xlim;yL=ylim;
text(0.99*xL(2),0.5*yL(1),'100-300 m','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize', 16)
text(0.99*xL(2),0.8*yL(2),'g','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize', 18)
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});xtickangle(45);
grid on
hold off


%saveas(gcf,'Fig2','epsc');saveas(gcf,'Fig2.tif');saveas(gcf,'Fig2.svg');saveas(gcf,'Fig2.fig');