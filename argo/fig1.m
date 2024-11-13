%% Figure 1 -  Brunt Vaisala frequency and Temperature anomalies
% modified 11/12/2024  by Mari Bif, marianabif@gmail.com

% To run this code, make sure fig1_brunt.txt and fig1_temp.txt are in the same
% directory of this .m file.

% This figure uses the dataset extracted from an Argo-base gridded product (Cheng et al., 2017). We selected data from Ocean Station Papa (50°N, 145°W) 
% to analyze temporal variability over the past decade. Temperature anomalies at each depth were calculated as deviations from 
% the climatological mean during the period 2004–2022, corresponding to the period matching BGC-Argo float and shipboard data. 
% The mixed layer depth (MLD) was defined as the depth where the temperature is 0.2°C lower than that at 10 meters, 
% which generally corresponds well with the depth at which the maximum Brunt-Väisälä frequency occurs. 
% The gridded product enabled Brunt-Väisälä estimates between 2010 and 2021.

% More information can be found on the original paper: Bif et al. Marine Heatwaves Modulate Food Webs and Carbon Transport Processes

% Ref:	Cheng, L., Trenberth, K.E., Fasullo, J., Boyer, T., Abraham, J. and Zhu, J., 2017. 
% Improved estimates of ocean heat content from 1960 to 2015. Science Advances, 3(3), p.e1601545.

% To reproduce the color scheme of these figures, please download
% ColorBrewer colorschemes:
% https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps

close all
clear all
clc

%% Brunt Vaisala time series
data = readtable('fig1_brunt.txt');
time_num = datenum(data.time);  
depth = data.depth;
N2 = data.N2;
MLD = data.MLD;
mld_indices = find(depth == MLD);

fig = figure;
fig.Position(4) = fig.Position(4) * 2;  

subplot(5,1,1);
scatter(time_num, depth, 36, N2, 'filled');  
colormap(brewermap(7, 'YlGnBu')); 
c = colorbar;  
c.Label.String = 'N2';
ylabel('Depth (m)');
set(gca, 'YDir', 'reverse');
ylim([0 100]);
datetick('x', 'yyyy', 'keeplimits');  
xtickangle(45);
xlim([min(time_num) max(time_num)]);  
set(gca, 'TickDir', 'out');  
set(gca, 'TickLength', [0.015 0.01]);
caxis([min(N2)*5 max(N2)]);
hold on;
plot(time_num(mld_indices), depth(mld_indices), 'k');  
hold off;
text(max(time_num), 90, 'a', 'FontWeight', 'bold', 'HorizontalAlignment', 'right');

%% Temperature anomalies: time series
temp_data = readtable('fig1_temp.txt');
tem_anomaly = temp_data.tem_anom;  
time_num_temp = datenum(temp_data.time);  
depth_temp = temp_data.depth;

subplot(5,1,2);
scatter(time_num_temp, depth_temp, 36, tem_anomaly, 'filled');  
colormap(subplot(5,1,2), brewermap(15, '*RdBu'));  
c = colorbar;  
c.Label.String = 'Temp (˚C)';
ylabel('Depth (m)');
set(gca, 'YDir', 'reverse');
ylim([0 1000]);
set(gca, 'TickDir', 'out');  
set(gca, 'TickLength', [0.015 0.01]);
datetick('x', 'yyyy', 'keeplimits');  
xtickangle(45);
xlim([min(time_num_temp) max(time_num_temp)]);  
text(max(time_num_temp), 900, 'b', 'FontWeight', 'bold', 'HorizontalAlignment', 'right');

%% Temperature anomalies at specified depths
depths = [50, 300, 500];
subplot_labels = {'c', 'd', 'e'};

for i = 1:length(depths)
    dep_indices = find(depth_temp == depths(i));
    dep_values = tem_anomaly(dep_indices);
    
    subplot(5,1,i+2);
    plot(time_num_temp(dep_indices), dep_values, 'LineWidth', 2, 'Color', 'k', ...
        'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 3);  
    hold on;
    yline(0, 'r--', 'LineWidth', 1);  % Horizontal line at y=0
    hold off;
    if i == 2
        ylabel('Temp Anomaly (˚C)');
    end
    datetick('x', 'yyyy', 'keeplimits');
    xtickangle(45);
    xlim([min(time_num_temp) max(time_num_temp)]);  
    grid on;
    
    text(max(time_num_temp), max(dep_values), subplot_labels{i}, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');  
    text(min(time_num_temp), max(dep_values), [num2str(depths(i)) 'm'], ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14);

%saveas(gcf,'Fig1','eps'); saveas(gcf,'Fig1','tiff'); saveas(gcf,'Fig1','svg'); saveas(gcf,'Fig1.fig');  