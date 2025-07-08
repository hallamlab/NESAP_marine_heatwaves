%% Script to make Fig.S11 from MHW paper. This combines data from POC processed, remote sensing, 
% gridded data and cruise data for a spatiotemporal analysis. - MBB 03/28/2025

close all
clear all
clc

cd /Users/mariana/Documents/2023_projects/Manuscript_MHW/Revision/Updated_figs_matlab

% Let's start by pulling the data for each analysis. First starting with
% S11a, which is a float x float comparison. I want to compare Chl stocks
% for the upper 100 m when two floats overlapped in time. The figure is a
% regression property-property analysis with Chl_float1 x Chl_float2

% The float dataset is a .mat file with processed data
fig11a = load("chla_NEPac_processed.mat");

% Extract data
chla_small = fig11a.chla1;
press1 = fig11a.press1;
datet = fig11a.datet;
floatn = fig11a.floatn;

% Find overlapping dates within 3 days for different floats
overlap_indices = [];
for i = 1:length(datet)
    for j = i+1:length(datet)
        if abs(datet(i) - datet(j)) <= 3 && floatn(i) ~= floatn(j)
            overlap_indices = [overlap_indices; i, j];
        end
    end
end

% Extract overlapping data
float1_data = chla_small(:, overlap_indices(:, 1));
float2_data = chla_small(:, overlap_indices(:, 2));
press1_data = press1(:, overlap_indices(:, 1));
press2_data = press1(:, overlap_indices(:, 2));

% Initialize lists to store integrated values
integrated_float1 = [];
integrated_float2 = [];

% Integrate float1_data and float2_data values for rows where press1_data < 50
for col = 1:size(press1_data, 2)
    rows_to_integrate_float1 = find(press1_data(:, col) < 100);
    rows_to_integrate_float2 = find(press2_data(:, col) < 100);
    integrated_float1 = [integrated_float1; trapz(float1_data(rows_to_integrate_float1, col))];
    integrated_float2 = [integrated_float2; trapz(float2_data(rows_to_integrate_float2, col))];
end

% Get the corresponding float numbers for labeling
float_num_1 = floatn(overlap_indices(1, 1));
float_num_2 = floatn(overlap_indices(1, 2));

% Create figure with four different plots
figure;

% Scatter plot 1: Integrated Float Data with linear fit
subplot(2, 2, 1);
scatter(integrated_float1, integrated_float2, 60, 'black', 'LineWidth', 1);
xlabel(['Chl (mg m^-^3) float #', num2str(float_num_1)]);
ylabel(['Chl (mg m^-^3) float #', num2str(float_num_2)]);
grid on;
hold on;

% Add linear fit line
coefficients = polyfit(integrated_float1, integrated_float2, 1);
poly_fit = polyval(coefficients, integrated_float1);
plot(integrated_float1, poly_fit, 'Color', [0.75, 0.75, 0], 'LineWidth', 2); % Dark yellow line
hold off;

clear all

%% Fig.S11b - aNCP per float in the region

fig11b = readtable('aNCP_export_spring_summer.txt', 'Delimiter', ' ');

% Extract relevant data
years = fig11b.year;
wmoids = fig11b.WMOID;
anpc_carbon_units = fig11b.aNCP_carbon_unit;

% Get unique WMOIDs
unique_wmoids = unique(wmoids);

% Create a color map for different WMOIDs
colors = lines(length(unique_wmoids)); % Using a different color scheme

% Plot year x aNCP_carbon_unit for each WMOID with different colors
subplot(2, 2, 2);
hold on;
for i = 1:length(unique_wmoids)
    wmoid_data = fig11b(fig11b.WMOID == unique_wmoids(i), :);
    scatter(wmoid_data.year, wmoid_data.aNCP_carbon_unit, 100, 'filled', 'MarkerFaceColor', colors(i, :), 'DisplayName', ['#', num2str(unique_wmoids(i))]);
end
ylabel('NCP (mol C m^-^2yr^-^1)');
legend('show', 'Location', 'northwest');
grid on;
hold off;

clear all
%% Fig.S11c - Remote sensing versus chla

% Load the third dataset
fig11c = readtable('climate_monthly_2008-2023.xlsx');

% Extract data
stations = fig11c.Station;
months = fig11c.Month;
chlor_a = fig11c.chlor_a;

% Get unique Stations
unique_stations = unique(stations);

% Create a color map for different Stations
colors = lines(length(unique_stations)); 

% Plot Month x chlor_a for each Station with different colors
subplot(2, 2, 3);
hold on;
for i = 1:length(unique_stations)
    station_data = fig11c(fig11c.Station == unique_stations(i), :);
    plot(station_data.Month, station_data.chlor_a, 'LineWidth', 2, 'Color', colors(i, :), 'DisplayName', ['Station ', num2str(unique_stations(i))]);
end
xlabel('Month');
ylabel('Chl (mg m^3)');
legend('show', 'Location', 'northwest');
grid on;
hold off;

set(gcf, 'Position', [100, 100, 1200, 800]); % Adjust figure size

%% Fig.S11d - Cruise data for chla, only warm months

% Load the fourth dataset
fig11d = readtable('forMarianaB_PhytoComposition_avg.xlsx');

% Extract data
dates = fig11d.Date;
stations = fig11d.Station;
tchl_a = fig11d.Tchl_a;

% Convert dates to datetime format and extract month and year
dates = datetime(dates, 'InputFormat', 'MM/dd/yy');
months = month(dates);
years = year(dates);

% Filter data for months between 6 and 9
filtered_data = fig11d(months >= 6 & months <= 9, :);

% Get unique Stations
unique_stations = unique(filtered_data.Station);

% Create a color map for different Stations
colors = lines(length(unique_stations)); 

% Plot Date x Tchl_a for each Station with different colors
subplot(2, 2, 4);
hold on;
for i = 1:length(unique_stations)
    station_data = filtered_data(strcmp(filtered_data.Station, unique_stations{i}), :);
    scatter(station_data.Date, station_data.Tchl_a, 50, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i, :), 'DisplayName', ['Station ', unique_stations{i}]);
    plot(station_data.Date, station_data.Tchl_a, 'LineWidth', 2, 'Color', colors(i, :),'HandleVisibility', 'off');
end
xlabel('Date');
ylabel('Chl');
legend('show', 'Location', 'northwest');
grid on;hold off;

% Add labels to each subplot
annotation('textbox', [0.1, 0.9, 0.1, 0.1], 'String', 'a', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');
annotation('textbox', [0.55, 0.9, 0.1, 0.1], 'String', 'b', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');
annotation('textbox', [0.1, 0.45, 0.1, 0.1], 'String', 'c', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');
annotation('textbox', [0.55, 0.45, 0.1, 0.1], 'String', 'd', 'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');

saveas(gcf, 'Fig.S11.tiff');saveas(gcf, 'Fig.S11.png');savefig(gcf, 'Fig.S11.fig');