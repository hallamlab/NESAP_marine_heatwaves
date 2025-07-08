%% Fig.1 – Temperature Anomalies:
% This script generates Fig. 1 for the manuscript:
% 	"Marine Heatwaves Modulate Food Webs and Carbon Transport Processes"
% Author: Mariana Bernardi Bif (MBB) | Date: 06/12/25
%
% Description:
%   - Visualizes upper-ocean thermal variability at Ocean Station Papa (50°N, 145°W)
%   - Combines three datasets: RG09 Argo climatology, Armor3D Copernicus product, and BGC-Argo float data
%   - Computes and compares monthly anomalies, Brunt–Väisälä frequency, and MLD
%
% Notes:
%   - All required CSV, mat processed and NetCDF files must be in the same folder as this script
%   - MLD is calculated using a 0.2°C potential temperature threshold
%
% Data and References:
%   - RG09 Climatology: Roemmich, D. and Gilson, J. (2009). The 2004–2008 mean and annual cycle of temperature, salinity, and steric height in the global ocean from the Argo Program. Progress in Oceanography, 82(2), 81–100. https://sio-argo.ucsd.edu/RG_Climatology.html
%   - Armor3D: Guinehut et al. (2012), Ocean Sci., 8(5):845–857; Mulet et al. (2012), Deep Sea Res. Part II, 77–80, 70–81. DOI: https://doi.org/10.48670/moi-00052
%   - The spatial domain used for Armor3D is 47°–53°N, 147°–137°W, matching the BGC-Argo float coverage.
%   - Color schemes from ColorBrewer (https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps)

close all            
clear all            
clc                  

%% Define file paths
n2_path = fullfile(pwd, 'tem+sal_1m_04-23');
temp_anom_path = fullfile(pwd, 'tem_anomaly');
float_csv = 'BGC_extracted_temperatures.csv';
armor_file = 'dataset-armor-3d-rep-monthly_1749743366953.nc';

%% Initialize storage arrays
brunt_vaisala_freq = [];   % Stores N² values over time and depth
temp_anomaly_data = [];    % Stores temperature anomalies
temperature_time_series = [];  % Stores temperature at 50 m
mld_data = [];             % Stores mixed layer depth

years = 2009:2023;
for year = years
    n2_file = fullfile(n2_path, sprintf('%02dTSO1m.csv', mod(year, 100)));
    temp_anom_file = fullfile(temp_anom_path, sprintf('anom_tem_%d_1m.csv', year));

    if isfile(n2_file) && isfile(temp_anom_file)
        n2_data = readtable(n2_file);
        anom_data = readtable(temp_anom_file);

        for month = 1:12
            sal_col = sprintf('sal%02d', month);
            tem_col = sprintf('tem%02d', month);
            anom_col = sprintf('anom_%02d', month);

            sal = n2_data.(sal_col);
            tem = n2_data.(tem_col);
            depth = n2_data.depth;
            tem_anomaly = anom_data.(anom_col);

            % Compute Brunt–Väisälä frequency
            [N2, ~] = sw_bfrq(sal, tem, depth);
            if length(N2) < length(depth)
                N2 = [N2; repmat(N2(end), length(depth) - length(N2), 1)];
            end

            time_num = datenum(year, month, 15);

            % Store processed values
            brunt_vaisala_freq = [brunt_vaisala_freq; repmat(time_num, length(depth), 1), depth, N2];
            temp_anomaly_data = [temp_anomaly_data; repmat(time_num, length(depth), 1), depth, tem_anomaly];
            temperature_time_series = [temperature_time_series; time_num, tem(depth == 50)];

            % Calculate mixed layer depth
            mld = ra_mld(sal, tem, depth, 0.2);
            mld_data = [mld_data; time_num, mld];
        end
    end
end

mld_table = array2table(mld_data, 'VariableNames', {'Time', 'ra_mld'});

%% Create interpolation grids
% Time-depth mesh for plotting N² (0–100 m) and anomalies (0–1000 m)
time_vals = linspace(min(brunt_vaisala_freq(:, 1)), max(brunt_vaisala_freq(:, 1)), 300);
depth_vals1 = linspace(0, 100, 100);
depth_vals2 = linspace(0, 1000, 150);
[T1, D1] = meshgrid(time_vals, depth_vals1);
[T2, D2] = meshgrid(time_vals, depth_vals2);

% Interpolate values across time-depth grid
F_N2 = scatteredInterpolant(brunt_vaisala_freq(:, 1), brunt_vaisala_freq(:, 2), brunt_vaisala_freq(:, 3), 'linear', 'none');
F_Temp = scatteredInterpolant(temp_anomaly_data(:, 1), temp_anomaly_data(:, 2), temp_anomaly_data(:, 3), 'linear', 'none');

N2_interp = F_N2(T1, D1);
temp_anom_interp = F_Temp(T2, D2);

%% Load and process Armor3D data
z = double(ncread(armor_file, 'depth'));
z(1) = 1;  % Set surface level to 1 m

% Load 4D temperature array (lon, lat, depth, time)
to = ncread(armor_file, 'to');
armor_time = datetime(2004,2,1) + calmonths(0:size(to, 4)-1);

% Average horizontally (lat, lon) for each time point
Tavg = nan(numel(z), size(to, 4));
for t = 1:size(to, 4)
    T = squeeze(mean(mean(to(:,:,:,t),1,'omitnan'),2,'omitnan'));
    Tavg(:,t) = T;
end

% Restrict to 2009–2023
mask = (armor_time >= datetime(2009,1,1));
armor_time_trim = armor_time(mask);
Tavg_trim = Tavg(:, mask);

% Identify target depth indices (e.g., 50 m)
target_depths = [50];
armor_idx = zeros(size(target_depths));
for i = 1:length(target_depths)
    [~, armor_idx(i)] = min(abs(z - target_depths(i)));
end
armor_temp_selected = Tavg_trim(armor_idx, :);

%% Load float CSV data and compute monthly statistics
float = readtable(float_csv);
float.date = datetime(float.date);

% Filter for 50 m depth and 2009 onward
idx = float.depth_dbar == 50 & float.date >= datetime(2009,1,1);
Tsub = float(idx, :);
Tsub.month = dateshift(Tsub.date, 'start', 'month');

% Group by month and compute mean ± std
[grp, G] = findgroups(Tsub.month);
T_mean = splitapply(@mean, Tsub.temperature_C, grp);
T_std  = splitapply(@std,  Tsub.temperature_C, grp);
T_time = G;

%% Generate multi-panel figure (Panels a–f)
fig = figure;
fig.Position(4) = fig.Position(4) * 1.6;

% Panel a: Brunt–Väisälä frequency with MLD overlay
subplot(7, 1, 1);
hold on;
pcolor(T1, D1, N2_interp);
colormap(brewermap(7, 'YlGnBu'));
shading flat;
c = colorbar;
c.Label.String = 'N^2 (1 s^{-2})';
plot(mld_table.Time, mld_table.ra_mld, 'k-', 'LineWidth', 1.5);
ylabel('Depth (m)');
set(gca, 'YDir', 'reverse');
ylim([0 100]);
datetick('x', 'yyyy', 'keeplimits');
xlim([min(time_vals), datenum(datetime(2023,12,31))]);
set(gca, 'XTickLabelRotation', 45);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.01 0.01]);
text(datenum(datetime(2023,12,31)), 10, 'a', 'FontWeight', 'bold', 'FontSize', 16, 'HorizontalAlignment', 'right');
hold off;

% Panel b: Temperature anomaly (0–1000 m)
subplot(7, 1, 2);
hold on;
pcolor(T2, D2, temp_anom_interp);
colormap(gca, brewermap(15, '*RdBu'));
shading flat;
c = colorbar;
c.Label.String = 'Temp Anomaly (°C)';
ylabel('Depth (m)');
set(gca, 'YDir', 'reverse');
ylim([0 1000]);
datetick('x', 'yyyy', 'keeplimits');
xlim([min(time_vals), datenum(datetime(2023,12,31))]);
set(gca, 'XTickLabelRotation', 45);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.01 0.01]);
text(datenum(datetime(2023,12,31)), 10, 'b', 'FontWeight', 'bold', 'FontSize', 16, 'HorizontalAlignment', 'right');
hold off;

% Panel c: Temperature time series at 50 m
subplot(7, 1, 3);
hold on;
t_seg = datenum(T_time);
y_seg = T_mean;
e_seg = T_std;
valid_idx = ~isnan(y_seg);
t_vals = t_seg(valid_idx);
d_vals = y_seg(valid_idx);
s_vals = e_seg(valid_idx);
temp_ts = sortrows(temperature_time_series);
all_y = [d_vals; armor_temp_selected(1,:)'; temp_ts(:,2)];
ymin = min(all_y, [], 'omitnan') - 0.5;
ymax = max(all_y, [], 'omitnan') + 0.5;

cut_start = datenum(datetime(2016,1,1));
cut_end   = datenum(datetime(2018,6,30));
segments = {t_vals < cut_start, t_vals > cut_end};
for seg = 1:2
    mask = segments{seg};
    t_valid = t_vals(mask);
    d_valid = d_vals(mask);
    s_valid = s_vals(mask);
    fill([t_valid; flipud(t_valid)], [d_valid - s_valid; flipud(d_valid + s_valid)], [0.1 0.4 0.9], 'FaceAlpha', 0.25, 'EdgeColor','none');
    h0 = scatter(t_valid, d_valid, 60, [0 0.3 1], '^', 'filled');
end
h1 = plot(datenum(armor_time_trim), armor_temp_selected(1,:), '-', 'Color', [0.6 0 0], 'LineWidth', 2);
h2 = plot(temp_ts(:,1), temp_ts(:,2), 'k-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'MarkerFaceColor', 'k');

ylabel('Temp (°C)');
datetick('x','yyyy','keeplimits');
xlim([min(time_vals), datenum(datetime(2023,12,31))]);
ylim([ymin ymax]);
set(gca, 'XTickLabelRotation', 45);
grid on;
text(datenum(datetime(2023,12,31)), ymax, 'c', 'FontWeight','bold','FontSize',16, 'HorizontalAlignment','right');
text(min(time_vals), ymax, '50m', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
hold off;

% Panels d–f: Anomaly time series at 50, 300, and 500 m
depths = [50, 300, 500];
subplot_labels = {'d', 'e', 'f'};
for i = 1:3
    dep_indices = find(temp_anomaly_data(:, 2) == depths(i));
    dep_values = temp_anomaly_data(dep_indices, 3);
    time_vals_sub = temp_anomaly_data(dep_indices, 1);
    subplot(7, 1, i + 3);
    hold on;
    if i == 3
        legend_handle = legend([h1 h2 h0], {'Armor3D', 'Argo','BGC-Argo'}, 'Orientation', 'horizontal');
        set(legend_handle, 'Units', 'normalized', 'Position', [0.75 0.147 0.108 0.060], 'FontSize', 16);

  
    end
    plot(time_vals_sub, dep_values, 'LineWidth', 2, 'Color', 'k', 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
    yline(0, 'r--', 'LineWidth', 1);
    if i == 2
        ylabel('Temp Anomaly (°C)');
    end
    datetick('x', 'yyyy', 'keeplimits');
    xlim([min(time_vals), datenum(datetime(2023,12,31))]);
    set(gca, 'XTickLabelRotation', 45);
    grid on;
    text(datenum(datetime(2023,12,31)), max(dep_values), subplot_labels{i}, 'FontWeight','bold', 'HorizontalAlignment','right', 'VerticalAlignment','top');
    text(min(time_vals), max(dep_values), [num2str(depths(i)) 'm'], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16);
hold off;

%% Save figure
% saveas(gcf, 'Fig1.png');
% saveas(gcf, 'Fig1.fig');
% saveas(gcf, 'Fig1.tiff');
% save('Fig1.mat');

%% Compute statistics for panel c comparisons
% Normalize timestamps to month-start
float_dates = dateshift(T_time, 'start', 'month');
tso_dates = dateshift(datetime(temperature_time_series(:,1), 'ConvertFrom', 'datenum'), 'start', 'month');
armor_dates = dateshift(armor_time_trim, 'start', 'month');

% FLOAT vs Argo
[~, ia, ib] = intersect(float_dates, tso_dates);
float_vals = T_mean(ia);
tso_vals = temperature_time_series(ib, 2);
rmse_ft = sqrt(mean((float_vals - tso_vals).^2, 'omitnan'));
bias_ft = mean(float_vals - tso_vals, 'omitnan');
R_ft = corr(float_vals, tso_vals, 'rows', 'complete');

% FLOAT vs Armor3D
[~, ia, ib] = intersect(float_dates, armor_dates);
armor_vals = armor_temp_selected(1, ib)';
float_vals2 = T_mean(ia);
rmse_fa = sqrt(mean((float_vals2 - armor_vals).^2, 'omitnan'));
bias_fa = mean(float_vals2 - armor_vals, 'omitnan');
R_fa = corr(float_vals2, armor_vals, 'rows', 'complete');

% Display results
fprintf('\nStatistics – FLOAT vs TSO1m (50 m):\n');
fprintf('RMSE = %.3f°C, Bias = %.3f°C, R = %.2f\n', rmse_ft, bias_ft, R_ft);

fprintf('\nStatistics – FLOAT vs Armor3D (50 m):\n');
fprintf('RMSE = %.3f°C, Bias = %.3f°C, R = %.2f\n', rmse_fa, bias_fa, R_fa);
