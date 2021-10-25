function plotAllScaled (cMap, cfg_pk, cfg_nm, outpath, stem)
% Plot 13C stacks with colored peaks.
% The time axis is normalized to the metabolic onset of cfg_pk.
% SEE Fig. 1A,C; Fig. S1
%
% Parameters:
%  - cMap is the colormap for the metabolite labels
%  - cfg_pk is the 1H chemical shift used for time axis calibration
%  - cfg_nm is the name of the chemical used for time axis calibration
%  - outpath is the path to which the output files are written
%  - stem is the filename stem of the spectra in NMRdata
%
% This script requires 1H and 13C spectra, given by <stem + "_1H.xlsx"> and
% <stem + "_13C.xlsx">, to be present in the NMRdata folder.
%
% Outputs:
%  - <stem + "_13Cn.png"> is the surface plot of the timecourse
%  - <stem + "_13Cspecs.svg"> is the final 1D spectrum with color-coded peaks
%
% Copyright 2021 Massachusetts Host-Microbiome Center
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%

r = 4;  % number of deviations for minimum peak prominence
tolerance = 0.5; % range around reference ppm
ppm_min = 0.0;   % minimum chemical shift
ppm_max = 200.0; % maximum chemical shift
t_min = 0;       % minimum time
t_max = 36;      % maximum time

% Scale time axis by glucose metabolic bounds. %
% these references numbers are the output from getTimescale run on the glucose dataset,
% using isocaproate (ppm=0.747) as the reference peak
y1 = 4.6548; % metabolism start time
y2 = 14.3857; % metabolism end time

% get metabolic time scale
[x1 x2] = getTimescale(stem + "_1H.xlsx", cfg_pk, cfg_nm, outpath + stem + "_1H_isc.png");
%ts_fcn = @(t) ((y2 - y1)/(x2 - x1)).*(t - x1) + y1;
ts_fcn = @(t) t - x1 + y1;

% load 13C spectra and scale time axis
fn13C = stem + "_13C.xlsx";
T = readmatrix(fn13C, 'Sheet', 1, 'Range', '2:2')';
T = ts_fcn(T(2:end)); % scale time axis
%T = T(2:end);
P = readmatrix(fn13C, 'Sheet', 1, 'Range', 'A:A', 'NumHeaderLines', 3);
Z = readmatrix(fn13C, 'Sheet', 1, 'Range', 'B4')';

% trim data bounds
Z = Z(T <= t_max & T >= t_min, P >= ppm_min & P <= ppm_max);
P = P(P >= ppm_min & P <= ppm_max);
T = T(T <= t_max & T >= t_min);

% normalize spectrum signal by modified signal-to-noise
meanDev = zeros(size(Z, 1), 1);
for i = 1:size(Z, 1)
    col = Z(i, :);
%    meanDev(i, 1) = sqrt(mean((col - mean(col)).^2)); % RMSE
    meanDev(i, 1) = mean(abs(col - mean(col))); % MAE
end
meanMeanDev = mean(meanDev);
Znorm = Z ./ meanMeanDev;

% handle colors
C = max(log2(abs(Znorm)), ones(size(Znorm, 1), size(Znorm, 2)));
C = C ./ max(C, [], 'all');
C2 = log2(abs(Znorm));
C2 = C2 ./ (max(C2(:, P < ppm_max & P > ppm_min), [], 'all') - 0.1);

% load reference ppm from cfg sheet
cfg_pks = readmatrix(fn13C, 'Sheet', 2, 'Range', 'A:A');
cfg_nms = readmatrix(fn13C, 'Sheet', 2, 'Range', 'B:B', 'OutputType', 'string');
unq_nms = unique(cfg_nms, 'stable');
cfg_ids = zeros(size(cfg_pks, 1), 1);
for i = 1:size(unq_nms, 1)
    choice = cfg_nms == unq_nms(i);
    cfg_ids = cfg_ids + i*(choice);
end
cpd_ids = unique(cfg_ids);

% initialize some variables for peak ID
all_ppm = [];   % vector of ppm values for all peaks
all_pks = [];   % vector of heights for all peaks
all_int = [];   % vector of integrals for all peaks
all_tim = [];   % vector timepoints for peaks
all_ind = [];   % compound id for each peak

% get peak data, iterating over time series
for i = 1:size(Znorm,1)
    Zvec = Znorm(i, :);
    sigs = zeros(size(cfg_pks, 1), 1); % peak signals at this timepoint
    ints = zeros(size(cfg_pks, 1), 1); % peak integrals at this timepoint
    incl = zeros(size(cfg_pks, 1), 1); % indicates whether each peak meets the signal
    thresh = r;
    for j = 1:size(cfg_pks,1)
        range = P > cfg_pks(j) - tolerance & P < cfg_pks(j) + tolerance;
        result = max(Zvec .* range');
        if result >= thresh
            sigs(j) = result;
            ints(j) = -1*trapz(P(range), Zvec(range));
            incl(j) = 1;
        end
    end

    incl = logical(incl);
    all_ppm = [all_ppm; cfg_pks(incl)];  % vector of ppm values for all peaks
    all_pks = [all_pks; sigs(incl)];     % vector of heights for all peaks
    all_int = [all_int; ints(incl)];
    all_tim = [all_tim; ones(sum(incl), 1)*T(i)];
    all_ind = [all_ind; cfg_ids(incl)];   % compound id for each peak
end

% PLOT STACK %
disp("Plotting...");
% 1. PLOT SURFACES %
% set figure resolution
wi = 1.9;   % width in inches
hi = 1;     % height in inches
f1 = figure('PaperUnits', 'inches', 'PaperPosition', [0 0 3.5 hi]);
f2 = figure('PaperUnits', 'inches', 'PaperPosition', [0 0 wi 0.95]);

% configure 1D panel
ax = subplot(1, 1, 1);
hold(ax, 'on');

% configure stack
set(0, 'currentfigure', f1);
ax = subplot(1, 1, 1);
hold(ax, 'on');
view(4.2, 30);
Zsurf = Znorm;

% configure colormap
nc = (4*ones(100, 3) + gray(100)) ./5;
for j = 1:size(unq_nms, 1)
    nc = [nc; (ones(100, 3) + gray(100)) .* 0.25 .* ([1 1 1] + cMap(j, :))];
end

% plotting loop, iterate over labeled peaks
for j = 1:size(cfg_pks,1)
    % select peak environs and erase in baseline plot
    rng = find(P > cfg_pks(j) - tolerance & P < cfg_pks(j) + tolerance);
    Zarea = Zsurf(:, rng);
    Zsurf(:, rng(2:size(rng, 1) -1)) = nan;
    % plot peak region in stack and panel
    p = surf(P(rng), T, Zarea, C(:, rng) + cfg_ids(j) - 0.01, 'EdgeColor', 'None');
    p.FaceAlpha = 0.15;
    p.HandleVisibility = 'off';
    set(0, 'currentfigure', f2);
    ax = subplot(1, 1, 1);
    if (cfg_pks(j) == 24.6) % special case for leucine isocap-isoval overlap
        plot(P(rng), Zarea(end, :), 'Color', "#8C8C8C", 'LineWidth', 0.5);
    else
        plot(P(rng), Zarea(end, :), 'Color', cMap(cfg_ids(j), :), 'LineWidth', 0.5);
    end
    set(0, 'currentfigure', f1);
end
% plot baselines
p = surf(P, T, Zsurf, C, 'EdgeColor', 'none');
set(0, 'currentfigure', f2);
ax = subplot(1, 1, 1);
plot(P, Zsurf(end, :), 'Color', "#8C8C8C", 'LineWidth', 0.5);
ts = [T(1) T(end)];

% 2. ADJUST PLOT PARAMETERS %
% panel
for i = [1]
    box on;
    xlim(ax, [ppm_min ppm_max]); % [0 200]
    ylim(ax, [-inf 2*max(Znorm(end, :), [], 'all')]);
    yticks(ax, []);
    ax.LineWidth = 0.5;
    ax.FontSize = 5;
    ax.XDir = 'reverse';
    set(ax, 'TickDir', 'out');
    chi = get(ax, 'Children');
    set(ax, 'Children', flipud(chi));
    set(ax, 'Layer', 'top');
end
han=axes(f2,'visible','off');

% stack
set(0, 'currentfigure', f1);
colormap(nc);
caxis([0 size(unq_nms, 1)+1]);
p.LineWidth = 0.2;
p.FaceAlpha = 1;
p.HandleVisibility = 'off';
Ax = gca;
Ax.ZAxis.Visible = 'off';
grid(Ax, 'off');
Ax.XDir = 'reverse';
Ax.Color = 'none';
Ax.LineWidth = 0.5;
Ax.FontSize = 5;
xlim([ppm_min ppm_max]);
ylim([t_min t_max]);
zlim([0 1.5*max(Znorm(end, :), [], 'all')]);
Ax.Clipping = 'off';
set(Ax, 'Layer', 'top');
yticks(t_min:12:t_max);
xticks(ppm_min:10:ppm_max);
shading(Ax, 'interp');

% 3. PLOT STEMS %
for i = 1:size(cpd_ids, 1)
    color = cMap(i, :);
    choice = all_ind == i;
    stem3(all_ppm(choice), all_tim(choice) , all_pks(choice), 'filled', 'Marker', '.', 'Color', color, 'MarkerSize', 1.5, 'LineWidth', 0.07); % max(Znorm, [], 'all')/max(all_pks)*
end

% 4. PLOT RIDGES %
for i = 1:size(cfg_pks)
    cpd = cfg_ids(i);
    choice = all_ppm == cfg_pks(i);
    if sum(choice) > 1
        rng = min(all_tim(choice), [], 'all'):0.1:max(all_tim(choice), [], 'all');
        fit_fcn = fit(all_tim(choice), all_pks(choice), 'smoothingspline');
        plot3(ones(size(rng, 2), 1)*cfg_pks(i), rng, fit_fcn(rng), 'Color', cMap(cpd, :), 'linewidth', 0.5);
    end
end

%view(4.2, 30); % 10, 30
print(f1, outpath + stem + "_13Cn", "-dpng", "-r1200");
saveas(f2, outpath + stem + "_13Cspecs.svg");
end
