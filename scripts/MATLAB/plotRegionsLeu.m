function plotRegionsLeu (fn, ppm_min, ppm_max, t_min, t_max, yaw, pitch, cmap, outf)
% Plot 13C stacks with colored peaks.
% The time axis is normalized to the metabolic onset of cfg_pk.
% With special considerations for plotting leucine due to peak overlap.
% SEE Fig. 1E
%
% Parameters:
%  - fn is the filename of the spectral data in NMRdata
%  - ppm_min, ppm_max are the chemical shift lower and upper bounds
%  - t_min, t_max are the time lower and upper bounds
%  - yaw, pitch are the camera axis rotations for the output
%  - cmap is the colormap for the metabolite labels
%  - outf is the filename of the output file
%
% This script requires a 13C spectrum given by <fn> in the NMRdata folder.
%
% Output:
%  - <outf> is the surface plot of the timecourse
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

r = 4;  % number of RMSE for minimum peak prominence = 3
tolerance = 0.45; % 0.5, range around reference ppm

% scale time axis by glucose metabolic bounds
y1 = 4.6548; % metabolism start time %0
y2 = 14.3857; % metabolism end time  %32
% get metabolic time bounds
[x1 x2] = getTimescale("20210322_13CLeu_1H.xlsx", 0.752, "Isocaproate", "20210322_13CLeu_1H_isc.png");
%ts_fcn = @(t) ((y2 - y1)/(x2 - x1)).*(t - x1) + y1;
ts_fcn = @(t) t - x1 + y1;

% load spectra
T = readmatrix(fn, 'Sheet', 1, 'range', '2:2')';
T = ts_fcn(T(2:end));
%T = T(2:end);
P = readmatrix(fn, 'Sheet', 1, 'range', 'A:A', 'NumHeaderLines', 3);
Z = readmatrix(fn, 'Sheet', 1, 'range', 'B4')';
Z = Z(T <= 36 & T >= 0, P >= 0 & P <= 200);
P = P(P >= 0 & P <= 200);
T = T(T <= 36 & T >= 0);

% normalize spectra
meanDev = zeros(size(Z, 1), 1);
for i = 1:size(Z, 1)
    col = Z(i, :);
%    meanDev(i, 1) = sqrt(mean((col - mean(col)).^2)); % RMSE
    meanDev(i, 1) = mean(abs(col - mean(col))); % MAE
end
meanMeanDev = mean(meanDev);
Znorm = Z ./ meanMeanDev;
C = max(log2(abs(Znorm)), ones(size(Znorm, 1), size(Znorm, 2)));
C = C ./ max(C, [], 'all');
C2 = log2(abs(Znorm));
C2 = C2 ./ (max(C2, [], 'all') - 0.1);

% load expected ppm for compounds
cfg_pks = readmatrix(fn, 'Sheet', 2, 'range', 'A:A');
cfg_nms = readmatrix(fn, 'Sheet', 2, 'range', 'B:B', 'OutputType', 'string');
unq_nms = unique(cfg_nms, 'stable');
cfg_ids = zeros(size(cfg_pks, 1), 1);
for i = 1:size(unq_nms, 1)
    choice = cfg_nms == unq_nms(i);
    cfg_ids = cfg_ids + i*(choice);
end
Cpds = [cfg_pks cfg_ids];

% plot surface
wi = 3.5;     % width in inches
hi = 1;     % height in inches
dpi = 300;  % dpi resolution
f1 = figure('PaperUnits', 'inches', 'PaperPosition', [0 0 wi hi]);
%f = figure('Position', [0 0 1536 864]);
hold on;

% initialize some variables for peak ID
cpd_ids = unique(cfg_ids);
cMap = cmap;

all_ppm = [];   % vector of ppm values for all peaks
all_pks = [];   % vector of heights for all peaks
all_int = [];   % vector of integrals for all peaks
all_tim = [];   % vector timepoints for peaks
all_ind = [];   % compound id for each peak
rel_conc = zeros(size(cpd_ids, 1), size(T, 1));  % rel conc over time
abs_conc = zeros(size(cpd_ids, 1), size(T, 1));  % abs conc over time

% identify and plot peaks
for i = 1:size(Znorm,1)
    Zvec = Znorm(i, :);
    pks = zeros(size(cfg_pks, 1), 1);
    ints = zeros(size(cfg_pks, 1), 1);
    incl = zeros(size(cfg_pks, 1), 1);
    thresh = r;%*meanDev(i);
    for j = 1:size(cfg_pks,1)
        rng = P > cfg_pks(j) - tolerance & P < cfg_pks(j) + tolerance;
        result = max(Zvec .* rng');
        if result >= thresh
            pks(j) = result;
            ints(j) = -1*trapz(P(rng), Zvec(rng));
            incl(j) = 1;
            if cfg_ids(j) == 1 & T(i) > 24
                incl(j) = 0;
            end
        end
    end
    % colors = zeros(size(pks, 1), 3);
    % indices = zeros(size(pks, 1), 1);

    incl = logical(incl);
    pks_prop = pks ./ sum(pks);
    all_ppm = [all_ppm; cfg_pks(incl)];  % vector of ppm values for all peaks
    all_pks = [all_pks; pks(incl)];   % vector of heights for all peaks
    all_int = [all_int; ints(incl)];
%    all_pks_prop = [all_pks_prop; pks_prop(incl)];
    all_tim = [all_tim; ones(sum(incl), 1)*T(i)];
    all_ind = [all_ind; cfg_ids(incl)];   % compound id for each peak

    for j = 1:size(rel_conc, 1)
        rel_conc(j, i) = sum(pks_prop(cfg_ids == j));
        abs_conc(j, i) = sum(ints(cfg_ids == j));
    end
end

% plot surfaces
Zsurf = Znorm;
nc = (4*ones(100, 3) + gray(100)) ./5;
for j = 1:size(unq_nms, 1)
    nc = [nc; (ones(100, 3) + gray(100)) .* 0.25 .* ([1 1 1] + cMap(j, :))];
end
for j = 1:size(cfg_pks,1)
    rng = find(P > cfg_pks(j) - tolerance & P < cfg_pks(j) + tolerance);
    Zarea = Zsurf(:, rng);
    Zsurf(:, rng(2:size(rng, 1) -1)) = nan;
    if j < size(cfg_pks,1)
        p = surf(P(rng), T, Zarea, C(:, rng) + cfg_ids(j) - 0.01, 'EdgeColor', 'None');
        p.FaceAlpha = 0.1;
        p.HandleVisibility = 'off';
    else
        p1 = surf(P(rng), T(T < 6), Zarea(T < 6, :), C(T < 6, rng) + 1 - 0.01, 'EdgeColor', 'None');
        p2 = surf(P(rng), T(T >= 6), Zarea(T >= 6, :), C(T >= 6, rng) + cfg_ids(j) - 0.01, 'EdgeColor', 'None');
        p1.FaceAlpha = 0.1;
        p1.HandleVisibility = 'off';
        p2.FaceAlpha = 0.1;
        p3.HandleVisibility = 'off';
    end
end
p = surf(P, T, Zsurf, C, 'EdgeColor', 'none');
colormap(nc); % gray
caxis([0 size(unq_nms, 1)+1]);
%p.LineWidth = 0.2;
p.FaceAlpha = 1;
p.HandleVisibility = 'off';
Ax = gca;
Ax.ZAxis.Visible = 'off';
grid(Ax, 'off');
Ax.XDir = 'reverse';
Ax.Color = 'none';
Ax.LineWidth = 0.5;
Ax.FontSize = 5;
%Ax.FontWeight = 'bold';
%xlabel('Chemical shift (ppm)');
%ylabel('time (h)');
xlim([0 200]); %[0 200] %[18 27] %[50 80]
ylim([0 36]);
yticks(0:12:36);
xticks(0:50:200);
%zlim([0 max(Znorm(end, :), [], 'all')]);
set(Ax, 'Layer', 'top');
%title('13C-NMR stack Leucine');
shading(Ax, 'interp');

% plot stems
for i = 1:size(cpd_ids, 1)
    color = cMap(i, :);
    choice = all_ind == i;
    if i < numel(cpd_ids)
        stem3(all_ppm(choice), all_tim(choice) , all_pks(choice), 'filled', 'Marker', '.', 'Color', color, 'MarkerSize', 1.5, 'LineWidth', 0.07); % max(Znorm, [], 'all')/max(all_pks)*
    else
        pbase = all_ppm(choice);
	tbase = all_tim(choice);
        kbase = all_pks(choice);
        stem3(pbase(tbase <= 6), tbase(tbase <= 6) , kbase(tbase <= 6), 'filled', 'Marker', 'd', 'Color', cMap(1, :), 'MarkerSize', 1.5, 'LineWidth', 0.07);
        stem3(pbase(tbase > 6), tbase(tbase > 6) , kbase(tbase > 6), 'filled', 'Marker', 'd', 'Color', color, 'MarkerSize', 1.5, 'LineWidth', 0.07);
    end
end

% plot ridges
for i = 1:size(cfg_pks)
    cpd = cfg_ids(i);
    choice = all_ppm == cfg_pks(i);
    if sum(choice) > 1
        fit_fcn = fit(all_tim(choice), all_pks(choice), 'smoothingspline');
        rng = min(all_tim(choice), [], 'all'):0.1:max(all_tim(choice), [], 'all');
        shf = ones(size(rng, 2), 1)*cfg_pks(i);
        ftf = fit_fcn(rng);
        if i < numel(cfg_pks)
            plot3(shf, rng, ftf, 'Color', cMap(cpd, :), 'linewidth', 0.5);
        else
            plot3(shf(rng <= 6), rng(rng <= 6), ftf(rng <= 6), 'Color', cMap(1, :), 'linewidth', 0.5);
            plot3(shf(rng > 6), rng(rng > 6), ftf(rng > 6), 'Color', cMap(cpd, :), 'linewidth', 0.5);
        end
    end
end

zlim([-3 inf]);

%legend(unq_nms, 'Location', 'northwest');
view(yaw, pitch); % 10, 30
print(f1, outf, "-dpng", "-r1200");

end
