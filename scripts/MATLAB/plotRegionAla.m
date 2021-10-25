function plotRegionsAla (fn, ppm_min, ppm_max, t_min, t_max, yaw, pitch, cmap, outf)
% Plot 13C waterfall with colored peaks at Ala-a13C.
% SEE Fig. 4B,D
%
% Parameters:
%  - fn is the filename of the spectral data
%  - ppm_min, ppm_max are the chemical shift lower and upper bounds
%  - t_min, t_max are the time lower and upper bounds
%  - yaw, pitch are the camera axis rotations for the output
%  - cmap is the colormap for the peak labels
%  - outf is the filename of the output file
%
% Output:
%  - <outf> is the waterfall plot of the timecourse
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

r = 0.4; %4;  % number of RMSE for minimum peak prominence = 3
tolerance = 0.035; %0.45; % 0.5; %0.015 range around reference ppm

disp("Loading and processing data...");
% load spectra
T = readmatrix(fn, 'Sheet', 1, 'range', '2:2')'; % 2:2
T = T(2:end);
P = readmatrix(fn, 'Sheet', 1, 'range', 'A:A', 'NumHeaderLines', 3); % 3
%P = P - 0.008;
Z = readmatrix(fn, 'Sheet', 1, 'range', 'B4')'; % B4
Z = Z(T >= t_min & T <= t_max, P >= ppm_min-0.1 & P <= ppm_max+0.1);
T = T(T >= t_min & T <= t_max);
P = P(P >= ppm_min-0.1 & P <= ppm_max+0.1);

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
C2 = C2 ./ (max(C2(:, P < ppm_max & P > ppm_min), [], 'all') - 0.1);
C2 = C;

% load expected ppm for compounds
cfg_pks = readmatrix(fn, 'Sheet', 2, 'range', 'A:A');
cfg_pks = cfg_pks + 0.01; % adjust % 0.015
cfg_nms = readmatrix(fn, 'Sheet', 2, 'range', 'B:B', 'OutputType', 'string');
unq_nms = unique(cfg_nms, 'stable');
cfg_ids = zeros(size(cfg_pks, 1), 1);
for i = 1:size(unq_nms, 1)
    choice = cfg_nms == unq_nms(i);
    cfg_ids = cfg_ids + i*(choice);
end
Cpds = [cfg_pks cfg_ids];

% plot surface
wi = 1.7;     % width in inches
hi = 1;     % height in inches
dpi = 300;  % dpi resolution
f = figure('PaperUnits', 'inches', 'PaperPosition', [0 0 wi hi]);
xlim([ppm_min ppm_max]); %[0 200] %[18 27] %[50 80]
ylim([t_min t_max]);
zlim([min(Znorm(:, P < ppm_max & P > ppm_min), [], 'all') max(Znorm(:, P < ppm_max & P > ppm_min), [], 'all')]);

hold on;

% initialize some variables for peak ID
cpd_ids = unique(cfg_ids);
cMap = cmap;

disp("Identifying peaks...");
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

disp("Plotting...");
% plot surfaces
Zsurf = Znorm;
nc = (2*ones(100, 3) + gray(100)) ./3;
for j = 1:size(unq_nms, 1)
    nc = [nc; ones(95, 3) .* cMap(j, :)];
end
for j = 1:size(cfg_pks,1)
    tol2 = tolerance;
    if j > 1 % NOTE: PEAKS MUST BE IN DESCENDING ORDER IN CFG SHEET
        if cfg_pks(j-1) < (cfg_pks(j) + 2.*tolerance)
            tol2 = (cfg_pks(j-1) - cfg_pks(j))/2.;
        end
    end
    rng = find(P > cfg_pks(j) - tolerance & P < cfg_pks(j) + tol2);
    C2(:, rng) = 0.5 + cfg_ids(j);
end
p = waterfall(P, T, Zsurf, C2); %'EdgeColor', 'none'
Z2 = Znorm(T >= t_min & T <= t_max, P >= ppm_min & P <= ppm_max);
plot3([53.8; 53.8].*ones(2, numel(T)), [T T]', [-5*ones(numel(T), 1) Z2(:, 1)]', 'Marker', 'none', 'LineWidth', 0.5, 'Color', "#C2C2C2");
plot3([53.0; 53.0].*ones(2, numel(T)), [T T]', [-5*ones(numel(T), 1) Z2(:, end)]', 'Marker', 'none', 'LineWidth', 0.5, 'Color', "#C2C2C2");
colormap(nc); % gray
caxis([0 size(unq_nms, 1)+1]);
p.LineWidth = 0.5;
p.FaceAlpha = 1;%0.35;
p.HandleVisibility = 'off';
Ax = gca;
Ax.ZAxis.Visible = 'off';
grid(Ax, 'off');
Ax.XDir = 'reverse';
Ax.Color = 'none';
Ax.LineWidth = 0.5;
Ax.FontSize = 5;
%Ax.FontWeight = 'bold';
yticks(t_min:12:t_max);
yticklabels(0:12:(t_max-t_min));
xlim([53 53.8]);
xticks([53 53.2 53.4 53.6 53.8]);
xticklabels({'53.0', '53.2', '53.4', '53.6', '53.8'});
%xlabel('chemical shift (ppm)');
Ax.Layer = 'top';

shading interp;

%legend(unq_nms, 'Location', 'northwest');
view(yaw, pitch); % 10, 30
print(f, outf, "-dpng", "-r1200");
disp("Done.");
end
