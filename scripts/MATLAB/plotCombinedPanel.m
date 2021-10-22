function plotCombinedPanel (stem, cfg_pk, cmap, outf)
% Plot a logistic fit of all metabolites in a single condition.
% SEE Fig. 1B,D,E.
%
% Parameters:
%  - stem is the filename stem of the spectra in NMRdata
%  - cfg_pk is the 1H chemical shift used for time axis calibration
%  - cmap is the colormap for the metabolite labels
%  - outf is the filename of the output file
%
% This script requires 1H and 13C spectra, given by <stem + "_1H.xlsx"> and
% <stem + "_13C.xlsx">, to be present in the NMRdata folder.
%
% Output:
%  - <outf> is the panel
%

r = 4;  % number of RMSE for minimum peak prominence = 3
t_max = 36; % timecourse length

% scale time axis by glucose metabolic bounds
y1 = 4.6548; % metabolism start time %0
y2 = 14.3857; % metabolism end time  %32
% get metabolic time bounds
[x1 x2] = getTimescale(stem + "_1H.xlsx", cfg_pk, "Isocaproate", stem + "_1H_isc.png");
%ts_fcn = @(t) ((y2 - y1)/(x2 - x1)).*(t - x1) + y1;
ts_fcn = @(t) t - x1 + y1;

% load spectra
fn = stem + "_13C.xlsx";
T = readmatrix(fn, 'Sheet', 1, 'Range', '2:2')';
%T = T(2:end);
T = ts_fcn(T(2:end));

P = readmatrix(fn, 'Sheet', 1, 'Range', 'A:A', 'NumHeaderLines', 3);
Z = readmatrix(fn, 'Sheet', 1, 'Range', 'B4')';
Z = Z(T <= 36 & T >= 0, :);
T = T(T <= 36 & T >= 0);

% normalize spectra
meanDev = zeros(size(Z, 1), 1);
for i = 1:size(Z, 1)
    col = Z(i, :);
    %meanDev(i, 1) = sqrt(mean((col(5538:8477) - mean(col(5538:8477))).^2)); % sqrt(mean((col - mean(col)).^2)); 
    meanDev(i, 1) = mean(abs(col - mean(col))); % MAE
end
%Znorm = Z ./ meanDev;
meanMeanDev = mean(meanDev);
Znorm = Z ./ meanMeanDev;

% load expected ppm for compounds
cfg_pks = readmatrix(fn, 'Sheet', 2, 'Range', 'A:A');
cfg_nms = readmatrix(fn, 'Sheet', 2, 'Range', 'B:B', 'OutputType', 'string');
unq_nms = unique(cfg_nms, 'stable');
cfg_ids = zeros(size(cfg_pks, 1), 1);
for i = 1:size(unq_nms, 1)
    choice = cfg_nms == unq_nms(i);
    cfg_ids = cfg_ids + i*(choice);
end

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
        range = P > cfg_pks(j) - 0.5 & P < cfg_pks(j) + 0.5;
        result = max(Zvec .* range');
        if result >= thresh
            pks(j) = result;
            ints(j) = -1*trapz(P(range), Zvec(range));
            incl(j) = 1;
        end
    end
    
    incl = logical(incl);
    pks_prop = pks ./ sum(pks);
    all_ppm = [all_ppm; cfg_pks(incl)];  % vector of ppm values for all peaks
    all_pks = [all_pks; pks(incl)];   % vector of heights for all peaks
    all_int = [all_int; ints(incl)];
    all_tim = [all_tim; ones(sum(incl), 1)*T(i)];
    all_ind = [all_ind; cfg_ids(incl)];   % compound id for each peak
    
    for j = 1:size(rel_conc, 1)
        rel_conc(j, i) = sum(pks_prop(cfg_ids == j));
        abs_conc(j, i) = sum(ints(cfg_ids == j));    % ints or pks
    end
end

if strcmp(stem, "20210519_13CGlc")
    ordr = [1 6 3 2 4 5 8 7 9];
else
    ordr = 1:size(unq_nms, 1);
end

%% PLOT %%
% set figure resolution
wi = 1.5;     % width in inches
hi = 0.95;     % height in inches
dpi = 300;  % dpi resolution
f2 = figure('PaperUnits', 'inches', 'PaperPosition', [0 0 wi hi]);
range = 0:0.5:48;
for i = ordr % plot order
    %yyaxis left;
    hold on
    loci = all_ind == i;
    times = all_tim(loci);
    assoc = all_ppm(loci);
    basecolor = cMap(i, :);
    num_pks = sum(cfg_ids == i);
    ylim([0 inf]);
    overall = abs_conc(i, :) ./ num_pks;
    if i == 1
        g = fittype( @(L, k, x0, c, x) (L./(1 + exp(-1*k*(x - x0)))) + c );
        if strcmp(stem, "20201230_13CPro")
            seed = [max(overall, [], 'all') -0.1 5 0.4*max(overall, [], 'all')];
        else
            seed = [max(overall, [], 'all') -0.1 5 0];
        end
        fit_fcn = fit(T, overall.', g, 'start', seed);
        plot(range, fit_fcn(range), '-', 'Color', basecolor, 'linewidth', 1);
        ymax = max(max(overall, [], 'all'), max(fit_fcn(range), [], 'all'));
    else
        g = fittype( @(L, k, x0, x) L./(1 + exp(-1*k*(x - x0))) );
        fit_fcn = fit(T, overall.', g, 'start', [5 0.3 12]);
        coeffs = coeffvalues(fit_fcn);
	coeffs
        if (coeffs(1) > 0)
            plot(range, fit_fcn(range), '-', 'Color', basecolor, 'linewidth', 1);
        end
    end
    disp(unq_nms(i) + ", appears at " + times(1));
    disp(fit_fcn);
end
for i = ordr % plot order
    num_pks = sum(cfg_ids == i);
    overall = abs_conc(i, :) ./ num_pks;
    basecolor = cMap(i, :);
    plot(T, overall,...
         'LineStyle', 'none',...
         'Marker', '.',...
         'MarkerSize', 6,...
         'MarkerEdgeColor', basecolor,...
         'HandleVisibility', 'off');
%         'MarkerEdgeColor', [0.8 0.8 0.8],...

end
xlim([0 36]);
ylim([0 ymax]);
%ylabel('Integrated signal (a.u.)');
xticks([0 12 24 36]);
%xlabel('Time (h)');
yticks([]);
ax = gca;
ax.LineWidth = 0.5;
ax.FontSize = 5;
ax.TickDir = 'out';
%ax.FontWeight = 'bold';
print(f2, outf, '-dsvg');

end