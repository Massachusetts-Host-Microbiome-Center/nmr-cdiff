function getScaledCoeffs (cMap, cfg_pk, cfg_nm, outpath, stem)
% Calculate the time-scaled logistic coefficients to be used dFBA analyses.
% The time axis is normalized to the metabolic onset and rate of cfg_pk.
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
%  - <stem + "_pks_scl.xlsx"> is an excel spreadsheet with logistic coefficients
%    and their 95% confidence interval given by the fit function
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

r = 1;  % deviations for minimum peak prominence
tolerance = 0.45; % range around reference ppm

% Scale time axis by glucose metabolic bounds. %
% these reference numbers are the output from getTimescale run on the glucose dataset,
% using isocaproate (ppm=0.747) as the reference peak
y1 = 4.6548; % cfg_pk metabolism start time
y2 = 14.3857; % cfg_pk metabolism end time

% get metabolic time scale
[x1 x2] = getTimescale(stem + "_1H.xlsx", cfg_pk, cfg_nm, outpath + stem + "_1H_isc.png");
ts_fcn = @(t) ((y2 - y1)/(x2 - x1)).*(t - x1) + y1;

% load 13C spectra and scale time axis
fn13C = stem + "_13C.xlsx";
T = readmatrix(fn13C, 'Sheet', 1, 'Range', '2:2')';
T = ts_fcn(T(2:end)); % scale time axis
P = readmatrix(fn13C, 'Sheet', 1, 'Range', 'A:A', 'NumHeaderLines', 3);
Z = readmatrix(fn13C, 'Sheet', 1, 'Range', 'B4')';

% normalize spectrum signal by modified signal-to-noise
meanDev = zeros(size(Z, 1), 1);
for i = 1:size(Z, 1)
    col = Z(i, :);
    meanDev(i, 1) = mean(abs(col - mean(col))); % MAE
end
meanMeanDev = mean(meanDev);
Znorm = Z ./ meanMeanDev;

% load reference ppm from cfg sheet
cfg_pks = readmatrix(fn13C, 'Sheet', 2, 'Range', 'A:A');
cfg_nms = readmatrix(fn13C, 'Sheet', 2, 'Range', 'B:B', 'OutputType', 'string');
unq_nms = unique(cfg_nms, 'stable');
cfg_ids = zeros(size(cfg_pks, 1), 1);
for i = 1:size(unq_nms, 1)
    choice = cfg_nms == unq_nms(i);
    cfg_ids = cfg_ids + i*(choice);
end

% initialize some variables for peak ID
all_ppm = [];   % vector of ppm values for all peaks
all_sig = [];   % vector of signal values for all peaks
all_int = [];   % vector of integrals for all peaks
all_tim = [];   % vector timepoints for peaks
all_ind = [];   % compound id for each peak

% get peak data, iterating over time series
for i = 1:size(Znorm,1)
    Zvec = Znorm(i, :);
    ints = zeros(size(cfg_pks, 1), 1); % peak integrals at this timepoint
    sigs = zeros(size(cfg_pks, 1), 1); % peak signals at this timepoint
    incl = zeros(size(cfg_pks, 1), 1); % indicates whether each peak meets the signal threshold
    thresh = r;
    for j = 1:size(cfg_pks,1)
        range = P > cfg_pks(j) - tolerance & P < cfg_pks(j) + tolerance;
        result = max(Zvec .* range');
        if result >= thresh
            ints(j) = -1*trapz(P(range), Zvec(range));
            sigs(j) = result;
            incl(j) = 1;
        end
    end

    % append prominent peaks to data vectors
    incl = logical(incl);
    all_ppm = [all_ppm; cfg_pks(incl)];  % vector of ppm values for all peaks
    all_sig = [all_sig; sigs(incl)];
    all_int = [all_int; ints(incl)];
    all_tim = [all_tim; ones(sum(incl), 1)*T(i)];
    all_ind = [all_ind; cfg_ids(incl)];   % compound id for each peak
end

% subplots
f2 = figure('Position', [10 10 1800 1800]); % output plot
range = 0:0.5:48;
num_pks = size(unique(all_ppm), 1);
M = zeros(num_pks, 5);
Q = zeros(num_pks, 5); % 95% confint lower bound
R = zeros(num_pks, 5); % 95% confint upper bound
count = 1;
t = tiledlayout('flow', 'TileSpacing', 'compact');
for i = 1:size(unq_nms, 1)
    nexttile;
    hold on
    loci = all_ind == i;
    times = all_tim(loci);
    peaks = all_int(loci);
    signl = all_sig(loci);
    assoc = all_ppm(loci);
    unq_ppm = unique(assoc);
    disp(compose("%s: %.2f", unq_nms(i), min(times(signl > 4), [], 'all')));
    basecolor = cMap(i, :);
    num_pks = sum(cfg_ids == i);
    factor = ([1 1 1] - basecolor)*(1.0/(num_pks + 1));
    for k = 1:size(unq_ppm, 1)
        incl = unq_ppm(k) == assoc;
        plot_pks = peaks(incl);
        scatter(times(incl), plot_pks, 30, 'filled', 'MarkerFaceColor', basecolor + (k-1)*factor);

        %% Calculate logistic coefficients using typical seed values for fit function. %%
        % Proline not completely consumed. Include case for (+C).
        if (strcmp(stem, "20201230_13CPro") & i == 1)
            g = fittype( @(L, k, x0, c, x) L./(1 + exp(-1*k*(x - x0))) + c);
            fit_fcn = fit(times(incl), plot_pks, g, 'start', [2.5 -0.2 12 0]);
            ymax = 1.1*max(peaks, [], 'all'); % for plot
        else
            g = fittype( @(L, k, x0, x) L./(1 + exp(-1*k*(x - x0))));
            % if input metabolite, seed with k < 1
            if (i == 1)
                fit_fcn = fit(times(incl), plot_pks, g, 'start', [12 -0.2 12]);
                ymax = 1.1*max(peaks, [], 'all');
            else
                fit_fcn = fit(times(incl), plot_pks, g, 'start', [12 0.2 25]);
                ymax = 1.1*max(all_int(all_ind ~= 1), [], 'all');
            end
        end
        plot(range, fit_fcn(range), '--', 'Color', basecolor + (k-1)*factor, 'linewidth', 2);

        % Assemble fit coefficients and confidence interval for output %
        coeffs = coeffvalues(fit_fcn);
	cfi = confint(fit_fcn);
        M(count, 1) = unq_ppm(k);
        Q(count, 1) = unq_ppm(k);
        R(count, 1) = unq_ppm(k);
        M(count, 2) = coeffs(1);
        Q(count, 2) = cfi(1, 1);
        R(count, 2) = cfi(2, 1);
        M(count, 3) = coeffs(2);
        Q(count, 3) = cfi(1, 2);
        R(count, 3) = cfi(2, 2);
        M(count, 4) = coeffs(3);
        Q(count, 4) = cfi(1, 3);
        R(count, 4) = cfi(2, 3);
        % handle case for (+C)
        if size(coeffs, 2) > 3
            M(count, 5) = coeffs(4);
            Q(count, 5) = cfi(1, 4);
            R(count, 5) = cfi(2, 4);
        else
            M(count, 5) = 0.0;
            Q(count, 5) = 0.0;
            R(count, 5) = 0.0;
        end
        count = count + 1;
    end
    yline(r, ':k');
    title(unq_nms(i));
    Ax = gca;
    Ax.LineWidth = 3;
    Ax.FontSize = 20;
    Ax.FontWeight = 'bold';
    xlim([0 48]);
    ylim([0 inf]); %ymax
    xticks(0:12:48);
end
title(t, "5/19 ¹³C-Glc, time-normalized", 'FontSize', 20, 'FontWeight', 'bold');
ylabel(t,"¹³C integrated signal", 'FontSize', 20, 'FontWeight', 'bold');
xlabel(t,"time (h)", 'FontSize', 20, 'FontWeight', 'bold');

% Assemble and write output tables %
Tb = array2table(M, 'VariableNames', {'Shift','L','k','x0','c'});
Tblb = array2table(Q, 'VariableNames', {'Shift','L','k','x0','c'});
Tbub = array2table(R, 'VariableNames', {'Shift','L','k','x0','c'});
%disp(Tb);
writetable(Tb, outpath + stem + "_pks_scl.xlsx", 'Sheet', 'coeffs');
writetable(Tblb, outpath + stem + "_pks_scl.xlsx", 'Sheet', 'coeffs_95lb');
writetable(Tbub, outpath + stem + "_pks_scl.xlsx", 'Sheet', 'coeffs_95ub');
%saveas(f2, outpath + stem + "_pks_sclz.png");  % uncomment this line to write figure
end
