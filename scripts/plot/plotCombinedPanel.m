function plotCombinedPanel (fn, outf)
% Plot a scatterplot of all metabolite trajectories in a single condition.
% For example, see Figure 1B.
% fn - the filename of the input spreadsheet
%    - see the StackSpec spreadsheets on GitHub for format examples 

r = 4;  % number of RMSE for minimum peak prominence = 3
t_max = 48; % timecourse length

% load spectra
T = readmatrix(fn, 'Sheet', 1, 'Range', '2:2')';
T = T(2:end);
P = readmatrix(fn, 'Sheet', 1, 'Range', 'A:A', 'NumHeaderLines', 3);
Z = readmatrix(fn, 'Sheet', 1, 'Range', 'B4')';

% normalize spectra
meanDev = zeros(size(Z, 1), 1);
for i = 1:size(Z, 1)
    col = Z(i, :);
    meanDev(i, 1) = sqrt(mean((col(5538:8477) - mean(col(5538:8477))).^2)); % sqrt(mean((col - mean(col)).^2)); 
end
Znorm = Z ./ meanDev;

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
cMap = [
    [0, 0.4470, 0.7410]; 
    [0.8500, 0.3250, 0.0980]; 
    [0.3010, 0.7450, 0.9330];
    [0.4, 0.8, 0.1];
    [0.9290, 0.6940, 0.1250]; 
    [0.8, 0.3, 0.8]; 
    [0.0, 0.0, 0.0]];
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
        abs_conc(j, i) = sum(pks(cfg_ids == j));    % ints or pks
    end
end

% subplots
f2 = figure('Position', [10 10 900 780]);
range = 0:0.5:48;
for i = 1:size(unq_nms, 1)
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
        fit_fcn = fit(T, overall.', g, 'start', [15 -0.1 15 3]);
        plot(range, fit_fcn(range), '-', 'Color', basecolor, 'linewidth', 5);
        ylabel('average peak signal');
    else
        g = fittype( @(L, k, x0, x) L./(1 + exp(-1*k*(x - x0))) );
        fit_fcn = fit(T, overall.', g, 'start', [5 0.3 12]);
        plot(range, fit_fcn(range), '-', 'Color', basecolor, 'linewidth', 5);
    end
    disp(unq_nms(i) + ", appears at " + times(1));
    disp(fit_fcn);
end
for i = 1:size(unq_nms, 1)
    num_pks = sum(cfg_ids == i);
    overall = abs_conc(i, :) ./ num_pks;
    basecolor = cMap(i, :);
    scatter(T, overall, 60, 'MarkerFaceColor', basecolor, 'MarkerEdgeColor', [0.8 0.8 0.8], 'HandleVisibility', 'off', 'Marker', 's');
end
xlim([0 48]);
ylabel('average peak signal');
xticks([0 12 24 36 48]);
xlabel('time (h)');
ylim([0 inf]);
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 42;
ax.FontWeight = 'bold';
saveas(f2, '/mnt/data/cctm/apavao/lsf/output/' + outf);

end