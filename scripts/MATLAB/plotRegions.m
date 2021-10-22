function plotRegions (fn, ppm_min, ppm_max, t_min, t_max, yaw, pitch, cmap, outf)
% fn is the filename of the spreadsheet
% surface plot with colored regions for reference ppms
% ridge tracing with smoothing spline fit

r = 4;  % number of RMSE for minimum peak prominence = 3
tolerance = 0.45; % 0.5, range around reference ppm

% load spectra
T = readmatrix(fn, 'Sheet', 1, 'range', '2:2')';
T = T(2:end);
P = readmatrix(fn, 'Sheet', 1, 'range', 'A:A', 'NumHeaderLines', 3);
Z = readmatrix(fn, 'Sheet', 1, 'range', 'B4')';

% normalize spectra
meanDev = zeros(size(Z, 1), 1);
for i = 1:size(Z, 1)
    col = Z(i, :);
    meanDev(i, 1) = sqrt(mean((col(5538:8477) - mean(col(5538:8477))).^2)); % sqrt(mean((col - mean(col)).^2)); 
end
Znorm = Z ./ meanDev;
C = log2(Znorm - min(Znorm, [], 'all'));
C = C ./ max(C, [], 'all');

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
f = figure('Position', [10 10 1536 864]);
hold on

% initialize some variables for peak ID
cpd_ids = unique(cfg_ids);
cMap = cmap;
cMap2 = [
    [0, 0.4470, 0.7410]; 
    [0.8500, 0.3250, 0.0980]; 
    [0.3010, 0.7450, 0.9330];
    [0.4, 0.8, 0.1]; 
    [0.9290, 0.6940, 0.1250];
    [0.8, 0.3, 0.8]; 
    [0.0, 0.0, 0.0]];
cMap3 = [
    [0, 0.4470, 0.7410]; 
    [0.8500, 0.3250, 0.0980]; 
    [0.3010, 0.7450, 0.9330];
    [0.4, 0.8, 0.1]; 
    [0.9290, 0.6940, 0.1250];
    [0.8, 0.3, 0.8]; 
    [0.59, 0.29, 0];
    [0.98, 0.98, 0.0];   
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

% plot surfaces
Zsurf = Znorm;
nc = (ones(100, 3) + gray(100)) .* 0.5;
for j = 1:size(unq_nms, 1)
    nc = [nc; (ones(100, 3) + gray(100)) .* 0.25 .* ([1 1 1] + cMap(j, :))];
end
for j = 1:size(cfg_pks,1)
    rng = find(P > cfg_pks(j) - tolerance & P < cfg_pks(j) + tolerance);
    Zarea = Zsurf(:, rng);
    Zsurf(:, rng(2:size(rng, 1) -1)) = nan;
    p = surf(P(rng), T, Zarea, C(:, rng) + cfg_ids(j) - 0.01, 'EdgeColor', 'None');
    p.FaceAlpha = 0.15;
    p.HandleVisibility = 'off';
end    
p = surf(P, T, Zsurf, C, 'EdgeColor', 'none');
colormap(nc); % gray
caxis([0 size(unq_nms, 1)+1]);
p.LineWidth = 0.2;
p.FaceAlpha = 0.35;
p.HandleVisibility = 'off';
Ax = gca;
Ax.ZAxis.Visible = 'off';
grid(Ax, 'off');
Ax.XDir = 'reverse';
Ax.Color = 'none';
Ax.LineWidth = 3;
Ax.FontSize = 42;
Ax.FontWeight = 'bold';
xlabel('¹³C chemical shift (ppm)');
ylabel('time (h)');
xlim([ppm_min ppm_max]); %[0 200] %[18 27] %[50 80]
ylim([t_min t_max]);
yticks([0 12 24 36 48]);
title('13C-NMR stack Leucine');
shading(Ax, 'interp');

% plot stems
for i = 1:size(cpd_ids, 1)
    color = cMap(i, :);
    choice = all_ind == i;
    stem3(all_ppm(choice), all_tim(choice) , all_pks(choice), 'filled', 'Marker', 'd', 'Color', color, 'MarkerSize', 7); % max(Znorm, [], 'all')/max(all_pks)*
end

% plot ridges
for i = 1:size(cfg_pks)
    cpd = cfg_ids(i);
    choice = all_ppm == cfg_pks(i);
    if sum(choice) > 1
        fit_fcn = fit(all_tim(choice), all_pks(choice), 'smoothingspline');
        rng = min(all_tim(choice), [], 'all'):0.1:max(all_tim(choice), [], 'all');
        plot3(ones(size(rng, 2), 1)*cfg_pks(i), rng, fit_fcn(rng), 'Color', cMap(cpd, :), 'linewidth', 4);
    end
end

zlim([-3 inf]);

%legend(unq_nms, 'Location', 'northwest');
view(yaw, pitch); % 10, 30
saveas(f, '/mnt/data/cctm/apavao/lsf/output/' + outf);

end