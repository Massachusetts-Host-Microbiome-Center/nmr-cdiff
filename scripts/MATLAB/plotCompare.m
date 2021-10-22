function plotCompare (fn1, fn2)
% Plots scatter with logistic fit of of 2 conditions superimposed.
% For example, see figure 2C.
% fn1, fn2 - the filenames of the spreadsheets

r = 4;  % number of signal-to-noise for minimum peak prominence = 3
t_max = 48; % hours

cMap = [
    [0, 0.4470, 0.7410]; 
    [0.8500, 0.3250, 0.0980]; 
    [0.9290, 0.6940, 0.1250]; 
    [0.4940, 0.1840, 0.5560]; 
    [0.4660, 0.6740, 0.1880]; 
    [0.3010, 0.7450, 0.9330]; 
    [0.0, 0.0, 0.0]];

filevec = [fn1 fn2];
for fileindex = [1 2]
    fn = filevec(fileindex);
    % load spectra
    T = readmatrix(fn, 'Sheet', 1, 'Range', '2:2')';
    T = T(2:end);
    P = readmatrix(fn, 'Sheet', 1, 'Range', 'A:A', 'NumHeaderLines', 3);
    Z = readmatrix(fn, 'Sheet', 1, 'Range', 'B4')';

    meanDev = zeros(size(Z, 1), 1);
    for i = 1:size(Z, 1)
        col = Z(i, :);
        meanDev(i, 1) = sqrt(mean((col(5538:8477) - mean(col(5538:8477))).^2)); %-1*mean(col(Nv(i, :)));
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
    Cpds = [cfg_pks cfg_ids];

    % initialize some variables for peak ID
    cpd_ids = unique(cfg_ids);
    
    all_ppm = [];   % vector of ppm values for all peaks
    all_pks = [];   % vector of heights for all peaks
    all_int = [];   % vector of integrals for all peaks
    all_pks_prop = [];  % peaks normalized over time series
    all_tim = [];   % vector timepoints for peaks
    all_col = [];   % colors for each peak
    all_ind = [];   % compound id for each peak
    all_asc = [];   % associated reference peak ppms
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
        % colors = zeros(size(pks, 1), 3);
        % indices = zeros(size(pks, 1), 1);

        incl = logical(incl);
        pks_prop = pks ./ sum(pks);
        all_ppm = [all_ppm; cfg_pks(incl)];  % vector of ppm values for all peaks
        all_pks = [all_pks; pks(incl)];   % vector of heights for all peaks
        all_int = [all_int; ints(incl)];
        all_tim = [all_tim; ones(sum(incl), 1)*T(i)];
        all_ind = [all_ind; cfg_ids(incl)];   % compound id for each peak

        for j = 1:size(rel_conc, 1)
            rel_conc(j, i) = sum(pks_prop(cfg_ids == j));
            abs_conc(j, i) = sum(pks(cfg_ids == j)); % pks or ints
        end
    end
    
    if fileindex == 1
        T1 = T;
        C1 = abs_conc ./ 5;
    else
        T2 = T;
        C2 = abs_conc ./ 5;
    end
end

figure();
% total signal overlay
for i = 1:size(abs_conc, 1)
    hold on
    C1n = C1(i, :);
    C2n = C2(i, :);
    col = cMap(i, :);
    z1 = T1(find(C1n, 1));
    z2 = T2(find(C2n, 1));
    if i == 1
        vals1 = C1n.' - min(C1n, [], 'all');
        vals2 = C2n.' - min(C2n, [], 'all');
        g = fittype( @(L, k, x0, x) L./(1 + exp(-1*k*(x - x0))) );
        g2 = fittype( @(k, x0, x) 70 ./(1 + exp(-1*k*(x - x0))) );
        fit1 = fit(T1, vals1, g2, 'start', [-0.2 6]);
        fit2 = fit(T2, vals2, g, 'start', [15 -0.1 15]);
        plot(0:0.2:48, fit1(0:0.2:48) + min(C1n, [], 'all'), 'Color', col, 'linewidth', 3);
        plot(0:0.2:48, fit2(0:0.2:48) + min(C2n, [], 'all'), 'Color', col + (1-col).*0.5, 'linewidth', 3);
    else
        g = fittype( @(L, k, x0, x) L./(1 + exp(-1*k*(x - x0))) );
        fit1 = fit(T1, C1n.', g, 'start', [15 0.1 15]);
        fit2 = fit(T2, C2n.', g, 'start', [15 0.1 15]);
        plot(0:0.2:48, fit1(0:0.2:48), 'Color', col, 'linewidth', 3);
        plot(0:0.2:48, fit2(0:0.2:48), 'Color', col + (1-col).*0.5, 'linewidth', 3);
    end
    times1 = T1(C1n > 0);
    times2 = T2(C2n > 0);
    disp("Condition 1, " + unq_nms(i) + " appears at " + times1(1));
    disp(fit1);
    disp(min(C1n, [], 'all'));
    disp("Condition 2, " + unq_nms(i) + " appears at " + times2(1));
    disp(fit2);
    disp(min(C2n, [], 'all'));
    plot(T1, C1n, '.', 'Color', col, 'MarkerSize', 20, 'HandleVisibility', 'off');
    plot(T2, C2n, '.', 'Color', col + (1-col).*0.5, 'MarkerSize', 20, 'HandleVisibility', 'off');
    xlim([0 48]);
    ylim([0 120]);
    xticks([0 12 24 36 48]);
    xlabel('time (h)');
    ylabel('average peak signal');
    Ax = gca;
    Ax.LineWidth = 2;
    Ax.FontSize = 20;
    Ax.FontWeight = 'bold';
end
legend(["Proline (–Se)", "Proline (+Se)", "5-aminovalerate (–Se)", "5-aminovalerate (+Se)"]);

end