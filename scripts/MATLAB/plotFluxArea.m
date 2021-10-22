function plotFluxArea (fn)
% fn is the filename of the spreadsheet
% stacked area plot, showing proportions of flux from a metabolite

T = readmatrix(fn, 'Sheet', 2, 'Range', 'A:A', 'NumHeaderLines', 2);
N = readmatrix(fn, 'Sheet', 2, 'Range', 'B2:K2', 'OutputType', 'string');
V = readmatrix(fn, 'Sheet', 2, 'Range', 'B3');

%vsum = sum(V, 2);
%A = V ./ vsum;

f = figure('Position', [10 10 900 720]);
area(T, V, 'LineWidth', 1.5);
xlim([0 48]);
%ylim([0 1]);
xticks([0 12 24 36 48]);
xlabel('time (h)');
ylabel('\fontsize{30}estimated pyruvate flux (mmol/g_{DW}/hr)');
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 36;
ax.FontWeight = 'bold';

% legend(N, 'Location', 'northeast');
newcolors = [
    0.96 0.79 0.24;
    0.66 0.84 1.00;
    0.57 0.7 1.00;
    0.433 0.500 1.00;
    0.500 0.342 1.00;
    0.400 0.000 0.90;
    0.265 0.650 0.2;
    0.8 0.8 0.8;
    0.6 0.6 0.6;
    0.4 0.4 0.4
    ];
colororder(newcolors);
saveas(f, '/mnt/data/cctm/apavao/lsf/output/pyrflux.png');
end
