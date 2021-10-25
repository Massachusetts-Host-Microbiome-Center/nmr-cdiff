function plotGroupedRxns (fn, sn, lbn, ubn, cmap, outf)
% Plot line plot with error regions, for select dFBA results.
% This function is not flexible. You must edit the code below to
%    add, remove, or recategorize reactions from tracking. All ordering
%    and assignment is done by the column indices in the input file.
% SEE Fig. 2A
%
% Parameters:
%  - fn is the filename of the dFBA data in NMRdata/dfba
%  - sn, lbn, ubn are the sheet names for the fluxes and FVA lower/upper bounds
%  - cfg_nm is the name of the chemical used for time axis calibration
%  - cmap is the colormap for the reaction labels
%  - outf is the filename of the output file
%
% Output:
%  - <outf> is the plot of flux trajectories
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

% LOAD DATA VALUES %
% load time scale (T), condition labels (L), z values (Z)
T = readmatrix(fn, 'Sheet', sn, 'range', 'A:A');
L = readmatrix(fn, 'Sheet', sn, 'range', '1:1', 'OutputType', 'string')';
Z = readmatrix(fn, 'Sheet', sn, 'range', 'B2')';

% LOAD STDEV VALUES %
% load time scale (Tl), condition labels (Ll), z values (Zl)
Tl = readmatrix(fn, 'Sheet', lbn, 'range', 'A:A');
Ll = readmatrix(fn, 'Sheet', lbn, 'range', '1:1', 'OutputType', 'string')';
Zl = readmatrix(fn, 'Sheet', lbn, 'range', 'B2')';
% load time scale (Tu), condition labels (Lu), z values (Zu)
Tu = readmatrix(fn, 'Sheet', ubn, 'range', 'A:A');
Lu = readmatrix(fn, 'Sheet', ubn, 'range', '1:1', 'OutputType', 'string')';
Zu = readmatrix(fn, 'Sheet', ubn, 'range', 'B2')';

alazu = Zu(7, :);
Zu(7, :) = Zl(7, :);
Zl(7, :) = alazu;

if ~isequal(T, Tl) | ~isequal(L, Ll) | ~isequal(T, Tu) | ~isequal(L, Lu)
    disp("Data and Stdev conditions not equivalent.");
    disp(T);
    disp(Tl);
    disp(Tu);
    disp(L);
    disp(Ll);
    disp(Lu);
    return
end

% REORDER TILES IF DESIRED
%tiles = [1 3 5 7 2 4 6 8];
tiles = 1:8;

% PLOT FLUXES
% set figure resolution
wi = 1.2; %1.1;     % width in inches
hi = 2.6; %2.25;     % height in inches
dpi = 300;  % dpi resolution
f = figure('PaperUnits', 'inches', 'PaperPosition', [0 0 wi hi]);
t = tiledlayout(8, 1, 'TileSpacing', 'none', 'Padding', 'none');

% set y-axis scale for panels
ymaxs = [5 3 2 2 1 1 1 1];

% define reaction groupings
for i = 1:size(L, 1)
    vid = 1;
    if ismember(i, [1 2 3 4])
        vid = 2;
    elseif ismember(i, [5])
        vid = 3;
    elseif ismember(i, [6])
        vid = 4;
    elseif ismember(i, [7])
        vid = 5;
    elseif ismember(i, [8])
        vid = 6;
    elseif ismember(i, [9 10])
        vid = 7;
    elseif ismember(i, [11 12 13])
        vid = 8;
    end

    ax = nexttile(tiles(vid));
    box(ax, 'on');
    hold(ax, 'on');
    ax.LineWidth = 0.5;
    ax.FontSize = 5;
    xlim([0 36]);
    xticks(0:12:36);
    xticklabels([]);
    ylim([0 ymaxs(vid)]);
    yticks([0 ymaxs(vid)]);

    color = cmap(mod(i - 1, size(cmap, 1)) + 1, :);
    % PLOT ERROR REGIONS %
    fill([T(1:end-1),T(1:end-1),T(2:end),T(2:end)]',...
         [Zl(i,1:end-1)',Zu(i,1:end-1)',Zu(i,2:end)',Zl(i,2:end)']',...
         (color + [1 1 1])/2,...
         'FaceAlpha', 0.5,...
         'LineStyle', 'none');
    % PLOT RIDGES %
    plot(T, Z(i, :),...
         'Color', color,...
         'LineWidth', 1.25,...
         'Marker', 'none');
end

print(f, outf, "-dsvg");

end
