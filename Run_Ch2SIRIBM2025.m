% This script runs the SIR IBM model and plots the results

% Clear workspace
close all;
clear all;
clc;

% Set up plotting
set(0, 'defaultaxesfontsize', 16)

% Define parameters
N = 1000; % Total population size
mu = 0.18; %per contact infection probability
gamma = 0.16; %recovery probability
C = 3; % Average contacts per day
MaxTime = 60; % Simulation finish time
Realizations = 100; % Number of realizations to achieve
OutbreakThreshold = 0.20*N; % Threshold to record a successful outbreak
group_size = 4; % Number of agents per group
num_external_neighbours_range = 1:6; % Range of external neighbours to simulate
% num_external_neighbours_range = 2; % Range of external neighbours to simulate
kappa = 0.6031; %fitted

% Initial state variables
ICs = struct('S', N-1, 'I', 1, 'R', 0);


%store parameters as structure (use this when using range of xs)
para = struct('mu', mu, 'gamma', gamma, 'N', N, 'n', group_size, ...
    'C', C, 'kappa', kappa,'MaxTime', MaxTime, 'Realizations', Realizations...
    , 'OutbreakThreshold', OutbreakThreshold,'exten',num_external_neighbours_range);


% Run the Synergistic IBM model without NPI
[I_counts_total, IcountNew, stored_t2, stored_x2] = simulateSIRIBM(para);%proposed model

% Save results as .mat files
% Define fixed folder name
folder_name = 'Saved_Results';

% Create the folder in the current working directory if it doesn't exist
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end
% Get current timestamp for unique filenames
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');

% Save outputs to .mat files in the fixed folder with timestamped names
save(fullfile(folder_name, ['I_counts_total_' char(timestamp) '.mat']), 'I_counts_total');
save(fullfile(folder_name, ['IcountNew_' char(timestamp) '.mat']), 'IcountNew');
save(fullfile(folder_name, ['stored_t2_' char(timestamp) '.mat']), 'stored_t2');
save(fullfile(folder_name, ['stored_x2_' char(timestamp) '.mat']), 'stored_x2');


% Run the Synergistic EBM model for the range of x without NPI
for idx = 1:length(num_external_neighbours_range)

    para.exten = idx;
    [Classes] = SynDETSIR(para,MaxTime,ICs);

    %store the infection results of the EBM for this exten
    I_countsEBM(idx, :) = Classes.I;
end


% Run the synergistic EBM SIR model  without NPIs
synpara = para;
synpara.exten = para.exten;
[SynClasses] = SynDETSIR(synpara,MaxTime,ICs);

%% Make plots
lineStyles = {'-','--',':','-.'}; % Different line styles

figure;
% Plot the Individual-Based Model results
hold on;
for e = 1:length(num_external_neighbours_range)
    h(e) = plot(0:MaxTime, I_counts_total(e, :), 'LineWidth', 4, ...
        'LineStyle', lineStyles{mod(e, length(lineStyles)) + 1}, 'DisplayName', sprintf('= %d', num_external_neighbours_range(e)));
end
hold off;
xlabel('Time (days)', 'Interpreter', 'latex');
ylabel('Mean Infected Population', 'Interpreter', 'latex');
ax = gca;
ax.ColorOrder = viridis(e);
lg = legend(h, 'Location', 'Best');
lg.Title.String = 'External Connections';
legend('show', 'Interpreter', 'latex');
title({'Infection dynamics for Varying','Average External Connections (x) -IBM'}, 'Interpreter', 'latex');
set(gcf, 'Position', [100,254,1080,446]);


figure('Position', [100, 100, 800, 600]);

% Create tiled layout with 2 rows, 1 column
tcl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot the Individual-Based Model results
nexttile;
hold on;
h = gobjects(1, length(num_external_neighbours_range)); % Initialize handle array
for e = 1:length(num_external_neighbours_range)
    h(e) = plot(0:MaxTime, I_counts_total(e, :), 'LineWidth', 4, ...
        'LineStyle', lineStyles{mod(e, length(lineStyles)) + 1}, ...
        'DisplayName', sprintf('x = %d', num_external_neighbours_range(e)));
end
hold off;

ylabel('Infections', 'FontSize', 26, 'Interpreter', 'latex');
ax = gca;
ax.ColorOrder = viridis(e);
title({'Infection dynamics for Varying', 'External Connections (x) - IBM'}, ...
    'FontSize', 26, 'Interpreter', 'latex');

% Plot the Equation-Based Model (theoretical) results
nexttile;
hold on;
for e = 1:length(num_external_neighbours_range)
    plot(0:MaxTime, I_countsEBM(e, :), 'LineWidth', 4, ...
        'LineStyle', lineStyles{mod(e, length(lineStyles)) + 1}, ...
        'HandleVisibility', 'off'); % Hide from legend to avoid duplication
end
hold off;

xlabel('Time (days)', 'FontSize', 26, 'Interpreter', 'latex');
ylabel('Infections', 'FontSize', 26, 'Interpreter', 'latex');
ax = gca;
ax.ColorOrder = viridis(e);
title({'Infection dynamics for Varying', 'External Connections (x) - EBM'}, ...
    'FontSize', 26, 'Interpreter', 'latex');

% Create single horizontal legend south of the figure
leg = legend(h, 'Orientation', 'horizontal', 'FontSize', 14, 'Interpreter', 'latex');
leg.Title.String = 'External Connections';
leg.Title.Interpreter = 'latex';
leg.Layout.Tile = 'south'; % Place legend south of all tiles

% Set overall figure title (optional, inspired by example)
title(tcl, 'Infection Dynamics for Varying External Connections (IBM and EBM)', ...
    'FontSize', 18, 'Interpreter', 'latex');

set(gcf, 'Color', 'w');

%% plot individual (95% C.I)-based versus the deterministic infection dynamics
figure
for e = 1:length(num_external_neighbours_range)
    subplot(3,2,e)
    fill(stored_t2{e}, stored_x2{e}, 'b', 'facealpha', 0.2, 'edgecolor', 'none', 'DisplayName', 'IBM- 95% C.I');
    hold on;
    plot(0:MaxTime, I_countsEBM(e, :), 'b', 'LineWidth', 3, 'DisplayName', 'EBM');
    hold off;

    if ismember(e, [1, 3, 5])
        ylabel('Number infected', 'Interpreter', 'latex');
    end

    if ismember(e, [5, 6])
        xlabel('Time (days)', 'Interpreter', 'latex');
    end
    % Single legend in subplot 2
    if e == 2
        legend('Location', 'northeast', 'FontSize', 14);
    end

    title(sprintf('$x = %d$', num_external_neighbours_range(e)), 'Interpreter', 'latex');
end
sgtitle({'IBM Vs EBM Comparison of Infection dynamics for Varying','External Connections (x)'}, 'FontSize', 26, 'Interpreter', 'latex');

set(gcf, 'Position', [100,73,1312,627]);


%% plot individual-based (mean) versus the deterministic infection dynamics
figure
for e = 1:length(num_external_neighbours_range)
    subplot(3,2,e) %make the first subplot
    %plot the IBM results
    plot(0:MaxTime, I_counts_total(e, :), 'bo', 'MarkerFaceColor','b', 'DisplayName', 'IBM'); hold on;
    % Plot the deterministic model results
    plot(0:MaxTime, I_countsEBM(e, :), 'r', 'LineWidth', 3, 'DisplayName', 'EBM');
    hold off;
    % Display ylabel only for subplots 1, 4, and 7
    if ismember(e, [1, 3, 5])
        ylabel('Number infected', 'Interpreter', 'latex');
    end

    % Single legend in subplot 2
    if e == 2
        legend('Location', 'northeast', 'FontSize', 14);
    end
    % Display xlabel only for subplots 7, 8, and 9
    if ismember(e, [5, 6])
        xlabel('Time (days)', 'Interpreter', 'latex');
    end

    title(sprintf('$x = %d$', num_external_neighbours_range(e)), 'Interpreter', 'latex');
end
% Adjust figure size and position
set(gcf, 'Position', [100,73,1312,627]);
sgtitle({'Comparison Infection dynamics for Varying','Average External Connections (x)'}, 'FontSize', 26, 'Interpreter', 'latex');


