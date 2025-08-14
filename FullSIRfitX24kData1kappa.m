% Parameter Estimation for Discrete SIR-type Model with Single Kappa
% This script fits a discrete-time SIR model to incidence data, estimating a
% single parameter kappa across all cluster sizes (n), external connections (x),
% and realizations (r). The fitted results represent incidence (new infections per time step).
% ---------------------------
% CLEAR WORKSPACE & DEFINE DATA
% ---------------------------
clear all; clc; close all;
% Load incidence data from .mat files, where each dataset contains incidence
% (new infections per time step) for different realizations
load('IcountNew_N1kn4.mat');
datasetN1kn4x1 = IcountNew{1};
datasetN1kn4x2 = IcountNew{2};
datasetN1kn4x3 = IcountNew{3};
datasetN1kn4x4 = IcountNew{4};
datasetN1kn4x5 = IcountNew{5};
datasetN1kn4x6 = IcountNew{6};
load('IcountNew_N1kn5.mat');
datasetN1kn5x1 = IcountNew{1};
datasetN1kn5x2 = IcountNew{2};
datasetN1kn5x3 = IcountNew{3};
datasetN1kn5x4 = IcountNew{4};
datasetN1kn5x5 = IcountNew{5};
datasetN1kn5x6 = IcountNew{6};
load('IcountNew_N1kn10.mat');
datasetN1kn10x1 = IcountNew{1};
datasetN1kn10x2 = IcountNew{2};
datasetN1kn10x3 = IcountNew{3};
datasetN1kn10x4 = IcountNew{4};
datasetN1kn10x5 = IcountNew{5};
datasetN1kn10x6 = IcountNew{6};
load('IcountNew_N1kn20.mat');
datasetN1kn20x1 = IcountNew{1};
datasetN1kn20x2 = IcountNew{2};
datasetN1kn20x3 = IcountNew{3};
datasetN1kn20x4 = IcountNew{4};
datasetN1kn20x5 = IcountNew{5};
datasetN1kn20x6 = IcountNew{6};
% Organize datasets into a cell array for different n values
all_data_sets = {{datasetN1kn4x1, datasetN1kn4x2, datasetN1kn4x3, datasetN1kn4x4, datasetN1kn4x5, datasetN1kn4x6}, ...
                 {datasetN1kn5x1, datasetN1kn5x2, datasetN1kn5x3, datasetN1kn5x4, datasetN1kn5x5, datasetN1kn5x6}, ...
                 {datasetN1kn10x1, datasetN1kn10x2, datasetN1kn10x3, datasetN1kn10x4, datasetN1kn10x5, datasetN1kn10x6}, ...
                 {datasetN1kn20x1, datasetN1kn20x2, datasetN1kn20x3, datasetN1kn20x4, datasetN1kn20x5, datasetN1kn20x6}};
dataset_labels = {'n=4', 'n=5', 'n=10', 'n=20'};
n_values = [4, 5, 10, 20];
x_values = 1:6; % External connections x ranges from 1 to 6
num_x = length(x_values);
num_realizations = 1000; % Number of realizations per dataset
% ---------------------------
% DEFINE FIXED MODEL PARAMETERS
% ---------------------------
N = 1000; % Total population
mu = 0.18; % Per-capita transmission probability
C = 4; % Average number of contacts per person
gamma = 0.16; % Recovery rate (proportion of infected individuals recovering per time step)
% --------------------------
% PROCESS ALL INDIVIDUAL REALIZATIONS INTO ONE BIG DATASET
% --------------------------
t_data = 1:60; % Time vector for 60 incidence time steps
I_data_all = {}; % Store all incidence series across all r, n, x
x_all = []; % Corresponding x for each series
dataset_index = []; % Corresponding dataset index (for n) for each series
n_all = []; % Corresponding n for each series
I_data_all_rep = cell(1, length(n_values) * num_x); % Store mean incidence for plotting
I_data_ci_lower = cell(1, length(n_values) * num_x); % Store 2.5th percentile
I_data_ci_upper = cell(1, length(n_values) * num_x); % Store 97.5th percentile
x_all_rep = zeros(1, length(n_values) * num_x);
dataset_index_rep = zeros(1, length(n_values) * num_x);
n_all_rep = zeros(1, length(n_values) * num_x);
valid_idx = true(1, length(n_values) * num_x); % Track valid datasets
for d = 1:length(all_data_sets)
    current_data = all_data_sets{d};
    current_n = n_values(d);
    for i = 1:num_x
        dataset = current_data{i}; % 1000 x 60 matrix for this n and x
        idx = (d-1) * num_x + i;
        % Store mean and CI for plotting
        I_data_all_rep{idx} = mean(dataset(:, 1:60), 1); % Mean incidence (1x60)
        I_data_ci_lower{idx} = quantile(dataset(:, 1:60), 0.025); % 2.5th percentile (1x60)
        I_data_ci_upper{idx} = quantile(dataset(:, 1:60), 0.975); % 97.5th percentile (1x60)
        x_all_rep(idx) = i;
        dataset_index_rep(idx) = d;
        n_all_rep(idx) = current_n;
        % Process all realizations for fitting
        for r = 1:num_realizations
            I_data = dataset(r, 1:60); % Get r-th realization
            I_data(I_data < 0) = 0; % Ensure non-negative
            if all(I_data == 0)
                fprintf('Warning: Skipping n=%d, x=%d, realization %d due to all-zero incidence\n', current_n, i, r);
                if r == 1
                    valid_idx(idx) = false;
                end
                continue;
            end
            I_data_all{end+1} = I_data;
            x_all(end+1) = i;
            dataset_index(end+1) = d;
            n_all(end+1) = current_n;
        end
    end
end
% Skip if no valid data
if isempty(I_data_all)
    error('No valid data available for fitting.');
end
% --------------------------
% OPTIMIZATION (Single kappa for all datasets across all r, n, x)
% --------------------------
paramGuess_single = 0.5;
lb_single = 0;
ub_single = Inf;
optionsLSQ = optimoptions('lsqnonlin', 'Display', 'iter', 'TolX', 1e-8, 'TolFun', 1e-8, ...
                          'MaxFunctionEvaluations', 10000, 'MaxIterations', 1000);
[best_kappa, resnorm] = lsqnonlin(@(p) SIR_discrete_residual_I_only(p, I_data_all, x_all, dataset_index, n_all, N, mu, C, gamma, false), ...
    paramGuess_single, lb_single, ub_single, optionsLSQ);
% Output the single fitted kappa
fprintf('\n--- PARAMETER ESTIMATION SUMMARY ---\n');
fprintf('Single kappa (across all datasets, n, x, and realizations) = %.6f\n', best_kappa);
fprintf('Residual norm = %.6f\n', resnorm);
% --------------------------
% SIMULATE MODEL WITH SINGLE KAPPA FOR PLOTTING
% --------------------------
kappa_est_ls = best_kappa; % Use single fitted kappa
Incidence_fits = cell(length(I_data_all_rep), 1); % Store fitted incidence
Pt_values = cell(length(I_data_all_rep), 1); % Store infection probabilities
S_values = cell(length(I_data_all_rep), 1); % Store susceptible values
initial_infected = zeros(length(I_data_all_rep), 1); % Store initial infected
for idx = 1:length(I_data_all_rep)
    if ~valid_idx(idx)
        continue; % Skip invalid datasets
    end
    x = x_all_rep(idx);
    dset = dataset_index_rep(idx);
    n = n_all_rep(idx);
    I_data = I_data_all_rep{idx}; % Mean observed incidence
    num_points = length(I_data) + 1;
    I_fit = zeros(1, num_points); % Infected compartment
    S_fit = zeros(1, num_points); % Susceptible compartment
    I_fit(1) = 1; % Initialize with 1 infected individual
    initial_infected(idx) = I_fit(1);
    S_fit(1) = N - I_fit(1); % Initial susceptible population
    Pt = zeros(1, num_points-1); % Infection probability per time step
    for t = 1:(num_points-1)
        % Calculate infection probability
        factor = 1 - exp(-kappa_est_ls * x * n / (n+1));
        ProbInf = (I_fit(t)/N) * factor;
        Pt(t) = 1 - (1 - mu * ProbInf)^C;
        % Update susceptible and infected compartments
        S_fit(t+1) = S_fit(t) - S_fit(t) * Pt(t);
        I_fit(t+1) = I_fit(t) + S_fit(t) * Pt(t) - gamma * I_fit(t);
    end
    % Store infection probabilities and susceptibles
    Pt_values{idx} = Pt;
    S_values{idx} = S_fit(1:end-1);
    % Compute incidence directly as new infections: S(t) * P(t)
    Incidence_fits{idx} = S_values{idx} .* Pt_values{idx};
end
% Calculate maximum y-value for consistent plotting
max_infected = 0;
for idx = 1:length(I_data_all_rep)
    if ~valid_idx(idx)
        continue;
    end
    max_infected = max(max_infected, max(I_data_ci_upper{idx}));
    max_infected = max(max_infected, max(Incidence_fits{idx}));
end
y_limit = [0, ceil(max_infected * 1.1)]; % Add 10% buffer
% --------------------------
% PLOT INCIDENCE RESULTS WITH 95% CI (4x6 layout) - Figure 1
% --------------------------
C = colororder("reef");
figure('Position', [100, 100, 1400, 1000]);
tcl = tiledlayout(4, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
p = gobjects(3,1);
n_vals = [4 5 10 20];
for idx = 1:length(I_data_all_rep)
    if ~valid_idx(idx)
        continue; % Skip invalid datasets
    end
    nexttile;
    % Plot 95% CI as shaded area
    t_fill = [t_data, fliplr(t_data)]; % 1x120
    ci_fill = [I_data_ci_lower{idx}, fliplr(I_data_ci_upper{idx})]; % 1x120
    p(1) = fill(t_fill, ci_fill, C(4,:), 'FaceAlpha', 0.3, 'EdgeColor', C(4,:), 'DisplayName', '95\% CI Data');
    hold on;
    % Plot fitted incidence
    p(3) = plot(t_data, Incidence_fits{idx}, '-', 'Color', C(6,:), 'LineWidth', 3, 'DisplayName', 'Fitted Incidence');
    % Plot mean incidence data
    p(2) = plot(t_data, I_data_all_rep{idx}, '--', 'Color', C(1,:), 'LineWidth', 3, 'DisplayName', 'Mean Data');
    x_val = x_all_rep(idx);
    dset = dataset_index_rep(idx);
    if idx < 7
        title(sprintf('x = %d', x_val), 'Interpreter', 'latex', 'FontSize', 18);
    end
    if idx > 18
        xlabel('Time', 'FontSize', 16, 'Interpreter', 'latex');
    end
    if mod(idx, 6) == 1
        ylabel({sprintf('{\\bf n = %d}', n_vals(dset)), '', 'Incidence'}, ...
               'FontSize', 16, 'Interpreter', 'latex');
    end
    xlim([0, max(t_data)]);
    ylim(y_limit);
    set(gca, 'Box', 'on');
    hold off;
end
leg = legend(p, {'95\% CI Data', 'Mean Data', 'Fitted Incidence'}, ...
             'Orientation', 'horizontal', ...
             'FontSize', 14, ...
             'Interpreter', 'latex');
leg.Layout.Tile = 'south';
title(tcl, sprintf('SIR Model Fit to Data with Varying Cluster Sizes n and External Connections x (Estimated $\\kappa = %.4f$)', kappa_est_ls), ...
      'FontSize', 18, 'Interpreter', 'latex');
set(gcf, 'Color', 'w');
% --------------------------
% GROUPED PLOTTING BY n VALUE WITH MODEL FITS (2x2, x in descending order) - Figure 2
% --------------------------
figure('Position', [200, 200, 1200, 950]);
tcl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
C = colororder("reef");
model_handles = [];
mean_handles = [];
model_labels = [];
mean_labels = [];
for d = 1:length(n_values)
    nexttile;
    hold on;
    if d == 1
        f_model = plot(-10, -10, 'k-', 'LineWidth', 3);
        f_mean = plot(-10, -10, '--', 'Color', C(1,:), 'LineWidth', 3);
    end
    indices = find(dataset_index_rep == d);
    % Sort indices by x in descending order (x=6 to x=1)
    [~, sort_order] = sort(x_all_rep(indices), 'descend');
    indices = indices(sort_order);
    for i = 1:length(indices)
        idx = indices(i);
        if ~valid_idx(idx)
            continue;
        end
        x_val = x_all_rep(idx);
        h_fit = plot(t_data, Incidence_fits{idx}, '-', 'Color', C(x_val,:), 'LineWidth', 3);
        h_mean = plot(t_data, I_data_all_rep{idx}, '--', 'Color', C(x_val,:), 'LineWidth', 3);
        if d == 1
            model_handles = [model_handles, h_fit];
            mean_handles = [mean_handles, h_mean];
            model_labels = [model_labels, {sprintf('x=%d', x_val)}];
            mean_labels = [mean_labels, {sprintf('x=%d', x_val)}];
        end
    end
    if ismember(d, [1,3])
        ylabel('Incidence', 'FontSize', 18, 'Interpreter', 'latex');
    end
    if d > 2
        xlabel('Time', 'FontSize', 18, 'Interpreter', 'latex');
    end
    title(sprintf('n = %d', n_values(d)), 'FontSize', 16, 'Interpreter', 'latex');
    xlim([0, max(t_data)]);
    ylim(y_limit);
    hold off;
end
% Sort legend labels to reflect x=6 to x=1
[~, sort_order] = sort(cellfun(@(x) str2double(x(3:end)), model_labels), 'descend');
model_labels = model_labels(sort_order);
mean_labels = mean_labels(sort_order);
model_handles = model_handles(sort_order);
mean_handles = mean_handles(sort_order);
all_handles = [f_model, model_handles, f_mean, mean_handles];
all_labels = [{'Fitted Incidence'}, model_labels, {'Mean Data'}, mean_labels];
lgd = legend(all_handles, all_labels, ...
             'Orientation', 'horizontal', ...
             'NumColumns', length(model_labels)+1, ...
             'FontSize', 14, ...
             'Interpreter', 'latex');
lgd.Layout.Tile = 'south';
title(tcl, sprintf('SIR Model Fits Grouped by Cluster Sizes n and Varying External Connections x (Estimated $\\kappa = %.4f$)', kappa_est_ls), ...
      'FontSize', 18, 'Interpreter', 'latex');
set(gcf, 'Color', 'w');

% GROUPED PLOTTING BY n VALUE WITH MODEL FITS AND 95% CI (2x2, x=1 to 4 in descending order) - Figure 3
figure('Position', [200, 200, 1200, 950]);
tcl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
C = colororder("reef");
model_handles = [];
ci_handles = [];
model_labels = [];
ci_labels = [];
for d = 1:length(n_values)
    nexttile;
    hold on;
    if d == 1
        f_model = plot(-10, -10, 'k-', 'LineWidth', 3);
        f_ci = fill([-10, -10, -10, -10], [-10, -10, -10, -10], C(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    indices = find(dataset_index_rep == d);
    % Filter indices for x=1 to x=4 and sort in descending order
    valid_x = x_all_rep(indices) >= 1 & x_all_rep(indices) <= 4;
    indices = indices(valid_x);
    [~, sort_order] = sort(x_all_rep(indices), 'descend');
    indices = indices(sort_order);
    for i = 1:length(indices)
        idx = indices(i);
        if ~valid_idx(idx)
            continue;
        end
        x_val = x_all_rep(idx);
        % Plot 95% CI as shaded area
        t_fill = [t_data, fliplr(t_data)];
        ci_fill = [I_data_ci_lower{idx}, fliplr(I_data_ci_upper{idx})];
        h_ci = fill(t_fill, ci_fill, C(x_val,:), 'FaceAlpha', 0.1, 'EdgeColor', C(x_val,:));
        % Plot fitted incidence
        h_fit = plot(t_data, Incidence_fits{idx}, '-', 'Color', C(x_val,:), 'LineWidth', 3);
        if d == 1
            model_handles = [model_handles, h_fit];
            ci_handles = [ci_handles, h_ci];
            model_labels = [model_labels, {sprintf('x=%d', x_val)}];
            ci_labels = [ci_labels, {sprintf('x=%d', x_val)}];
        end
    end
    if ismember(d, [1,3])
        ylabel('Incidence', 'FontSize', 18, 'Interpreter', 'latex');
    end
    if d > 2
        xlabel('Time', 'FontSize', 18, 'Interpreter', 'latex');
    end
    title(sprintf('n = %d', n_values(d)), 'FontSize', 16, 'Interpreter', 'latex');
    xlim([0, max(t_data)]);
    ylim(y_limit);
    hold off;
end
% Sort legend labels to reflect x=4 to x=1
[~, sort_order] = sort(cellfun(@(x) str2double(x(3:end)), model_labels), 'descend');
model_labels = model_labels(sort_order);
ci_labels = ci_labels(sort_order);
model_handles = model_handles(sort_order);
ci_handles = ci_handles(sort_order);
all_handles = [f_model, model_handles, f_ci, ci_handles];
all_labels = [{'Fitted Incidence'}, model_labels, {'95\% CI Data'}, ci_labels];
lgd = legend(all_handles, all_labels, ...
    'Orientation', 'horizontal', ...
    'NumColumns', length(model_labels)+1, ...
    'FontSize', 14, ...
    'Interpreter', 'latex');
lgd.Layout.Tile = 'south';
title(tcl, sprintf('SIR Model Fits Grouped by Cluster Sizes n and External Connections x=1 to 4 (Estimated $\\kappa = %.4f$)', kappa_est_ls), ...
    'FontSize', 18, 'Interpreter', 'latex');
set(gcf, 'Color', 'w');

% --------------------------
% RESIDUAL FUNCTION FOR OPTIMIZATION
% --------------------------
function residuals = SIR_discrete_residual_I_only(params, I_data_all, x_all, dataset_index, n_all, N, mu, C, gamma, use_multiple_kappa)
    kappa = params; % Single kappa
    all_residuals = [];
    for idx = 1:length(I_data_all)
        x = x_all(idx);
        dset = dataset_index(idx);
        n = n_all(idx);
        I_data = I_data_all{idx};
        num_points = length(I_data) + 1;
        I_model = zeros(1, num_points);
        S_model = zeros(1, num_points);
        I_model(1) = 1; % Initial infected
        S_model(1) = N - I_model(1); % Initial susceptible
        Pt = zeros(1, num_points-1); % Infection probabilities
        for t = 1:(num_points-1)
            factor = 1 - exp(-kappa * x * n / (n+1));
            ProbInf = (I_model(t)/N) * factor;
            Pt(t) = 1 - (1 - mu * ProbInf)^C;
            S_model(t+1) = S_model(t) - S_model(t) * Pt(t);
            I_model(t+1) = I_model(t) + S_model(t) * Pt(t) - gamma * I_model(t);
        end
        % Compute model incidence directly as S(t) * P(t)
        Incidence_model = S_model(1:end-1) .* Pt;
        all_residuals = [all_residuals, I_data - Incidence_model];
    end
    residuals = all_residuals;
end