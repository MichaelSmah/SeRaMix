% Parameter Estimation for Discrete SIR-type Model with Single Mean Kappa for Plots
% This script fits a discrete-time SIR model to incidence data, estimating the
% parameter kappa for different cluster sizes (n) and external connections (x).
% The fitted results represent incidence (new infections per time step).

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
% PROCESS ALL INDIVIDUAL REALIZATIONS
% --------------------------
t_data = 1:60; % Time vector for 60 incidence time steps
kappa_estimates = zeros(num_realizations, length(n_values)); % Store kappa for each n
kappa_single_estimates = zeros(num_realizations, 1); % Store single kappa for all datasets
I_data_all_rep = cell(1, length(n_values) * num_x); % Store mean incidence for plotting
I_data_ci_lower = cell(1, length(n_values) * num_x); % Store 2.5th percentile
I_data_ci_upper = cell(1, length(n_values) * num_x); % Store 97.5th percentile
x_all_rep = zeros(1, length(n_values) * num_x);
dataset_index_rep = zeros(1, length(n_values) * num_x);
n_all_rep = zeros(1, length(n_values) * num_x);
valid_idx = true(1, length(n_values) * num_x); % Track valid datasets
for r = 1:num_realizations
    I_data_all = {};
    x_all = [];
    dataset_index = [];
    n_all = [];
    % Extract data for current realization
    for d = 1:length(all_data_sets)
        current_data = all_data_sets{d};
        current_n = n_values(d);
        for i = 1:num_x
            I_data = current_data{i}(r, 1:60); % Get r-th realization, first 60 columns
            x_all(end+1) = i;
            dataset_index(end+1) = d;
            n_all(end+1) = current_n;
            % Ensure incidence is non-negative
            I_data(I_data < 0) = 0;
            % Skip if incidence is all zeros
            if all(I_data == 0)
                fprintf('Warning: Skipping n=%d, x=%d, realization %d due to all-zero incidence\n', current_n, i, r);
                if r == 1
                    idx = (d-1) * num_x + i;
                    valid_idx(idx) = false;
                end
                continue;
            end
            % Store incidence directly
            I_data_all{end+1} = I_data;
            % Store mean and CI for first realization (for plotting)
            if r == 1
                idx = (d-1) * num_x + i;
                I_data_all_rep{idx} = mean(current_data{i}(:, 1:60), 1); % Mean incidence (1x60)
                I_data_ci_lower{idx} = quantile(current_data{i}(:, 1:60), 0.025); % 2.5th percentile (1x60)
                I_data_ci_upper{idx} = quantile(current_data{i}(:, 1:60), 0.975); % 97.5th percentile (1x60)
                x_all_rep(idx) = i;
                dataset_index_rep(idx) = d;
                n_all_rep(idx) = current_n;
            end
        end
    end
    % Skip if no valid data
    if isempty(I_data_all)
        fprintf('Warning: Realization %d has no valid data\n', r);
        continue;
    end
    % --------------------------
    % OPTIMIZATION (Separate kappa for each n)
    % --------------------------
    paramGuess = 0.5 * ones(1, length(n_values)); % Initial guess for kappa_n4, kappa_n5, kappa_n10, kappa_n20
    lb = zeros(1, length(n_values));
    ub = Inf(1, length(n_values));
    optionsLSQ = optimoptions('lsqnonlin', 'Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-8);
    [bestParams_ls, ~] = lsqnonlin(@(p) SIR_discrete_residual_I_only(p, I_data_all, x_all, dataset_index, n_all, N, mu, C, gamma, true), ...
        paramGuess, lb, ub, optionsLSQ);
    kappa_estimates(r, :) = bestParams_ls;
    % --------------------------
    % OPTIMIZATION (Single kappa for all datasets)
    % --------------------------
    paramGuess_single = 0.5;
    lb_single = 0;
    ub_single = Inf;
    [bestParams_single, ~] = lsqnonlin(@(p) SIR_discrete_residual_I_only(p, I_data_all, x_all, dataset_index, n_all, N, mu, C, gamma, false), ...
        paramGuess_single, lb_single, ub_single, optionsLSQ);
    kappa_single_estimates(r) = bestParams_single;
end

% Summarize results
kappa_mean = mean(kappa_estimates);
kappa_std = std(kappa_estimates);
kappa_ci = [quantile(kappa_estimates, 0.025); quantile(kappa_estimates, 0.975)];
kappa_median = median(kappa_estimates);
kappa_single_mean = mean(kappa_single_estimates);
kappa_single_median = median(kappa_single_estimates);
kappa_single_std = std(kappa_single_estimates);
kappa_single_ci = [quantile(kappa_single_estimates, 0.025), quantile(kappa_single_estimates, 0.975)];

% Compare single mean kappa with mean of median kappas
mean_median_kappa = mean(kappa_median);
[h, p_value] = ttest(kappa_single_estimates, mean_median_kappa, 'Alpha', 0.05);

fprintf('\n--- PARAMETER ESTIMATION SUMMARY ---\n');
fprintf('Single Mean kappa (all datasets) = %.6f, Std. Dev. = %.6f, 95%% CI = [%.6f, %.6f]\n', ...
    kappa_single_mean, kappa_single_std, kappa_single_ci(1), kappa_single_ci(2));
fprintf('Mean of Median kappa (across n) = %.6f\n', mean_median_kappa);
fprintf('T-test p-value (Single Mean kappa vs Mean of Median kappa) = %.6f\n', p_value);
if p_value < 0.05
    fprintf('Result: Single Mean kappa is significantly different from Mean of Median kappa (p < 0.05)\n');
else
    fprintf('Result: No significant difference between Single Mean kappa and Mean of Median kappa (p >= 0.05)\n');
end
for d = 1:length(n_values)
    fprintf('n=%d: Mean kappa = %.6f, Std. Dev. = %.6f, Median = %.6f, 95%% CI = [%.6f, %.6f]\n', ...
        n_values(d), kappa_mean(d), kappa_std(d), kappa_median(d), kappa_ci(1,d), kappa_ci(2,d));
end

% Plot histogram of kappa estimates
figure('Position', [300, 300, 1500, 700]);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- First plot (single kappa histogram) ---
nexttile([2, 1]);
histogram(kappa_single_estimates, 20, 'FaceColor', 'b', 'EdgeColor', 'k');
hold on
plot([0.6031 0.6031], [0 180],'k')
text(0.61,177, 'Global \kappa estimate','fontsize',18)
xlabel('$\kappa$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Frequency', 'Interpreter', 'latex', 'FontSize', 18);
title(sprintf('$\\kappa$ estimated for each Dataset \n (Mean=%.4f, Median=%.4f)', kappa_single_mean,kappa_single_median), 'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'Box', 'on');

% --- Remaining plots (per-n kappa histograms) ---
for d = 1:length(n_values)
    nexttile;
    histogram(kappa_estimates(:, d), 20, 'FaceColor', 'b', 'EdgeColor', 'k');
    xlabel(sprintf('$\\kappa$'), 'Interpreter', 'latex', 'FontSize', 18);
    % ylabel('Frequency', 'Interpreter', 'latex', 'FontSize', 18);
    title(sprintf('$\\kappa$ for each data set with n=%d \n (Mean=%.4f,Median=%.4f)', n_values(d), kappa_mean(d), kappa_median(d)), 'Interpreter', 'latex', 'FontSize', 16);
    set(gca, 'Box', 'on');
end

% Overall figure title
sgtitle('Distribution of $\kappa$ estimates (global mean from joint estimation across all datasets $\kappa$=0.6031)', ...
     'Interpreter', 'latex', 'FontSize', 18);
set(gcf, 'Color', 'w');

% Simulate model with single mean kappa to compute fitted incidence
kappa_est_ls = kappa_single_mean; % Use single mean kappa for fitting
Incidence_fits = cell(length(I_data_all_rep), 1); % Store fitted incidence
Pt_values = cell(length(I_data_all_rep), 1); % Store infection probabilities
S_values = cell(length(I_data_all_rep), 1); % Store susceptible values for incidence
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
        % - S_fit(t) * Pt(t): New infections (susceptible to infected)
        % - gamma * I_fit(t): Recoveries (infected to recovered)
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

% Plot incidence results with 95% CI (4x6 layout) - Figure 2
C = colororder("reef");
figure('Position', [100, 100, 1400, 1000]);
tcl = tiledlayout(4, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
p = gobjects(3,1); % Updated to include mean data
n_vals = [4 5 10 20];
for idx = 1:length(I_data_all_rep)
    if ~valid_idx(idx)
        continue; % Skip invalid datasets
    end
    nexttile;
    % Plot 95% CI as shaded area
    t_fill = [t_data, fliplr(t_data)]; % 1x120
    ci_fill = [I_data_ci_lower{idx}, fliplr(I_data_ci_upper{idx})]; % 1x120
    p(1) = fill(t_fill, ci_fill, C(4,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '95\% CI Data');
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


% Grouped plotting by n value with only model fits (2x2, x in descending order) - Figure 4
figure('Position', [200, 200, 1200, 950]);
tcl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
C = colororder("reef");
model_handles = [];
mean_handles = []; % Added for mean data
model_labels = [];
mean_labels = []; % Added for mean data
for d = 1:length(n_values)
    nexttile;
    hold on;
    if d == 1
        f_model = plot(-10, -10, 'k-', 'LineWidth', 3);
        f_mean = plot(-10, -10, '--', 'Color', C(1,:), 'LineWidth', 3); % Dummy for mean data
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
        h_mean = plot(t_data, I_data_all_rep{idx}, '--', 'Color', C(x_val,:), 'LineWidth', 3); % Plot mean data
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

% Residual function for optimization
function residuals = SIR_discrete_residual_I_only(params, I_data_all, x_all, dataset_index, n_all, N, mu, C, gamma, use_multiple_kappa)
    kappa = params; % Single kappa or array of kappas
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
            if use_multiple_kappa
                factor = 1 - exp(-kappa(dset) * x * n / (n+1));
            else
                factor = 1 - exp(-kappa * x * n / (n+1));
            end
            ProbInf = (I_model(t)/N) * factor;
            Pt(t) = 1 - (1 - mu * ProbInf)^C;
            S_model(t+1) = S_model(t) - S_model(t) * Pt(t);
            I_model(t+1) = I_model(t) + S_model(t) * Pt(t) - gamma * I_model(t);
            % - S_model(t) * Pt(t): New infections
            % - gamma * I_model(t): Recoveries
        end
        % Compute model incidence directly as S(t) * P(t)
        Incidence_model = S_model(1:end-1) .* Pt;
        all_residuals = [all_residuals, I_data - Incidence_model];
    end
    residuals = all_residuals;
end