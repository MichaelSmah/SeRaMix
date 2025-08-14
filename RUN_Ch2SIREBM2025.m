% Main Run File for Semi-Random Mixing (SeRaMix) Epidemic Model
% Author: Michael Smah
% Institute: Institute for Global Pandemic Planning, Warwick Medical School,
%            & System Biology and Infectious Disease Epidemiology (SBIDER)
% Institution: University of Warwick, United Kingdom
% Project: Doctoral Project (2021-2025)
% Supervisors: Professor Kat Rock, Professor Anna Seale

%% Clear workspace and set up plotting defaults
clear all
close all
clc
% Define colormap and default plotting settings
cmap = colormap(parula(10));
set(gca, 'fontsize', 18)
set(0, 'defaultaxesfontsize', 18)
set(0, 'defaultlinelinewidth', 3)

%% Model setup: Initial conditions, parameters, and simulation settings
tic
% Population and initial conditions
N = 1000; % Total population size
S_in = N - 1; % Initial susceptible individuals
E_in = 0; % Initial exposed individuals
I_in = 1; % Initial infected individuals
R_in = 0; % Initial recovered individuals

% Clustering and contact parameters
NC = 200; % Number of clusters
clustersize = N / NC; % Average cluster size
Contacts = clustersize - 1; % Average daily contacts per individual
exten = 2; % Average number of external connections

% Epidemiological parameters
mu = 0.18; % Per-contact infection probability
gamma = 0.16; % Recovery probability
sigma = 0.15; % Probability of loss of latency period
kappa = 0.6031; % Fitted kappa value
C = Contacts; % Store contacts for model use
Maxtime = 80; % Simulation duration (days)

% Parameter structure for synergistic model
para = struct('mu', mu, 'gamma', gamma, 'N', N, 'sigma', sigma, ...
              'C', C, 'n', clustersize, 'K', NC, 'exten', exten, ...
              'kappa', kappa);

% Parameter structure for classical model (n = N, K = 1)
paraclassic = struct('mu', mu, 'gamma', gamma, 'N', N, 'sigma', sigma, ...
                     'n', N, 'K', 1, 'exten', N-1, 'kappa', kappa);

% Compute beta and contacts for classical model to match synergistic R0
paraclassic.beta = para.mu * (calculate_R0(para) * gamma) / ...
                   (mu * (1 - (exp(-(kappa * N * (N-1) / (1+N))))));
paraclassic.C = (calculate_R0(para) * gamma) / ...
                (mu * (1 - (exp(-(kappa * N * (N-1) / (1+N))))));
R0classic = paraclassic.beta / paraclassic.gamma; % Classical model R0

% Initial conditions structure
ICs_SIR = struct('S', S_in, 'I', I_in, 'R', R_in);

%% Run SIR models (SeRaMix and classical)
[ClassesSIR] = DETSIR(para, Maxtime, ICs_SIR);
[ClassesSIRClassic] = DETSIRClassic(paraclassic, Maxtime, ICs_SIR);

%% Simulate Non-Pharmaceutical Interventions (NPIs)
ChangeTime = 15; % Time to introduce NPIs (during exponential growth)

% Run models up to ChangeTime without NPIs
[ClassesSIR_preNPI] = DETSIR(para, ChangeTime, ICs_SIR);
[ClassesClassicCONT] = DETSIRClassic(paraclassic, ChangeTime, ICs_SIR);

% Store final states as initial conditions for NPI phase
ICs_postNPI = struct('S', ClassesSIR_preNPI.S(end), ...
                     'I', ClassesSIR_preNPI.I(end), ...
                     'R', ClassesSIR_preNPI.R(end));
ICsClassicNPI = struct('S', ClassesClassicCONT.S(end), ...
                       'I', ClassesClassicCONT.I(end), ...
                       'R', ClassesClassicCONT.R(end));

% NPI settings: 50% reduction in parameters
NPIsFactor = 0.5;

% NPI 1: Reduce contact rate (synergistic model)
para_NPI1 = para;
para_NPI1.C = para.C * (1 - NPIsFactor);
[ClassesSIR_NPIS] = DETSIR(para_NPI1, (Maxtime - ChangeTime), ICs_postNPI);

% NPI 2: Reduce external connections (synergistic model)
para_NPI2 = para;
para_NPI2.exten = (1 - NPIsFactor) * para.exten;
[SynClassesNPI2] = DETSIR(para_NPI2, (Maxtime - ChangeTime), ICs_postNPI);

% NPI 3: Reduce cluster size (synergistic model)
para_NPI3 = para;
para_NPI3.n = (1 - NPIsFactor) * para.n;
[SynClassesNPI3] = DETSIR(para_NPI3, (Maxtime - ChangeTime), ICs_postNPI);

% NPI 4: Reduce both cluster size and external connections
para_NPI3 = para;
para_NPI3.n = (1 - NPIsFactor) * para.n;
para_NPI3.exten = para.exten * (1 - NPIsFactor);
[SynClassesNPI4] = DETSIR(para_NPI3, (Maxtime - ChangeTime), ICs_postNPI);

% NPI: Reduce contact rate (classical model)
parac_NPI1 = paraclassic;
parac_NPI1.C = paraclassic.C * (1 - NPIsFactor);
[ClassesSIR_NPIC] = DETSIRClassic(parac_NPI1, (Maxtime - ChangeTime), ICsClassicNPI);

% Classical model with synergistic model’s C and mu
paraclassicx = paraclassic;
paraclassicx.C = para.C;
[ClassesSIRClassicx] = DETSIRClassic(paraclassicx, Maxtime, ICs_SIR);

%% Plot infection dynamics for varying cluster sizes
parancx = para;
lineStyles = {'-', '--', ':', '-.'}; % Line styles for plotting
colors = viridis(9); % Colormap for different cluster sizes
cluster_sizes = [2 4 5 10 20 25 40 50 N]; % Cluster sizes to simulate

figure;
hold on;
h = zeros(length(cluster_sizes), 1); % Preallocate handle array
for i = 1:length(cluster_sizes)
    CS = cluster_sizes(i);
    parancx.n = CS;
    if parancx.n == N
        parancx.exten = N - 1;
    else
        parancx.exten = para.exten;
    end
    [ClassesSIRCL] = DETSIR(parancx, Maxtime, ICs_SIR);
    h(i) = plot(ClassesSIRCL.t, ClassesSIRCL.I, ...
                'LineWidth', 4, ...
                'LineStyle', lineStyles{mod(i-1, length(lineStyles)) + 1}, ...
                'Color', colors(i,:));
end

% Plot classical model
h_classic = plot(ClassesSIRClassicx.t, ClassesSIRClassicx.I, ...
                 '--b', 'LineWidth', 4);

% Formatting plot
hold off;
ylabel('Infections', 'FontSize', 26, 'Interpreter', 'latex')
xlabel('Time (Days)', 'FontSize', 26, 'Interpreter', 'latex')
title('Infection dynamics for varying cluster sizes ($n$)', ...
      'FontSize', 26, 'Interpreter', 'latex')

% Create legend
legend_labels = arrayfun(@(x) sprintf('$n = %d$', x), cluster_sizes, ...
                         'UniformOutput', false);
legend_labels{end+1} = 'Classical model';
legend([fliplr(h); h_classic], legend_labels, ...
       'Interpreter', 'latex', 'FontSize', 20, 'Location', 'best');
set(gcf, 'position', [357 290 1054 805]);
set(gca, 'FontSize', 26);

%% Plot NPI effects
C = colororder("reef");
figure
plot(ClassesSIR_NPIS.t + ChangeTime, ClassesSIR_NPIS.I, ...
     'color', [62 150 81]./255, 'LineWidth', 5, ...
     'DisplayName', 'SeRaMix-EBM (contact (C) reduction)');
hold on
plot(ClassesSIR_NPIC.t + ChangeTime, ClassesSIR_NPIC.I, ...
     'color', [0.4940 0.1840 0.5560], 'LineWidth', 5, 'LineStyle', '--', ...
     'DisplayName', 'Classic-EBM ($\beta$ reduction)');
plot(SynClassesNPI2.t + ChangeTime, SynClassesNPI2.I, ...
     'Color', C(6,:), 'LineWidth', 5, ...
     'DisplayName', 'SeRaMix-EBM (external connection (x) reduction)');
plot(SynClassesNPI3.t + ChangeTime, SynClassesNPI3.I, ...
     'Color', C(1,:), 'LineWidth', 5, ...
     'DisplayName', 'SeRaMix-EBM (average cluster size (n) reduction)');
plot(SynClassesNPI4.t + ChangeTime, SynClassesNPI4.I, ...
     'Color', C(2,:), 'LineWidth', 5, ...
     'DisplayName', 'SeRaMix-EBM (all connections (n,x) reduction)');
plot(ClassesSIR.t, ClassesSIR.I, ...
     'color', [0.5 0.5 0.5], 'LineWidth', 5, ...
     'DisplayName', 'SeRaMix-EBM (without NPI)');
plot(ClassesSIRClassic.t, ClassesSIRClassic.I, ...
     'b', 'LineWidth', 5, 'LineStyle', '--', ...
     'DisplayName', 'Classic-EBM (without NPI)');

% Add intervention marker
xline_position = ChangeTime;
xline(xline_position, 'k--', 'LineWidth', 4, 'DisplayName', 'Start of Intervention');
text(xline_position - 2, 150, 'Start of Intervention', ...
     'FontSize', 24, 'Color', 'k', 'Rotation', 90, 'Interpreter', 'latex');

% Plot formatting
ylabel('Number infected', 'FontSize', 24, 'Interpreter', 'latex')
xlabel('Time (Days)', 'FontSize', 24, 'Interpreter', 'latex')
title({'Effect of $50\%$ reduction in contact parameters', ...
       'on transmission dynamics (EBM)'}, ...
      'Interpreter', 'latex', 'FontSize', 20)
legend('Location', 'southoutside', 'Orientation', 'horizontal', ...
       'NumColumns', 2, 'Interpreter', 'latex', 'FontSize', 18);
set(gcf, 'position', [357 290 1054 805]);

%% Compute R0 and epidemic thresholds
% Calculate R0 for base parameters
R0 = calculate_R0(para);

% R0 as a function of cluster size (n)
R0paran = para;
R0n = zeros(1, 5);
for n = 1:5
    R0paran.n = n;
    R0n(n) = calculate_R0(R0paran);
end

% R0 as a function of external connections (x)
R0parax = para;
R0x = zeros(1, 5);
for x = 1:5
    R0parax.exten = x;
    R0x(x) = calculate_R0(R0parax);
end

% R0 as a function of cluster size (n) and external connections (x)
R0paranx = para;
R0nx = zeros(5, 6);
for n = 1:5
    R0paranx.n = n;
    for x = 1:6
        R0paranx.exten = x;
        R0nx(n, x) = calculate_R0(R0paranx);
    end
end

% R0 as a function of number of clusters (K) and contacts (C)
K_valuesx = [1, 100, 200, 250, 300, 400, 500];
paraKmn = para;
R0mn = zeros(length(K_valuesx), 10);
for NumClusters = 1:length(K_valuesx)
    paraKmn.n = para.N / K_valuesx(NumClusters);
    for avcontact = 1:10
        paraKmn.C = avcontact;
        R0mn(NumClusters, avcontact) = calculate_R0(paraKmn);
    end
end

% Epidemic threshold as a function of contacts (C)
parath = para;
epidemicthresholdx = zeros(1, 10);
for c = 1:10
    parath.C = c;
    epidemicthresholdx(c) = epidemicthreshold(parath);
end

%% Plot R0 heatmap
figure
imagesc(R0nx) % Plot R0 as a function of n and x
hold on
axis('xy') % Ensure y-axis is positive upward
hold off
yticks(1:5)
yticklabels({'1', '2', '3', '4', '5'})
ylabel('average cluster size (n)', 'FontSize', 26, 'Interpreter', 'latex')
xlabel('Average external connection (x)', 'FontSize', 26, 'Interpreter', 'latex')
title({'Color map representing the basic reproduction number', ...
       'with varying average cluster sizes and external connections'}, ...
      'FontSize', 26, 'Interpreter', 'latex')
colormap("viridis");
cb = colorbar('southoutside');
ylabel(cb, 'The basic reproduction number, $R_0$', ...
       'FontSize', 26, 'Interpreter', 'latex')
set(gcf, 'position', [357 290 1054 805])

% Finalize timing
toc