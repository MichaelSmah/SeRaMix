clear all
close all
clc

    % Define parameters
    N = 200; % Total population size
    mu = 0.2094; % per contact infection probability
    gamma = 0.1428; % recovery probability
    C = 3; % Average contacts per day
    e = 2; %number of external potential contacts
    K = 20; %number of clusters
    MaxTime = 50; %number to stop simulation

    %store parameters as a structure
    para = struct('N',N, 'mu',mu,'gamma',gamma,'C',C,'K',K,'e',e);

    
    %generate Potentially-Recurrent Contact Network (PRCN) social circle
    hypergraph = PRCN(100, 10, 1);

    %run the epidemic simulation ON Potentially-Recurrent Contact Network (PRCN)
    EpidemicPRCN2(para,MaxTime);