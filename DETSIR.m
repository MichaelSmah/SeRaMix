
%This function computes the deterministic simulation SIR of the Semi-Random Mixing
% Epidemic Model with clusters
%Authour: Michael Smah
%Institute for Global Pandemic Planning, Warwick Medical School, & System
%Biology and Infectious Disease Epidemiology (SBiDER)
%University of Warwick, United Kingdom
%Doctoral Project (2021-2025)

function [Classes] = DETSIR(para,Maxtime,ICs)

% Starting point of the simultion from the ICs copied to a new structure called Classes
Classes = ICs;
Classes.t = 0;

% Specify current states 
S = ICs.S; 
I = ICs.I; 
R = ICs.R; 
n = para.n;
N = para.N;
exten = para.exten;
kappa = para.kappa;
t = 0; % Initial time
while (t < Maxtime) %Condition to run the model 

    CX = 1 - exp(-(kappa*exten * n) / (1+n));%probability that groups connect to each other 


    ProbInf = (I/ N) .* CX;
    
    % compute the probability of being infected after making Contacts in the population
    P_t = 1 - (1 - para.mu*ProbInf).^para.C; 


    % Compute the proportion of S, I and R individuals at time t
    Snew = S - S.*P_t;
    Inew = S.*P_t + (1 - para.gamma).*I;
    Rnew = para.gamma*I + R;


    % Update each quantity
    t = t+1;
    S = Snew;
    I = Inew;
    R = Rnew;

    % Save information in the Classes structure by extending the matrices of the
    % model state and the associated time
    Classes.t = [Classes.t t];
    Classes.S = [Classes.S S];
    Classes.I = [Classes.I I];
    Classes.R = [Classes.R R];
end

% Calculate the epidemic size
Classes.Patch_Epsize = sum(R,2);
% Calculate the global epidemic size
Classes.Global_Epsize = sum(Classes.Patch_Epsize);






