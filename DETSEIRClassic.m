
%This function computes the deterministic simulation SEIR of the Synergistic
%Disadvantage Epidemic Model with clusters
%Authour: Michael Smah
%Institute for Global Pandemic Planning, Warwick Medical School, & System
%Biology and Infectious Disease Epidemiology (SBiDER)
%University of Warwick, United Kingdom
%Doctoral Project (2021-2025)

function [Classes] = DETSEIRClassic(para,Maxtime,ICs)

% Starting point of the simultion from the ICs copied to a new structure called Classes
Classes = ICs;
Classes.t = 0;

% Specify current states 
S = ICs.S; 
E = ICs.E;
I = ICs.I; 
R = ICs.R; 
N = para.N;

t = 0; % Initial time
while (t < Maxtime) %Condition to run the model 

     %uncomment this part to run the linearised FOI
    % ProbInf = (I/ N);
    % beta = para.beta;
    % % compute the probability of being infected after making Contacts in the population
    % P_t = beta*ProbInf; 



    ProbInf = (I/ N);%Per-contact probability of choosing infected person

    % compute the probability of being infected after making C contacts in the population
    P_t = 1 - (1 - para.mu*ProbInf).^para.C;


    % Compute the proportion of S, I and R individuals at time t
    Snew = S - S.*P_t;
    Enew = (1 - para.sigma).*E + S.*P_t;
    Inew = (1 - para.gamma).*I + para.sigma.*E;
    Rnew = para.gamma*I + R;


    % Update each quantity
    t = t+1;
    S = Snew;
    E = Enew;
    I = Inew;
    R = Rnew;

    % Save information in the Classes structure by extending the matrices of the
    % model state and the associated time
    Classes.t = [Classes.t t];
    Classes.S = [Classes.S S];
    Classes.E = [Classes.E E];
    Classes.I = [Classes.I I];
    Classes.R = [Classes.R R];
end

% Calculate the epidemic size
Classes.Patch_Epsize = sum(R,2);
% Calculate the global epidemic size
Classes.Global_Epsize = sum(Classes.Patch_Epsize);






