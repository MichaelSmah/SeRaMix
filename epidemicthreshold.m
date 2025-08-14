%Afunction for computing the numerical value of the epidemic threshold
function Thresholdmu = epidemicthreshold(para)
    n = para.n;
    exten = para.exten;
    C = para.C;
    kappa = para.kappa;
    gamma = para.gamma;

    part2 = C*(1 - exp(-(kappa*exten * n) / (1+n)));
   Thresholdmu =  gamma / part2;
end
