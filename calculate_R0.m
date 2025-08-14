%Afunction for computing the numerical vale of the basic reproduction
%number R0
function R0 = calculate_R0(para)
    n = para.n;
    exten = para.exten;
    mu = para.mu;
    C = para.C;
    kappa = para.kappa;
    gamma = para.gamma;
      
    part1 = (mu*C)/gamma;
    part2 = 1 - exp(-kappa*exten .* n./ (n+1));%probability that groups connect to each other 
    R0 =  part1 * part2;
end

