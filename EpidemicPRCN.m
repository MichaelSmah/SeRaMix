function EpidemicPRCN(para, MaxTime)

    % Define parameters
    N = para.N; % Total population size
    mu = para.mu; % per contact infection probability
    gamma = para.gamma; % recovery probability
    C = para.C; % Average contacts per day
    K = para.K; %number of clusters
    p = para.e; %number of external potential contacts


    % Ensure N is evenly distributed across K groups
    if mod(N, K) ~= 0
        error('N must be evenly divisible by K');
    end
    
    n = N / K; % Group size
    adjacencyMatrix = zeros(N, N); % Initialize adjacency matrix
    
    % Assign people to groups
    groups = reshape(1:N, [n, K]);
    
    % Connect people within the same group
    for k = 1:K
        members = groups(:, k);
        for i = 1:n
            for j = i+1:n
                adjacencyMatrix(members(i), members(j)) = 1;
                adjacencyMatrix(members(j), members(i)) = 1;
            end
        end
    end
    
    % Connect each person to p random people from other groups
    for k = 1:K
        members = groups(:, k);
        otherMembers = setdiff(1:N, members); % People not in the same group

        for i = 1:n
            person = members(i);
            % randomConnections = randsample(otherMembers, e, false); % Pick e unique connections
            if K > 1 && ~isempty(otherMembers)
                randomConnections = randsample(otherMembers, min(p, length(otherMembers)), false);
                for j = 1:length(randomConnections)
                    adjacencyMatrix(person, randomConnections(j)) = 1;
                    adjacencyMatrix(randomConnections(j), person) = 1;
                end
            end

        end
    end

    % Generate node positions to keep groups separate
    nodePositions = zeros(N, 2);
    for k = 1:K
        theta = (k - 1) * (2 * pi / K); % Assigning each group an angle
        groupCenter = [cos(theta), sin(theta)] * 10; % Spacing groups apart
        
        for i = 1:n
            angleOffset = (i - 1) * (2 * pi / n); % Distribute within group
            nodePositions(groups(i, k), :) = groupCenter + [cos(angleOffset), sin(angleOffset)];
        end
    end
    
   
    states = zeros(N, 1); % 0: Susceptible (blue), 1: Infected (red), 2: Recovered (green)
    
    % Randomly infect one person at time 0
    initial_infected = randi(N);
    states(initial_infected) = 1;
    
    % Record numbers
    susceptible_count = zeros(MaxTime, 1);
    infected_count = zeros(MaxTime, 1);
    recovered_count = zeros(MaxTime, 1);
    
    % Simulation loop
    figure;
    G = graph(adjacencyMatrix);
    for t = 1:MaxTime
        new_states = states;
        for person = 1:N
            neighbors = find(adjacencyMatrix(person, :) == 1);
            infected_neighbors_count = sum(states(neighbors) == 1);
            infected_contacts = binornd(C, infected_neighbors_count / length(neighbors));
            infection_prob = 1 - (1 - mu)^infected_contacts;
            
            if states(person) == 0 && rand < infection_prob
                new_states(person) = 1;
            elseif states(person) == 1 && rand < gamma
                new_states(person) = 2;
            end
        end
        states = new_states;
        
        % Count states
        susceptible_count(t) = sum(states == 0);
        infected_count(t) = sum(states == 1);
        recovered_count(t) = sum(states == 2);
        
        % Update graph visualization
        colorMap = [0 0 1; 1 0 0; 0 1 0]; % Blue (susceptible), Red (infected), Green (recovered)
        nodeColors = colorMap(states + 1, :);


        plot(G, 'XData', nodePositions(:,1), 'YData', nodePositions(:,2), 'NodeColor', nodeColors, 'MarkerSize', 6);
        title(['Epidemic Simulation at Time ', num2str(t)]);
        set(gcf,'position',[357         290        1054         805])
        pause(0.3);
    end
    
    % Plot epidemic trends
    figure;
    plot(1:MaxTime, susceptible_count, 'b', 1:MaxTime, infected_count, 'r', 1:MaxTime, recovered_count, 'g','LineWidth',4);
    legend('Susceptible', 'Infected', 'Recovered');
    xlabel('Time'); ylabel('Population Count');
    title('Epidemic Spread Over Time');
end
