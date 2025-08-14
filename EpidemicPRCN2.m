function EpidemicPRCN2(para, MaxTime)
    % Define parameters
    N = para.N; % Total population size
    mu = para.mu; % per contact infection probability
    gamma = para.gamma; % recovery probability
    C = para.C; % Average contacts per day
    K = para.K; % number of clusters
    p = para.e; % number of external potential contacts

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
    
    % Create a single figure with two subplots
    figure;
    set(gcf, 'position', [357 290 1054 805]);
    
    % Simulation loop
    G = graph(adjacencyMatrix);
    for t = 1:MaxTime
        new_states = states;
        for person = 1:N
            neighbors = find(adjacencyMatrix(person, :) == 1);
            if ~isempty(neighbors)
                % Randomly select C neighbors without replacement
                selected_neighbors = randsample(neighbors, min(para.C, length(neighbors)), false);
                % Count infected neighbors among the selected
                infected_contacts = sum(states(selected_neighbors) == 1);
                infection_prob = mu * infected_contacts;
            else
                infection_prob = 0; % No neighbors, no infection
            end
            
            if states(person) == 0 && rand < infection_prob
                new_states(person) = 1;
            elseif states(person) == 1 && rand < gamma
                new_states(person) = 2;
            end
        end
        
        % Update states
        states = new_states;
        
        % Count states
        susceptible_count(t) = sum(states == 0);
        infected_count(t) = sum(states == 1);
        recovered_count(t) = sum(states == 2);
        
        % Clear the figure to prevent overlapping plots
        clf;
        
        % Subplot 1: Network visualization
        subplot(1, 2, 1);
        colorMap = [0 0 1; 1 0 0; 0 1 0]; % Blue (susceptible), Red (infected), Green (recovered)
        nodeColors = colorMap(states + 1, :);
        plot(G, 'XData', nodePositions(:,1), 'YData', nodePositions(:,2), 'NodeColor', nodeColors, 'MarkerSize', 6);
        title(['Epidemic Simulation at Time ', num2str(t)]);
        
        % Subplot 2: Epidemic trends (real-time update)
        subplot(1, 2, 2);
        plot(1:t, susceptible_count(1:t), 'b', 1:t, infected_count(1:t), 'r', 1:t, recovered_count(1:t), 'g', 'LineWidth', 4);
        legend('Susceptible', 'Infected', 'Recovered');
        xlabel('Time'); ylabel('Population');
        title('Epidemic Spread Over Time');
        sgtitle('Semi Random Mixing Epidemic Model','FontSize', 24)
        % Pause to visualize the update
        pause(0.3);
        
        % Check if there are no infected people after visualization
        if infected_count(t) == 0
            break;
        end
    end
end