function hypergraph = PRCN(N, K, p)
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
            randomConnections = randsample(otherMembers, p, false); % Pick p unique connections
            for j = 1:p
                adjacencyMatrix(person, randomConnections(j)) = 1;
                adjacencyMatrix(randomConnections(j), person) = 1;
            end
        end
    end
    
    % Output the adjacency matrix as the network representation
    hypergraph = adjacencyMatrix;
    
    % Generate node positions to keep groups separate
    nodePositions = zeros(N, 2);
    for k = 1:K
        theta = (k - 1) * (2 * pi / K); % Assigning each group an angle
        groupCenter = [cos(theta), sin(theta)] * 5; % Spacing groups apart
        
        for i = 1:n
            angleOffset = (i - 1) * (2 * pi / n); % Distribute within group
            nodePositions(groups(i, k), :) = groupCenter + [cos(angleOffset), sin(angleOffset)];
        end
    end
    
    % Visualize the graph with adjusted layout
    figure;
    G = graph(adjacencyMatrix);
    plot(G, 'XData', nodePositions(:,1), 'YData', nodePositions(:,2), 'MarkerSize', 15, 'NodeColor', 'b', 'LineWidth', 1, 'NodeLabel', []);
    title('Potentially Recurrent Contacts Network', 'Interpreter', 'latex');
    set(gcf, 'Position', [100, 100, 800, 600]);  % Adjust window size if needed

    figure;

    G = graph(adjacencyMatrix);

    % First plot with white (invisible) node color
    h = plot(G, ...
        'XData', nodePositions(:,1), ...
        'YData', nodePositions(:,2), ...
        'NodeCData', ones(numnodes(G),1), ...  % Dummy color data
        'MarkerSize', 15, ...
        'NodeColor', 'w', ...                  % Node fill color white
        'LineWidth', 1, ...
        'NodeLabel', []);                      % No labels

    % Add a ring-style node by using circle outlines
    % We simulate edge-only node look using a second plot on top
    hold on;
    scatter(nodePositions(:,1), nodePositions(:,2), ...
        200, ...             % Marker size (adjust as needed)
        'b', ...             % Edge color
        'o', ...             % Circle marker
        'LineWidth', 1.5, ...
        'MarkerFaceColor', 'none');  % No fill

    % Set tighter axis limits to reduce empty space
    xMargin = 0.1; % Small margin for x-axis
    yMargin = 0.1; % Small margin for y-axis
    xLimits = [min(nodePositions(:,1))-xMargin, max(nodePositions(:,1))+xMargin];
    yLimits = [min(nodePositions(:,2))-yMargin, max(nodePositions(:,2))+yMargin];
    axis([xLimits, yLimits]);

    axis off;

    title('Potentially Recurrent Contacts Network', 'Interpreter', 'latex');
    set(gcf, 'Position', [440,63,867,635]);  % Adjust window size if needed

end