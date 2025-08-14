function [I_counts_total, IcountNew, stored_t2, stored_x2] = simulateSIRIBM(para)
   
%redefine parameters
N = para.N;
C = para.C;
mu = para.mu;
gamma = para.gamma;
MaxTime = para.MaxTime;
Realizations = para.Realizations;
OutbreakThreshold = para.OutbreakThreshold;
group_size = para.n;
num_external_neighbours_range = para.exten;


% Arrays to store all S, I, R counts over all num_external_neighbours
    S_counts_total = [];
    I_counts_total = [];
    R_counts_total = [];
    
    % Arrays to store S, I, R counts over time for successful realizations
    S_counts_all = [];
    I_counts_all = [];
    R_counts_all = [];
    
    %Array to store new infections over all realisations
    Inew_counts_all = [];
    
    % Initialize storage for confidence intervals
    stored_t2 = cell(length(num_external_neighbours_range), 1);
    stored_x2 = cell(length(num_external_neighbours_range), 1);
    
    %Arrays to store new infections
    IcountNew = cell(length(num_external_neighbours_range), 1);
    
    % Loop over different num_external_neighbours and create network connectivity
    for idx = 1:length(num_external_neighbours_range)
        num_external_neighbours = num_external_neighbours_range(idx); %current number of external neighbours
        
        % Generate group assignments to assign agents to their clusters
        group_assignments = reshape(randperm(N), group_size, [])';
        
        % Preallocate neighbor lists
        primary_neighbours = cell(N, 1);
        external_neighbours_per_agent = cell(N, 1);
        potential_neighbours = nan(N, group_size - 1 + num_external_neighbours);
        
        % Assign primary neighbors within each group
        for i = 1:size(group_assignments, 1)
            members = group_assignments(i, :);
            for j = 1:group_size
                agentA = members(j);
                primary_neighbours{agentA} = members(members ~= agentA); %skip self
            end
        end
        
        % Track available slots for external connections
        available_slots = repmat(num_external_neighbours, N, 1);  % Each agent can have 'num_external_neighbours' links
        
        % External neighbor assignment (bidirectional)
        while any(available_slots > 0)  % Continue until all slots are filled
            % Find agents that still need neighbors
            open_agents = find(available_slots > 0);
            
            % If too few open agents remain, pair them together
            if length(open_agents) < 5 && length(open_agents) > 1
                for j = 1:length(open_agents)-1
                    agentA = open_agents(j);
                    agentB = open_agents(j+1);
                    
                    if length(external_neighbours_per_agent{agentA}) < num_external_neighbours && ...
                            length(external_neighbours_per_agent{agentB}) < num_external_neighbours
                        external_neighbours_per_agent{agentA} = [external_neighbours_per_agent{agentA}, agentB];
                        external_neighbours_per_agent{agentB} = [external_neighbours_per_agent{agentB}, agentA];
                        
                        available_slots(agentA) = max(available_slots(agentA) - 1, 0);
                        available_slots(agentB) = max(available_slots(agentB) - 1, 0);
                    end
                end
                break; % Stop further iterations
            end
            
            % Ensure at least two agents remain open
            if length(open_agents) < 2
                break;
            end
            
            % Randomly shuffle open agents
            open_agents = open_agents(randperm(length(open_agents)));
            
            % Process pairs of open agents
            i = 1;
            while i < length(open_agents)  % Ensure we pick pairs
                agentA = open_agents(i);
                
                % Get list of valid external neighbors
                excluded_agents = unique([agentA, primary_neighbours{agentA}, external_neighbours_per_agent{agentA}]);
                excluded_agents = excluded_agents(excluded_agents >= 1 & excluded_agents <= N);
                
                % Get possible neighbors
                possible_neighbors = setdiff(1:N, excluded_agents);
                possible_neighbors = possible_neighbors(available_slots(possible_neighbors) > 0);
                
                % If no valid neighbors exist, move to the next agent
                if isempty(possible_neighbors)
                    i = i + 1;
                    continue;
                end
                
                % Pick a random external neighbor
                agentB = possible_neighbors(randi(length(possible_neighbors)));
                
                % Assign bidirectional external connection
                if length(external_neighbours_per_agent{agentA}) < num_external_neighbours && ...
                        length(external_neighbours_per_agent{agentB}) < num_external_neighbours
                    external_neighbours_per_agent{agentA} = [external_neighbours_per_agent{agentA}, agentB];
                    external_neighbours_per_agent{agentB} = [external_neighbours_per_agent{agentB}, agentA];
                    
                    % Reduce available slots
                    available_slots(agentA) = max(available_slots(agentA) - 1, 0);
                    available_slots(agentB) = max(available_slots(agentB) - 1, 0);
                end
                
                % Ensure that both agents are fully assigned before removing them
                if available_slots(agentA) == 0 && available_slots(agentB) == 0
                    open_agents(open_agents == agentA) = [];
                    open_agents(open_agents == agentB) = [];
                end
                
                % Break if all slots are filled
                if all(available_slots == 0)
                    break;
                end
                
                i = i + 1;  % Move to next agent
            end
        end
                
        % Ensure all agents have the correct number of external neighbors
        for agentA = 1:N
            % Get the total number of neighbors (primary + external)
            total_neighbors = length(primary_neighbours{agentA}) + length(external_neighbours_per_agent{agentA});
            
            % If the number of neighbors is less than expected, pad with zeros
            if total_neighbors < num_external_neighbours
                % Calculate how many zeros are needed
                zeros_needed = num_external_neighbours - total_neighbors;
                
                % Pad the external neighbors list with zeros
                external_neighbours_per_agent{agentA} = [external_neighbours_per_agent{agentA}, zeros(1, zeros_needed)];
            end
            % Compute the total number of neighbors
            total_neighbors = length(primary_neighbours{agentA}) + length(external_neighbours_per_agent{agentA});
            
            % Ensure potential_neighbours can hold all neighbors (primary + external)
            if total_neighbors > size(potential_neighbours, 2)
                % Dynamically resize row for agentA if needed
                potential_neighbours(agentA, 1:total_neighbors) = [primary_neighbours{agentA}, external_neighbours_per_agent{agentA}];
            elseif total_neighbors < size(potential_neighbours, 2)
                % If fewer neighbors than expected, pad with zeros
                potential_neighbours(agentA, :) = [primary_neighbours{agentA}, external_neighbours_per_agent{agentA}, ...
                    zeros(1, size(potential_neighbours, 2) - total_neighbors)];
            else
                % Assign neighbors normally
                potential_neighbours(agentA, :) = [primary_neighbours{agentA}, external_neighbours_per_agent{agentA}];
            end
        end
        
        % Loop over realizations
        successful_realizations = 0; % Counter for successful outbreaks
        while successful_realizations < Realizations %run simulation until the number of realisations is achieved
            % Display progress
            fprintf('Outside neighbours = %d, Attempting realization %d of %d\n', num_external_neighbours_range(idx)...
                , successful_realizations + 1, Realizations);
            
            % Initialize health states (S=0, I=1, R=2)
            health_matrix = zeros(MaxTime+1, N); % Store health states over time
            initial_infected = randi(N); % Pick one random agent to be infected
            health_matrix(1, initial_infected) = 1; % Set initial infected agent
            
            % Arrays to store S, I, R counts over time
            S_counts = zeros(1, MaxTime+1);
            I_counts = zeros(1, MaxTime+1);
            R_counts = zeros(1, MaxTime+1);

            INew_counts = zeros(1, MaxTime+1); % Initialize INew_counts
            
            % Count initial S, I, R
            S_counts(1) = sum(health_matrix(1, :) == 0); % Susceptible
            I_counts(1) = sum(health_matrix(1, :) == 1); % Infected
            R_counts(1) = sum(health_matrix(1, :) == 2); % Recovered
            
            % Simulation loop
            for t = 2:MaxTime+1

                % Uncomment this part to implement NPI (contact reduction)
                % if t>=15
                %     C = 0.5*para.C;
                % end

                % Copy the previous time step's health states
                health_matrix(t, :) = health_matrix(t-1, :);
                
                inew_counter = 0; % Initialize before the agent loop
                % Update health states based on infection and recovery rates
                for agent = 1:N
                    if health_matrix(t-1, agent) == 0  % Susceptible
                        % Get the neighbours of the current agent
                        neighbours = potential_neighbours(agent, :);
                        
                        % Remove any invalid indices (e.g., zeros or out-of-bounds values)
                        neighbours = neighbours(neighbours > 0 & neighbours <= N);
                        
                        %% Use this part (actual selection) when using integer C.
                        % Randomly select C unique contacts from the neighbours (without replacement)
                        selected_indices = randperm(length(neighbours), C);
                        selected_neighbours = neighbours(selected_indices);

                        % Count how many of the C selected contacts are infected
                        infected_in_sample = sum(health_matrix(t, selected_neighbours) == 1);

                        % Compute the probability of infection from contact
                        infection_prob = mu*infected_in_sample;

                        % %% Use this part (approximate) when using non-integer C
                        % % Count infected neighbours
                        % infected_neighbours_count = sum(health_matrix(t, neighbours) == 1);
                        % 
                        % % Probability of infection
                        % infection_prob = 1 - (1 - mu*(infected_neighbours_count / length(neighbours)))^C;
                        
                        %% Determine if the agent becomes infected or recovered
                        if rand < infection_prob
                            health_matrix(t, agent) = 1; % Become infected
                            inew_counter = inew_counter+1; %update new infections count
                        end
                    elseif health_matrix(t-1, agent) == 1 % Infected
                        if rand < gamma
                            health_matrix(t, agent) = 2; % Recover
                        end
                    end
                end
                
                % Count S, I, R at the current time step
                S_counts(t) = sum(health_matrix(t, :) == 0); % Susceptible
                I_counts(t) = sum(health_matrix(t, :) == 1); % Infected
                R_counts(t) = sum(health_matrix(t, :) == 2); % Recovered
                INew_counts(t) = inew_counter; % Newly Infected
            end
            
            if idx <= group_size
                OutbreakThreshold = group_size+num_external_neighbours-1;
            else
                OutbreakThreshold = OutbreakThreshold;
            end
            
            % Check if the outbreak is successful
            if R_counts(end) >= OutbreakThreshold
                % Increment successful realizations counter
                successful_realizations = successful_realizations + 1;
                % Store the results for this successful realization
                S_counts_all(successful_realizations, :) = S_counts;
                I_counts_all(successful_realizations, :) = I_counts;
                R_counts_all(successful_realizations, :) = R_counts;
                Inew_counts_all(successful_realizations, :) = INew_counts;
            else
                fprintf('Outbreak failed (Final Size < %d), discarding this realization.\n', OutbreakThreshold);
            end
        end
        
        % Calculate mean populations over valid realizations for this num_external_neighbours_range
        S_counts_total(idx, :) = mean(S_counts_all, 1);
        I_counts_total(idx, :) = mean(I_counts_all, 1);
        R_counts_total(idx, :) = mean(R_counts_all, 1);

        %store new infections for this external connections across all
        %realisations
        IcountNew{idx} = Inew_counts_all;
        
        %Times to plot shaded (flip) plot of the IBM Confidence Interval (C.I)
        t = 0:MaxTime;
        %The times forward then backward as a row vector
        t2 = [t flip(t)];
        
        x = I_counts_all;
        %The upper interval forward and lower 95% backward as a row vector
        x2 = [quantile(x,0.025) flip(quantile(x,0.975))];
        
        %Store t2, x2, for each e in the loop
        stored_t2{idx} = t2;
        stored_x2{idx} = x2;
    end
end