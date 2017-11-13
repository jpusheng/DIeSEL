function [estimate] = diesel(init, anchors,...
 distances, nodesV, number_of_anchors, time_window, time_step, MAXITER, stopping_criteria)
    
    % the anchors give us the dimension, [a(t=1); a(t=2) ...], a(t) = [a1(t), a2(2), ...]
    % the distances are stacked, horizontal distance matrices [D1 D2 D3 ...]
    % we assume that no edge is ever created or destroyed, the same set of edges is transported from the first time instant
    % time_window is our time window (1 = 1 timestep, etc)    
    % nodesV are time_window*dimension times columns of size n concatenated horizontally, the last column is the most recent
    %%%% corresponding to the displacements within time-instants, we still need to integrate distances across time, done in acc_nodes
    % assumptions: there's always at least one anchor-node edge
    % 
    % now we don't have the anchor speeds, take them off from the deltaVanchors
    % anchors now have the following structure:
    % pos1t=1 pos1t=2
    % pos2t=1 pos2t=2, etc
    
    NODE_NODE = 1;
    ANCHOR_NODE = 2; % for the projection function
    dimension = size(anchors,1)/number_of_anchors; %dimension
   
    d = dimension;
    number_of_points = size(distances,1);
    unknowns = number_of_points - number_of_anchors;
    
    check_matrix = distances(:, 1:number_of_points) > 0; %one snapshot
    
    %%%% code to generate Lipschitz constant, from sphereProjection
    
    nd = sum(check_matrix(number_of_anchors+1:end,number_of_anchors+1:end),2);
    dmax = max(nd);
  %% max |A_i|
    maxAi = max(sum(check_matrix(1:number_of_anchors,number_of_anchors+1:end),1));

    L = (2*dmax + maxAi) * time_window + 2; % adapted
    %%%% end of Lipschitz
    
    %generate total edge array
    
    edge_array = [];
    
    %counter = 0; % this is to separate the different type of edges, node-node, node-edge
    for ii = 1:number_of_points
        for jj = ii + 1:number_of_points
            if check_matrix(ii,jj)
                edge_array = [edge_array; [ii, jj]];
            end
        end
    end
    
  
    edge_array_nodes = []; %only between unknown nodes
    edge_array_nodes_anchors = []; % only between nodes and anchors, the unknowns are always the second index
    
    for ii = 1:size(edge_array,1) 
        if edge_array(ii,1) > number_of_anchors && edge_array(ii,2) > number_of_anchors
            edge_array_nodes = [edge_array_nodes; edge_array(ii,:)];
        end
        
        if edge_array(ii,1) <= number_of_anchors && edge_array(ii,2) > number_of_anchors
            edge_array_nodes_anchors = [edge_array_nodes_anchors; edge_array(ii,:)];
        
        elseif edge_array(ii,1) > number_of_anchors && edge_array(ii,2) <= number_of_anchors
            edge_array_nodes_anchors = [edge_array_nodes_anchors; edge_array(ii,:)];
        end
        
    end

    number_of_edges_nodes = size(edge_array_nodes,1);
    number_of_edges_anchors = size(edge_array_nodes_anchors, 1);
    
    
    C = zeros(number_of_edges_nodes, unknowns); %arc-node incidence matrix, only between nodes
    for ii = 1:size(edge_array_nodes,1);
        C(ii, edge_array_nodes(ii,1)- number_of_anchors) = 1;
        C(ii, edge_array_nodes(ii,2)- number_of_anchors) = -1;    
    end
  
    
    Id = eye(dimension);
    
    A = kron(C, Id);

    D_aux = []; % the vertical concatenation of identity matrices
    
    for ii = 1:time_window
       D_aux = [D_aux; eye(number_of_edges_nodes)];
    end
    
    D = kron(D_aux, Id); % to match size o A


    anchor_connectivity_array = sum(check_matrix(1:number_of_anchors, number_of_anchors + 1:end),1); 
        
    aux = [];
    
    for ii = 1:time_window
        aux = [aux; eye(sum(anchor_connectivity_array))]; % the equivalent of D matrix for anchors
    end
  
    aux = kron(aux,Id); % same procedure as line 89
    
    X = zeros(sum(anchor_connectivity_array) ,unknowns); % the selector
    index_counter = 0;
    for ii = 1:unknowns
        for jj = 1:anchor_connectivity_array(ii)
            X(index_counter + 1, ii) = 1;
            index_counter = index_counter + 1;
        end
    end
    
    X = kron(X, Id);
    
    E = aux*X; 

    %all matrices completed, now moving on to stacks:
    %we need to integrate the distances across the time
    
    %distances traversed against the fluid by nodes and current;
     
    nodesD = nodesV * time_step;
    
    %we assume that the last column is the most recent
    
    acc_nodes = [];

    
    for ii = 1:time_window
        auxD = zeros(unknowns,dimension);
        
        for jj = 1:ii
            auxD = auxD + nodesD(:,dimension*(jj-1) + 1:dimension*jj);        
        end     
        
        acc_nodes = [acc_nodes, auxD];
    end
    %acc_nodes gives us the integrations of speed along all the instances
    %of our window

    %now we make the appropriate subtractions
    deltaVnodes = [];
    deltaVanchors = [];
    
    for ii = 1:time_window
        aux = [];
        for jj = 1:number_of_edges_nodes
            index1 = edge_array_nodes(jj,1) - number_of_anchors;
            index2 = edge_array_nodes(jj,2) - number_of_anchors; % to offset
            mult1 = C(jj, index1);
            mult2 = -mult1;         
            aux = [aux; (mult1*acc_nodes(index1,(ii-1)*dimension + 1: ii*dimension) + mult2*acc_nodes(index2,(ii-1)*dimension + 1: ii*dimension))'];
        end
        deltaVnodes = [deltaVnodes; aux];
    end
    
 	% the stacks are completed, now moving to the iterations
    
    %matrices which may prove to be useful
    Iy = eye(time_window*dimension*number_of_edges_nodes); % corresponding to range correction vectors
    Iw = eye(time_window*dimension*size(edge_array_nodes_anchors,1));
    zerosy = zeros(size(Iy));
    zerosw = zeros(size(Iw));

    
    B = A'*D'*D*A + E'*E;
    P = [];
    
    %initialization of P
    for ii = 1:unknowns
        P = [P; init(:,ii)]; %rows of position columns
    end
    
    Y_init = [];
    W_init = [];
     
    %initialization of Y
    for ii = 1:number_of_edges_nodes

        indexp1 = edge_array_nodes(ii,1) - number_of_anchors;
        indexp2 = edge_array_nodes(ii,2) - number_of_anchors;
        %the first one is always + and the second one always -
        Y_init = [Y_init; (init(:,indexp1) - init(:,indexp2))];

    end
    
    Y_acc = [];
    for ii = 1:time_window
        Y_acc = [Y_acc; Y_init];
    
   end
    
    Y_init = Y_acc;
    
    Y_init = Y_init + deltaVnodes; % to account for the relative movements 
                        % during the time window
    
    % our edges with anchors-nodes are not the most suitable to continue
    % the algorithm, australian_edges stack them according to node index,
    % as opposed to anchor index
    australian_edges = [];
    for jj = 1:unknowns
        for ii = 1:size(edge_array_nodes_anchors,1)
            indexn = edge_array_nodes_anchors(ii,2) - number_of_anchors;
            
            if indexn == jj
                australian_edges = [australian_edges; edge_array_nodes_anchors(ii,:)];
            end
        end
    end
    
    for ii = 1:time_window 
        aux = [];
        for jj = 1:number_of_edges_anchors
            index_node = australian_edges(jj,2) - number_of_anchors;
            index_anchor = australian_edges(jj,1);
            
            aux = [aux; (acc_nodes(index_node, (ii-1)*dimension + 1: ii*dimension))'];
        end
        deltaVanchors = [deltaVanchors; aux];
    end
    % Initialization of W, working!
    
    for ii = 1:time_window % 
    
        for jj = 1:size(australian_edges,1)
            
            indexa = australian_edges(jj,1);
            indexn = australian_edges(jj,2) - number_of_anchors;
            
            W_init = [W_init; (init(:,indexn) - anchors(dimension*(indexa-1)+1:dimension*indexa,ii))];
        end
    end
    
    W_ = W_init;
    % same procedure as line 206
    W_init = W_init + deltaVanchors;
    
    % see details about projection function below
    Y = projection(Y_init, edge_array_nodes, distances, time_window, number_of_points, dimension);
    W = projection(W_init, australian_edges, distances, time_window, number_of_points, dimension);
   
    
    alpha = []; %stack of anchor positions
    
     %initialization of alpha, according to australian_edges   
    for ii = 1:size(anchors,2)
        for jj = 1:size(australian_edges)
            indexa = australian_edges(jj,1);
            alpha = [alpha; anchors(dimension*(indexa-1) + 1:dimension*indexa , ii)];
        
        end
    end
    
    alpha = alpha - deltaVanchors;
    
    
    M1 = [D*A, -eye(time_window*number_of_edges_nodes*dimension, time_window*number_of_edges_nodes*dimension), zeros(time_window*dimension*number_of_edges_nodes,time_window*dimension*number_of_edges_anchors)];
    M2 = [E, -zeros(dimension*number_of_edges_anchors*time_window,dimension*number_of_edges_nodes*time_window), -Iw];
    
    M = M1'*M1 + M2'*M2;
    
    b = M2'*alpha - M1'*deltaVnodes; 

    % finally the algorithm


    for tt = 1:MAXITER
        P_new = (eye(size(B)) - (B/L))*P + (A'*D'/L)*(Y - deltaVnodes) + (E'/L)*(W + alpha );
        Y_new = ((L-1)/L)*Y + (D*A/L)*P + deltaVnodes/L;
        W_new = ((L-1)/L)*W + (E/L)*P - (1/L)*(alpha );
        z_old = [P;Y;W];
        P = P_new;
        Y = projection(Y_new, edge_array_nodes, distances, time_window, number_of_points, dimension);
        W = projection(W_new, australian_edges, distances, time_window, number_of_points, dimension);
        z_new = [P;Y;W];  
        
    end



    estimate = [];
    
    %add the contribution of speeds to the estimate
    for ii = 1:unknowns
        estimate = [estimate, P(dimension*(ii-1) + 1:dimension*ii, :) + acc_nodes(ii, dimension*(time_window - 1) + 1:end)'];
    end
    
   
    
    function [projected_array] = projection(previous_array, partial_edge_array, distances, time_window, number_of_points, dimension)
       %projects according to the type of array, partial_edge_array gives the indexes with which to normalize 
        array_dim = size(partial_edge_array,1);
        for ii = 1:time_window
            for jj = 1:array_dim
                indexr1 = partial_edge_array(jj,1);
                indexr2 = partial_edge_array(jj,2);
                actual_y = previous_array(dimension*array_dim*(ii-1) + dimension*(jj-1) + 1:...
                    dimension*array_dim*(ii-1) + dimension*jj,:);
                try
                    factor = distances(indexr1, (ii-1)*number_of_points + indexr2)/norm(actual_y);
                catch
                    factor = 0;
                end
                previous_array(dimension*array_dim*(ii-1) + dimension*(jj-1) + 1:...
                    dimension*array_dim*(ii-1) + dimension*jj,:) = factor*previous_array(dimension*array_dim*(ii-1) + dimension*(jj-1) + 1:...
                    dimension*array_dim*(ii-1) + dimension*jj,:);
            end
        end
        projected_array = previous_array; % this should work
    end
end