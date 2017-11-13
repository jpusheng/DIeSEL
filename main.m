clear all
clc

TIME_WINDOW = 5; % window width
time_window = TIME_WINDOW;
MAXITER = 400; 
pos_acc = []; % will contain all the exact locations for anchors & nodes
est_acc = []; % will contain all the estimated locations for the nodes
n_iter = 10; % size of the simulation

range_deviation = 0.5;
speed_deviation = 0.01;
angle_deviation =  2;

%script with the sample trajectory generation
sample_trajectory;

%add noise here for the speeds
stopping_criteria = unknowns*1e-4;

vel_nodes_n = zeros(size(vel_nodes));
vel_anchors_n = zeros(size(vel_anchors));


for ii = 1:n_iter
    for jj = 1:unknowns
        value = vel_nodes((jj-1)*dimension + 1:jj*dimension,ii);
        coeff = norm(value) + speed_deviation*randn();
        noisy_speed = (coeff/norm(value))*value;
        vel_nodes_n((jj-1)*dimension + 1:jj*dimension,ii) = noisy_speed;
    end
end

vel_nodes_n(:,1) = [];

initial_position_noise_deviation = 2;
pos_init = pos_nodes(:,1) + initial_position_noise_deviation * randn(size(pos_nodes(:,1)));

% we initiate our position stacks with the first instance
anchor_positions_stack = anchors(:,1);
position_estimates_stack = pos_init*ones(1,time_window);
distance_matrices_stack = [];

est_acc = pos_init;


%produce all the distance matrices (with noise) for all the instances,
%stacked horizontally
for ii = 1:n_iter
    points = [];
    
    for jj = 1:number_of_anchors
        points = [points anchors((jj-1)*dimension + 1:jj*dimension, ii)];
    end
    
    for jj = 1:unknowns
        points = [points pos_nodes((jj-1)*dimension + 1:jj*dimension, ii)];
    end
    
    d_aux = zeros(number_of_points, number_of_points);
    
    for kk = 1:number_of_points
        for hh = kk+1:number_of_points
            d = norm(points(:,kk) - points(:,hh)) + range_deviation*randn();
            if d <= 0
                d = 1e-5; % to avoid issues in the algorithm
            end
            d_aux(hh,kk) = d;
            d_aux(kk,hh) = d;
        end
    end
    
    distance_matrices_stack = [distance_matrices_stack d_aux]; 
end

for ii = 1:n_iter

    if ii <= TIME_WINDOW
        time_window = ii;
        %anchor_positions_stack = [anchor_positions_stack anchors(:,ii)];
    else
        time_window = TIME_WINDOW;
        %anchor_positions_stack = [anchor_positions_stack(:,2:end) anchors(:,ii)];
    end
    %produce the stuff to be submitted for the algorithm
    %we need time_window speeds, and ranges
    nodesV = [];

    %generate the speeds in the time window to send to the algorithm
    % one row per node, speeds stacked according to their time instant
    % at beginning we add zero speeds
    for jj = 1:unknowns
        aux = [];
        for kk = 1:time_window 
            try 
                toadd = vel_nodes_n((jj-1)*dimension+1:jj*dimension, ii+(kk-time_window) -1); %the snapshot; 
            catch
                toadd = zeros(1,dimension)'; %in case the index is negatives, which                                            %happens in the beggining
            end
            aux = [aux; toadd];    
        end
        nodesV = [nodesV; aux'];
    end
    
    % provide the stack of time_window distances to be sent to the algorithm,
    % stacked horizontally
    distances = [];
    
    if ii >= time_window 
        distances = distance_matrices_stack(:,(number_of_points*(ii-time_window )) + 1:number_of_points*(ii));
    else
        number_of_delays = time_window - ii ;
        for cc = 1:number_of_delays
            distances = [distances, distance_matrices_stack(:,1:number_of_points)];
        end 
        distances = [distances, distance_matrices_stack(:,1:number_of_points*(time_window - number_of_delays))];
    end

    
    %initialization to be sent to the algorithm, works as a revolving door,
    % the one on the left is the oldest in the time window and will be used for
    % initialization of the nodes
    init = [];
    for yy = 1:unknowns
        init = [init, position_estimates_stack(dimension*(yy-1) + 1:dimension*yy,1)];
    
    end
    
    Vstack = nodesV;

    [est] = diesel(init, anchor_positions_stack,...
        distances, Vstack, number_of_anchors, time_window, time_step, MAXITER);


   % PP = [PP P];
    toadd = [];  
    for rr = 1:unknowns
        toadd = [toadd; est(:,rr)];  
        
    end
    
    %revolving door pt2: we add the new estimate to the right
    position_estimates_stack = [position_estimates_stack(:,2:end), toadd];

    est_acc = [est_acc, toadd];  
    

    if ii < TIME_WINDOW && ii <n_iter
        time_window = ii;
        anchor_positions_stack = [anchor_positions_stack anchors(:,ii+1)];
    elseif ii < n_iter;
        time_window = TIME_WINDOW;
        anchor_positions_stack = [anchor_positions_stack(:,2:end) anchors(:,ii+1)];
    end
   
   %end_EKF
   
end
  
    