
dimension = 2;
number_of_points = 4;
number_of_anchors = 2;
unknowns = number_of_points - number_of_anchors;
neighbors = number_of_points - 1;


previous_mode = 1;
previous_center_of_rotation = 0;
previous_delta_angle = 0;
time_step = 1;
vel_acc = [];

PP = [];

omega = pi/200;
rot_rad = 0.4;
rr = 29;
rl = 43;
separation = 8;
initial_pos = [50;70];
offset1 = [8;8];
offset2 = [8;-8];

[pos_acc, vel_acc] = trajectoryLaps(n_iter, rr, rl, separation, ...
                initial_pos, offset1, offset2, time_step);



%we have our trajectory now
%for simplification, we'll have three points doing the same trajectory, one of 
%them is the anchor and it'll have zero offset

%add current
initial_pos2 = initial_pos + [16;0];
[pos_acc2, vel_acc2] = trajectoryLaps(n_iter, rr, rl, separation, ...
                initial_pos2, offset1, offset2, time_step);
                
anchor2 = [pos_acc2(1:dimension, :)];
vel2 = vel_acc2(1:dimension,:);
current_speed = [0;0];[0.01; -0.02];
vel_acc = vel_acc; + current_speed;

anchors = [pos_acc(1:2,:); anchor2];
vel_anchors = [vel_acc(1:dimension*(number_of_anchors-1),:); vel2];

pos_nodes = pos_acc(dimension*(number_of_anchors-1) +1:end,:); % take away + 1 hot fix
vel_nodes = vel_acc(dimension*(number_of_anchors-1) +1 :end,:);