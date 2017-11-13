function [pos_acc, vel_acc] = trajectoryLaps(n_iter, rr, rl, separation, ...
                initial_pos, offset1, offset2, time_step);

pos_acc = [];

% I can generalize this, but let's to this some other time

%initial config:
previous_position = initial_pos;
previous_mode = 1;
previous_center_of_rotation = 0;
previous_delta_angle = 0;
vel_acc = [];
radius_left = rl;
radius_right = rr;
mm = [];

rl2 = rl + separation;
rr2 = rr + separation;

rl3 = rl - separation;
rr3 = rr - separation;
p2 = initial_pos + offset1;
p3 = initial_pos + offset2;

m2 = 1;
pcr2 = 0;
pda2 = 0;
m3 = 1;
pcr3 = 0;
pda3 = 0;
original_s = 1;

for ii = 1:n_iter
    
    [position, mode, velocity, center_of_rotation, delta_angle] = ...
    movementDataLaps(previous_position, previous_mode, previous_center_of_rotation, ...
        previous_delta_angle, time_step, radius_left, radius_right, radius_left, radius_right,original_s);
    [p2, m2, v2, pcr2, pda2] =  movementDataLaps(p2, m2, pcr2, pda2, time_step, rl2, rr2, radius_left, radius_right , original_s);
    [p3, m3, v3, pcr3, pda3] =  movementDataLaps(p3, m3, pcr3, pda3, time_step, rl3, rr3, radius_left, radius_right, original_s);
    vel_acc = [vel_acc,[ velocity; v2; v3]];
    pos_acc = [pos_acc, [position; p2;p3]];
    previous_position = position;
    previous_mode = mode;
    previous_center_of_rotation = center_of_rotation;
    previous_delta_angle = delta_angle;
    mm = [mm mode];
end

%n_iter = 500;
%time_step = 1;
%speed_2d = 0.02;
%speed_z = -0.01;
%radius = 0.5;
%initial_z = 1;
%initial_phase = 0;
%center_2d = [0; 0];
%mode = 1;
%[pos_acc, vel_acc] = movementDataHelix(n_iter, time_step,...
%         speed_2d, speed_z, radius, initial_z, initial_phase, center_2d, mode);