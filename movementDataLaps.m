function [position, mode, velocity, center_of_rotation, delta_angle] =...
 movementDataLaps(previous_position, previous_mode, previous_center_of_rotation, ...
 previous_delta_angle, time_step, radius_left, radius_right, ref_left,ref_right, original_s)

% center_of_rotation is only used if the next mode was rotation, previous_center_of_rotation is only used 
% in case the previous mode was rotation 
% 4 modes:
% 1 - going left
% 2 - going right
% 3 - turn counterclockwise
% 4 - turn clockwise ///also counterclockwise, but up

% bounds: we can only move left until 0.3, move right until 0.75
% the radii for counterclockwise mode is 0.2, for clockwise is 0.15
% the turn modes can only turn at an maximum of half a circumference
% the speed is 0.05 units/sec
% we can only move on the following order: 1 - 3 - 2 - 4 - 1
% only one transition is permitted! else the code fails



left_mode = 1;
right_mode = 2;
counterclockwise_mode = 3;
clockwise_mode = 4;

bound_left = 0;
bound_right = 100;
speed = original_s;

v_rot_coeff_r = radius_right / ref_right;
v_rot_coeff_l = radius_left / ref_left;

switch previous_mode
    case left_mode
    
        if previous_position(1) - speed*time_step > bound_left % i.e. we don't screw up te bounds   

            position = [previous_position(1) - speed*time_step; previous_position(2)]; 
            mode = previous_mode; % we keep going left
            center_of_rotation = 0; % we don't need, it's a random number
            velocity = (position-previous_position)/time_step;
            delta_angle = 0; % random number
        else
        
            mode = counterclockwise_mode;
            extra_time = time_step - (previous_position(1) - bound_left)/speed; % the time left
            speed = v_rot_coeff_l * speed;                                                                    % after we reached the bound
            center_of_rotation = [bound_left; previous_position(2) - radius_left];
            omega = speed/radius_left; %rotation speed
            delta_angle = omega * extra_time; %traversed angle
            position = [center_of_rotation(1) - radius_left*sin(delta_angle); center_of_rotation(2) + radius_left*cos(delta_angle)];
            velocity = (position-previous_position)/time_step;
        end   
    case right_mode
        if previous_position(1) + speed*time_step < bound_right % i.e. we don't screw up te bounds
            position = [previous_position(1) + speed*time_step; previous_position(2)]; 
            mode = previous_mode; % we keep going left
            center_of_rotation = 0; % we don't need, it's a random number
            velocity = (position-previous_position)/time_step;
            delta_angle = 0; % random number
        else
            mode = clockwise_mode;
            extra_time = time_step - (bound_right - previous_position(1))/speed; % the time left
                                                                              % after we reached the bound
            speed = v_rot_coeff_r * speed;
            center_of_rotation = [bound_right; previous_position(2) + radius_right];
            omega = speed/radius_right; %rotation speed
            delta_angle = omega * extra_time; %traversed angle
            position = [center_of_rotation(1) + radius_right*sin(delta_angle); center_of_rotation(2) - radius_right*cos(delta_angle)];
            velocity = (position-previous_position)/time_step;
        end 
    
    case counterclockwise_mode
        speed = v_rot_coeff_l * speed;
        omega = speed/radius_left;
        if previous_delta_angle + time_step*omega < pi
            mode = previous_mode;
            delta_angle = previous_delta_angle + time_step*omega;
            center_of_rotation = previous_center_of_rotation;
            position = [center_of_rotation(1) - radius_left*sin(delta_angle); center_of_rotation(2) + radius_left*cos(delta_angle)];
            velocity = (position-previous_position)/time_step;
        else
            extra_time = (previous_delta_angle + time_step*omega - pi)/omega;
            mode = right_mode;
            center_of_rotation = 0;
            speed = original_s;
            position = [previous_center_of_rotation(1) + extra_time*speed; previous_center_of_rotation(2) - radius_left];
            velocity = (position-previous_position)/time_step;
            delta_angle = 0;        
        end
        
        
    case clockwise_mode
        speed = v_rot_coeff_r * speed;
        omega = speed/radius_right;
        if previous_delta_angle + time_step*omega < pi
            mode = previous_mode;
            delta_angle = previous_delta_angle + time_step*omega;
            center_of_rotation = previous_center_of_rotation;
            position = [center_of_rotation(1) + radius_right*sin(delta_angle); center_of_rotation(2) - radius_right*cos(delta_angle)];
            velocity = (position-previous_position)/time_step;
        else
            extra_time = (previous_delta_angle + time_step*omega - pi)/omega;
            mode = left_mode;
            center_of_rotation = 0;
            speed = original_s;
            position = [previous_center_of_rotation(1) - extra_time*speed; previous_center_of_rotation(2) + radius_right];
            velocity = (position-previous_position)/time_step;
            delta_angle = 0;        
        end
end 
