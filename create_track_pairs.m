function [x, Y, Z, path_genuine, path_spoofed] = create_track_pairs(T, ...
    track_dynamics, spoof_parameters, track_conditions)

% [x, Y, Z, path_genuine, path_spoofed] = create_track_pairs(T, ...
% track_dynamics, spoof_parameters, track_conditions)
% 
% This function creates a pair of tracks for T time steps according to
% track_dynamics: parameters constant velocity model and the angualar bias
% spoof_parameters: spoofing parameters
% track_conditions: given locations and velocities
% 
% OUTPUTS:
% 
% x: Tx1 vector of states for the spoofing/genuine status of the target
% Y: Radar track in polar coordinates
% Z: AIS track in polar coordinates
% path_genuine: genuine track in cartesian coordinates
% path_spoofed: spoofed track in cartesian coordinates

spoof_create_mode = spoof_parameters.spoof_create_mode;

mu_0 = track_dynamics.mu0;
Cov_0 = track_dynamics.Cov0;
H = track_dynamics.H;
U = track_dynamics.U;
d = 2;
I_d = eye(d);
angle_comp = track_dynamics.angle_comp;
theta = track_dynamics.theta;
sigmaNoise = track_dynamics.sigmaNoise;

% Step A. Generate spoofing/no spoofing process
if spoof_create_mode == 1
    % This is according to the given transition matrix
    F = spoof_parameters.F;
    K = size(F, 1);
    x = zeros(1, T);
    for t = 1:T
        if t == 1
            x(t) = 1;
        else
            x(t) = randsample(1:K, 1, 'true', F(x(t-1), :));
        end
    end
else
    % Create spoofing once per target and with a certain probability
    spoof_time_coeff_vec = spoof_parameters.spoof_time_coeff_vec;
    probSpoofedTraj = spoof_parameters.probSpoofedTraj;
    if rand < probSpoofedTraj
        % determine the spoofing period
        cond = 0;
        while cond == 0
            t_temp = gamrnd(spoof_time_coeff_vec, 1);
            cond = t_temp(2)*T > spoof_parameters.min_spoof_duration;
        end
        cumsum_t_temp = cumsum(t_temp);
        t_spoof_begin = ceil(T*cumsum_t_temp(1)/cumsum_t_temp(end));
        t_spoof_end = ceil(T*cumsum_t_temp(2)/cumsum_t_temp(end));
        x = ones(1, T);
        x(t_spoof_begin:t_spoof_end) = 2;
    else
        x = ones(1, T); % set x as all genuine
    end
end

% Step B. Generate the tracks
% B.1 Generate the genuine track in cartesian coordinates
path_genuine = conditional_sampling(mu_0, Cov_0, H, U, T, track_conditions);

% B.2 Now generate the spoofed track
t_chp_b = find(diff(x) ~= 0 & x(1:end-1) == 1); % Identify change points in the trajectory
t_chp_e = find(diff(x) ~= 0 & x(1:end-1) == 2);
if x(end) == 2
    t_chp_e = [t_chp_e T];
end
n_chp = length(t_chp_b);
% initialize
path_spoofed = path_genuine;
for i = 1:n_chp
    StartPoint = path_genuine([1 3], t_chp_b(i));
    EndPoint = path_genuine([1 3], t_chp_e(i));
    X_middle_spoofed = point_mediator(StartPoint, EndPoint, 0.7);
    T_temp = t_chp_e(i) - t_chp_b(i)+1;
    spoofed_track_conditions.dim_given = {[1, 2, 3, 4], [1, 3], [1, 2, 3, 4]};
    spoofed_track_conditions.t_given = [1 floor(T_temp/2) T_temp];
    spoofed_track_conditions.X_given = {path_genuine(:, t_chp_b(i)), ...
        X_middle_spoofed, path_genuine(:, t_chp_e(i))};
    path_spoofed_section = conditional_sampling(mu_0, Cov_0, H, U, T_temp, ...
        spoofed_track_conditions);
    path_spoofed(:, t_chp_b(i):t_chp_e(i)) = path_spoofed_section;
end

% B.3 AIS and radar tracks in polar coordinates
y = [atan2(path_genuine(1, :), path_genuine(3, :));...
    hypot(path_genuine(1, :), path_genuine(3, :))];

z = [atan2(path_spoofed(1, :), path_spoofed(3, :));...
    hypot(path_spoofed(1, :), path_spoofed(3, :))];

% Get the radius and angle components only
Y = (y + theta*I_d(:, angle_comp))'; % This is Radar
Z = (z  + sqrt(sigmaNoise).*randn(2, T))'; % This is AIS