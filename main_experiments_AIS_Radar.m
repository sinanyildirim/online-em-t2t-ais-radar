% Reproduce the synthetic AIS/radar experiment from the paper
% "Online EM for T2T Association and Spoofing Detection Using AIS and Radar Data".
%
% The script generates paired radar/AIS tracks with intermittent AIS spoofing,
% runs the online EM track-to-track association and spoofing-detection method,
% compares it with an offline EM baseline, and writes the resulting workspace
% and paper-style figures under the local results/ folder.

clear; clc; close all; fc = 0;
assert_prerequisites;

rng(1);
saveoutputs = true;
plotresults = true;
plotresultstogether = true;

resultsfolder = fullfile(pwd, 'results');
figuresavefolder = fullfile(resultsfolder, 'figures');
outputsavefolder = fullfile(resultsfolder, 'outputdata');
if ~exist(figuresavefolder, 'dir')
    mkdir(figuresavefolder);
end
if ~exist(outputsavefolder, 'dir')
    mkdir(outputsavefolder);
end
T = 5000; % number of time steps
d = 2; % dimension of the observations
angle_comp = 1; % the component where the angle information is stored
theta = 0.1; % angular bias
N = 40; % number of tracks for AIS and Radar
K = 2; % number of modes

F = [0.999 0.001;
    0.005 0.995];
I_d = eye(2);
sigmaNoise = 0.25*rand(2, 1);

%Sigma = (0.5*(1:K)').^2*rand(1, d);
H = [[1, 0.1, 0, 0]; [0, 0.99, 0, 0]; [0, 0, 1, 0.1]; [0, 0, 0, 0.99]];
U = diag([0.00001, 0.001, 0.00001, 0.001]);

Y0 = [0, 0.1, 0, 0.1]'; % tracks start with this value
YT = [100, 0.1, 100, 0.1]'; % tracks end with this value
mu_0 = Y0;
Cov_0 = 0.1*eye(4);

T_mid = floor(T/2);
t_given_genuine = [1 T_mid T];
dim_given_genuine = {[1, 2, 3, 4], [1, 3], [1, 2, 3, 4]};
Y_middle_genuine = 0.5.*([Y0(1) Y0(3)]+[YT(1) YT(3)])';
Y_given_genuine = {Y0, Y_middle_genuine, YT};

X = cell(1, N);
Y = cell(1, N);
Z = cell(1, N);
path_genuine = cell(1, N);
path_spoofed = cell(1, N);
spoofed_tracks = zeros(1, N);

spoof_parameters.spoof_create_mode = 2;
spoof_parameters.F = F;
spoof_parameters.spoof_time_coeff_vec = [1, 1, 1];
spoof_parameters.probSpoofedTraj = 0.2;
spoof_parameters.min_spoof_duration = 500;

track_conditions.t_given = t_given_genuine;
track_conditions.dim_given = dim_given_genuine;
track_conditions.X_given = Y_given_genuine;

track_dynamics.mu0 = mu_0;
track_dynamics.Cov0 = Cov_0;
track_dynamics.H = H;
track_dynamics.U = U;
track_dynamics.angle_comp = angle_comp;
track_dynamics.theta = theta;
track_dynamics.sigmaNoise = sigmaNoise;

for n = 1:N
    % Generate the (AIS, Radar) pair
    [X{n}, Y{n}, Z{n}, path_genuine{n}, path_spoofed{n}] = create_track_pairs(T, track_dynamics, spoof_parameters, track_conditions);
    spoofed_tracks(n) = sum(X{n} == 2) > 0;
end

%% Estiation of the parameters according using the complete data
% MLE_complete_data()
Y_concat = cell2mat(Y')';
Z_concat = cell2mat(Z')';
X_concat = cell2mat(X);

E_concat = Z_concat - (Y_concat - theta*I_d(:, angle_comp));
sigmaNoise_MLE = zeros(d, K);
for k = 1:K
    sigmaNoise_MLE(:, k) = var(E_concat(:, (X_concat == k)), [], 2);
end
% Estimate the F matrix
S = zeros(K); % Initialize the matrix of counted transitions
for n = 1:N
    for t = 2:T
        S(X{n}(t-1), X{n}(t)) = S(X{n}(t-1), X{n}(t)) + 1;
    end
end
% Estimate the transition probabilities
F_MLE = S ./ sum(S, 2);

%% Estimation, tracking, and association from partially observed data
B = 5; % the number of best assignments to be considered
a_EM = 0.65; % power for the coefficient of stochastic approximation
burnin_EM = 0.01*T; % burn-in time for EM
Sigma0 = [0.2*ones(K, 1) (10.^(0:K-1))']; % Initial noise parameter first row is for angle
theta0 = 0.2; % initial angular bias
F0 = (0.1/(K-1))*ones(K)+(0.9-0.1/(K-1))*eye(K); % initial transition probabilities
L_smooth = 0;

% Run the algorithm
[Parameter_est_online, Assoc_est_online, FiltX_online, SmoothX_online] = T2TA_onEM(Z, Y, K, F0, Sigma0, theta0, angle_comp, B, a_EM, burnin_EM, L_smooth);

EM_iter = 10;
[Parameter_est_offline] = T2TA_offEM(Z, Y, K, F, Sigma0, theta0, angle_comp, B, EM_iter);
F_est = Parameter_est_offline.Fs(:, :, end);
Sigma_est = Parameter_est_offline.Sigmas(:, :, end);
theta_est = Parameter_est_offline.Thetas(end);
[~, Assoc_est_offline, FiltX_offline, SmoothX_offline] = T2TA_onEM(Z, Y, K, F_est, Sigma_est, theta_est, angle_comp, B, a_EM, T, L_smooth);

%%
thr_vec = 0.0:0.001:1.00;
[performance_results_online] = eval_performance(X, Assoc_est_online, FiltX_online, SmoothX_online, thr_vec);
[performance_results_offline] = eval_performance(X, Assoc_est_offline, FiltX_offline, SmoothX_offline, thr_vec);

if saveoutputs == 1
    save(fullfile(outputsavefolder, 'T2TResults.mat'));
end

%%
if plotresults == 1
    Parameter_est = Parameter_est_online;
    Assoc_est = Assoc_est_online;
    FiltX = FiltX_online;
    SmoothX = SmoothX_online;
    performance_results = performance_results_online;
    plot_results_for_paper;
end

%%
if plotresultstogether == 1
   plot_online_offline_results; 
end
