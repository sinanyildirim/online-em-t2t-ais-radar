function [Parameter_est, Assoc_est, FiltX, SmoothX] = T2TA_onEM(Z, Y, ...
    K, F, Sigma, theta, angle_comp, B, a_EM, burnin_EM, L_smooth)

% [Parameter_est, Assoc_est, FiltX, SmoothX] = T2TA_onEM(Z, Y, K, F, ...
% Sigma, theta, angle_comp, B, a_EM, burnin_EM, L_smooth)
% 
% This function performs Online EM for parameter estimation and 
% online spoofing detection with track-to-track assocaition.
% 
% INPUTS
% 
% Z: AIS tracks
% Y: radar tracks
% K: number of states for X_t: the spoofing/genuine state
% F: initial transition matrix for X_t
% Sigma: initial covariance matrix for the difference between radar and AIS
% signals
% theta: initial angular bias
% angle_comp: the index number for the angular component
% B: number of best associations to be considered
% a_EM: the power parameter for the stochastic approximation in online EM
% burnin_EM: burn-in time for online EM - prior to which no EM updates
% L_smooth: the lag of real-time (fixed-lag) backward smoothing.
%
% OUTPUTS
% 
% Parameter_est: Online EM estimates over time
% Assoc_est: Association estimates
% FiltX: Filtering probabilities for X_t
% SmoothX: Fixed-lag (L) smoothing probabilities for X_t

% get the sizes
[T, d] = size(Y{1});
M = length(Y);
N = length(Z);

% These matrices will be needed
I_K = eye(K);
I_d = eye(d);
theta_vec = I_d(angle_comp, :)*theta;

% initialize the assignment matrix
L = zeros(M, N);
% initialize the forward filtering and prediction probabilities
% T = number of time steps; K = Number of modes of x
Alpha = repmat({zeros(T, K)}, M, N);
Beta = repmat({zeros(T, K)}, M, N);
Gamma = repmat({zeros(T, K)}, M, N);

% initialize the expected statistics of EM
BS = repmat({{zeros(1, K), zeros(1, K), zeros(1, K), zeros(d, K), zeros(K)}}, M, N); 

% initialize the intermediate functions
T_f = repmat({{repmat({zeros(1, K)}, 1, K), repmat({zeros(1, K)}, 1, K), ...
    repmat({zeros(1, K)}, 1, K), repmat({zeros(1, K, d)}, 1, K), repmat({zeros(1, K, K)}, 1, K)}}, M, N);

Thetas = zeros(1, T);
Fs = zeros(K, K, T);
Sigmas = zeros(K, d, T);

Best_Assoc = cell(1, T);
w_assoc = zeros(B, T);

FiltX = repmat({zeros(T,K)},1,M);
SmoothX = repmat({zeros(T,K)},1,M);

for t = 1:T
    if mod(t, 100) == 0
        disp(t);
    end
    gamma_t = t^(-a_EM); % coefficient for the stochastic approximation
    for m = 1:M
        for n = 1:N
            % Forward filtering
            if t == 1
                alpha_prev = [1 zeros(1, K-1)];
            else
                alpha_prev = Alpha{m, n}(t-1, :);
            end
            [Alpha{m, n}(t, :), Beta{m, n}(t, :), L(m, n)] ...
                = ffpfu(alpha_prev, L(m, n), Y{n}(t, :), Z{m}(t, :), F, theta_vec, Sigma); 
            
            % Backward smoothing
            Gamma{m, n}(t, :) = Alpha{m, n}(t, :);
            for tau = t-1:-1:max(t-L_smooth, 1)
                temp_1 = Gamma{m, n}(tau+1, :)./Beta{m, n}(tau+1, :);
                temp_2 = Alpha{m, n}(tau, :)'.*Fs(:, :, tau);
                Gamma{m, n}(tau, :) = sum(temp_1.*temp_2, 2)';
            end

            % prediction distribution
            beta_vec = Beta{m, n}(t, :);
            alpha_vec = Alpha{m, n}(t, :);
            Alpha_F_mtx = alpha_prev'.*F;

            if burnin_EM < T % If burnin_EM = T, no statistic is needed
                for i = 1:K
                    % Calculate the additive statistics S^1,... ,S^5
                    % I_K : matrix conta
                    s{1} = repmat(I_K(:, i), 1, K);
                    s{2} = repmat(I_K(i, :), K, 1);
                    s{3} = (Y{n}(t, angle_comp) - Z{m}(t, angle_comp))*s{2};
                    s{4} = tensorprod(s{2}, ((Y{n}(t, :) - theta_vec -  Z{m}(t, :) ).^2)');
                    s{5} = squeeze(tensorprod(I_K(:, i), I_K)); %=S^2
                    for is = 1:5
                        % Update the intermediate function of forward
                        % smoothing
                        Tpluss_mtx = (1 - gamma_t)*pagetranspose(T_f{m, n}{is}{i}) + gamma_t*s{is};
                        numer_mtx = Tpluss_mtx.*Alpha_F_mtx;
                        numer_mtx_sum = sum(numer_mtx, 1); %=values of the sum on j for different values of i                    
                        T_f{m, n}{is}{i} = numer_mtx_sum./beta_vec;
                        % calculate the expectations
                        BS{m, n}{is}(:, i) = squeeze(pagemtimes(T_f{m, n}{is}{i}, alpha_vec'));
                    end
                end
            end
        end
    end

    % B-best with assignkbest
    %L is the log likelihood of the associations used as a cost for the associations
    [Us, ~, ~, costs]= assignkbest(-L, 2*M*(max(-L(:))-min(-L(:))), B); 
    % This code returns B assignements possibilities, the first is the one
    % with the smallest total cost
    % (approximate) posterior probabilities (positive for the B best assignments)
    ls = -costs;
    % posterior probabilities of the B best assignments in terms of
    % likelihood
    post_prob_assgn = exp(ls - log_sum_exp(ls));

    BS_avg = {zeros(1, K), zeros(1, K), zeros(1, K), zeros(d, K), zeros(K, K)};
    for b = 1:B
        % the assignment vector for the b'th solution
        inds_to_get = sub2ind([M, N], Us{b}(:, 1), Us{b}(:, 2));
        for is = 1:5
            for js = 1:M
                BS_avg{is} = BS_avg{is} + post_prob_assgn(b)*BS{inds_to_get(js)}{is};
            end
        end
    end
    
    estAssoc = cell(1,B);
    for b = 1:B
        estAssoc{b} = double(Us{b}(:,2));
    end

    % M-step (after sufficiently many time steps so that the statistics
    % have 'matured').
    if t > burnin_EM
        Sigma = (BS_avg{4}./BS_avg{2})';
        F = BS_avg{5}./BS_avg{1}';
        theta = (BS_avg{3}*Sigma(:, angle_comp))/(BS_avg{2}*Sigma(:, angle_comp));
        theta_vec = I_d(angle_comp, :)*theta;

        % Postprocess F
        [~, temp_b] = sort(Sigma(:, 2));
        for i = 1:K
            for j = 1:K
                if abs(temp_b(i)-temp_b(j))>1
                    F(i, j) = 0;
                end
            end
        end
        F = F./sum(F, 2);
    end
    
    % Store the estimates
    Thetas(t) = theta;
    Fs(:, :, t) = F;
    Sigmas(:, :, t) = Sigma;
    Best_Assoc{t} = estAssoc;
    w_assoc(:, t) = post_prob_assgn;
    
    for m = 1:M
        % Filtering
        FiltX{m}(t,:) = zeros(1,K);
        SmoothX{m}(t:-1:max(t-L_smooth, 1),:) = 0;
        for b = 1:B
            n = Us{b}(Us{b}(:, 1) == m, 2);
            FiltX{m}(t,:) = FiltX{m}(t,:) + Alpha{m,n}(t,:)*post_prob_assgn(b);
            % Smoothing 
            for tau = t:-1:max(t-L_smooth, 1)
                SmoothX{m}(tau,:) = SmoothX{m}(tau,:) + Gamma{m, n}(tau,:)*post_prob_assgn(b);
            end
        end
    end
end

%% Store the outputs
Parameter_est.Thetas = Thetas;
Parameter_est.Fs= Fs;
Parameter_est.Sigmas= Sigmas;
Assoc_est.Best_Assoc = Best_Assoc;
Assoc_est.w_assoc = w_assoc;