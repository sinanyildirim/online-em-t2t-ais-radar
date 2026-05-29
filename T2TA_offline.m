function [Parameter_est] = T2TA_offline(Z, Y, K, F, Sigma, theta, angle_comp, B, EM_iter)

% [Parameter_est] = T2TA_offline(Z, Y, K, F, Sigma, theta, angle_comp, B, EM_iter)

% get the sizes
[T, d] = size(Y{1});
M = length(Y);
N = length(Z);

% These matrices will be needed
I_K = eye(K);
I_d = eye(d);
theta_vec = I_d(angle_comp, :)*theta;

Thetas = zeros(1, EM_iter);
Fs = zeros(K, K, EM_iter);
Sigmas = zeros(K, d, EM_iter);

S1_ready = cell(1, K);
S2_ready = cell(1, K);
S5_ready = cell(1, K);
for k = 1:K
    S1_ready{k} = repmat(I_K(:, k), 1, K);
    S2_ready{k} = repmat(I_K(k, :), K, 1);
    S5_ready{k} = squeeze(tensorprod(I_K(:, k), I_K));
end


for iter = 1:EM_iter
    disp(iter);

    % First, run the algorithm for every pari of AIS and radar tracks and
    % just store their log-likelihoods
    L = zeros(M, N);
    for m = 1:M
        for n = 1:N
            disp([m, n])
            l_mn = 0;
            for t = 1:T
                ztm = Z{m}(t, :);
                ytn = Y{n}(t, :);
                % Forward filtering
                if t == 1
                    alpha_prev = [1 zeros(1, K-1)]; %Modification to check
                else
                    alpha_prev = alpha_vec; %alpha_(t-1)
                end
                % MaJ de alpha et beta par l'algo de filtrage
                [alpha_vec, ~, l_mn] = ffpfu(alpha_prev, l_mn, ytn, ztm, F, theta_vec, Sigma);
            end
            L(m, n) = l_mn;
        end
    end

    % Now, find the best B associations
    [Us, ~, ~, costs]= assignkbest(-L, 2*M*(max(-L(:))-min(-L(:))), B);
    ls = -costs;
    post_prob_assgn = exp(ls - log_sum_exp(ls));

    % For each association, run the forward filtering again for the pairs
    % corresponding to that association

    % initialize the expected statistics of EM
    BS_avg = {zeros(1, K), zeros(1, K), zeros(1, K), zeros(d, K), zeros(K)};

    for b = 1:B

        m_inds = Us{b}(:, 1);
        n_inds = Us{b}(:, 2);
        % the assignment vector for the b'th solution
        for m0 = 1:M
            m = m_inds(m0);
            n = n_inds(m0);

            % initialize the intermediate functions
            T_f = {repmat({zeros(1, K)}, 1, K), repmat({zeros(1, K)}, 1, K), ...
                repmat({zeros(1, K)}, 1, K), repmat({zeros(1, K, d)}, 1, K), repmat({zeros(1, K, K)}, 1, K)};


            % Now run the forward filtering according to this pair
            for t = 1:T
                ztm = Z{m}(t, :);
                ytn = Y{n}(t, :);
                % Forward filtering
                if t == 1
                    alpha_prev = [1 zeros(1, K-1)]; %Modification to check
                else
                    alpha_prev = alpha_vec; %alpha_(t-1)
                end
                % MaJ de alpha et beta par l'algo de filtrage
                [alpha_vec, beta_vec, l_mn] = ffpfu(alpha_prev, l_mn, ytn, ztm, F, theta_vec, Sigma);

                % prediction distribution
                Alpha_F_mtx = alpha_prev'.*F;

                for i = 1:K
                    % Calculate the additive statistics S^1,... ,S^5
                    % I_K : matrix conta
                    s{1} = S1_ready{i};
                    s{2} = S2_ready{i};%=S^4
                    s{3} = (ytn(angle_comp) - ztm(angle_comp))*s{2}; %=S^5
                    s{4} = tensorprod(s{2}, ((ytn - theta_vec -  ztm).^2)'); %=S^1
                    s{5} = S5_ready{i};
                    for is = 1:5
                        % Update T
                        Tpluss_mtx = pagetranspose(T_f{is}{i}) + s{is};
                        numer_mtx = Tpluss_mtx.*Alpha_F_mtx;
                        numer_mtx_sum = sum(numer_mtx, 1); % = values of the sum on j for different values of i
                        T_f{is}{i} = numer_mtx_sum./beta_vec;
                    end
                end
            end

            for is = 1:5
                for i = 1:K
                    SS = squeeze(pagemtimes(T_f{is}{i}, alpha_vec'));
                    BS_avg{is}(:, i) = BS_avg{is}(:, i) + post_prob_assgn(b)*SS;
                end
            end
        end
    end

    % M-step (after sufficiently many time steps so that the statistics
    % have 'matured').
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

    % Store the estimates
    Thetas(iter) = theta;
    Fs(:, :, iter) = F;
    Sigmas(:, :, iter) = Sigma;
end

%% Store the outputs
Parameter_est.Thetas = Thetas;
Parameter_est.Fs = Fs;
Parameter_est.Sigmas = Sigmas;