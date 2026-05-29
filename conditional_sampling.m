function X = conditional_sampling(mu_0, Cov_0, F, U, T, track_conditions)

% Generates a Linear and Markov chain that is conditional on 
% 
% X(t_given(i)) = X_given(i)
% 
% given time points indicated in t_given


t_given = track_conditions.t_given;
dim_given = track_conditions.dim_given;
X_given = track_conditions.X_given;

obs_vec = zeros(1, T);
obs_vec(t_given) = 1;

d = length(mu_0);

%% Forward pass:
mu_filts = zeros(d, T);
Sigma_filts = zeros(d, d, T);
c = 0;

for t = 1:T
    if t == 1
        mu_pred = mu_0;
        Sigma_pred = Cov_0;
    else
        mu_pred = F*mu_filt;
        Sigma_pred = F*Sigma_filt*F' + U;
        Sigma_pred = (Sigma_pred + Sigma_pred')/2;
    end
    % Filtering
    if obs_vec(t) == 1
        c = c + 1;
        d_given  = dim_given{c};
        d_not_given = setdiff(1:d, d_given);
        l_d = length(d_given);
        G = zeros(l_d, d);
        for i = 1:l_d
            G(l_d, d_given(i)) = 1;
        end
        
        Sigma_filt = zeros(d);
        mu_filt = zeros(d, 1);
        Sigma_filt(d_given, d_given) = zeros(l_d);
        mu_filt(d_given) = X_given{c};
        Sigma_filt(d_given, d_given) = zeros(l_d);
        a =  X_given{c};
        mu2 = mu_pred(d_given);
        mu1 = mu_pred(d_not_given);
        
        Sigma11 = Sigma_pred(d_not_given, d_not_given);
        Sigma12 = Sigma_pred(d_not_given, d_given);
        Sigma21 = Sigma_pred(d_given, d_not_given);
        Sigma22 = Sigma_pred(d_given, d_given);
        
        mu_filt(d_not_given) = mu1 + Sigma12*(Sigma22\(a - mu2));
        Sigma_filt(d_not_given, d_not_given) = Sigma11 - Sigma12*(Sigma22\Sigma21);
    else
        mu_filt  = mu_pred;
        Sigma_filt = Sigma_pred;
    end

    Sigma_filt = (Sigma_filt + Sigma_filt')/2;
    % Store the filtering moments
    mu_filts(:, t) = mu_filt;
    Sigma_filts(:, :, t) = Sigma_filt;    
end

%% Backward pass
for t = T:-1:1
    % sample X
    if t == T
        mu_smooth = mu_filts(:, T);
        Sigma_smooth = Sigma_filts(:, :, T);
        Sigma_smooth = (Sigma_smooth + Sigma_smooth')/2;
    else
        mu_filt = mu_filts(:, t);
        Sigma_filt = Sigma_filts(:, :, t);
        
        z = X(:, t+1) - F*mu_filt;
        S = F*Sigma_filt*F' + U;
        K = Sigma_filt*(F'/S);
        mu_smooth = mu_filt + K*z; %Updated (a posteriori) state estimate
        Sigma_smooth = (eye(d) - K*F)*Sigma_filt; %Updated (a posteriori) estimate covariance
        Sigma_smooth = (Sigma_smooth + Sigma_smooth')/2;
    end
    X(:, t) = mvnrnd(mu_smooth, Sigma_smooth);
end
