function [performance_results] = eval_performance(X, Assoc_est, FiltX, SmoothX, thr_vec)

%% Compute the results
Best_Assoc = Assoc_est.Best_Assoc;
w_assoc = Assoc_est.w_assoc;

[B, T] = size(w_assoc);
N = length(FiltX);
K = size(FiltX{1}, 2);

% Calculate the average association error vs time
Error_mtx = zeros(B, T);
for t = 1:T
    Error_mtx(:, t) = arrayfun(@(b) sum(Best_Assoc{t}{b} ~= (1:N)'), 1:B);
end
Assoc_err = sum(w_assoc.*Error_mtx, 1)/N;

%% C. Spoof detection
% This is plotting the average asolute error. This is relevant for only K = 2
filt_spoof_prob = zeros(T, N);
smooth_spoof_prob = zeros(T, N);
Spoofdetectionerror_mtx_filt = zeros(T, N);
Spoofdetectionerror_mtx_smooth = zeros(T, N);
for n = 1:N
    filt_spoof_prob(:, n) = FiltX{n}*[0; 1];
    smooth_spoof_prob(:, n) = SmoothX{n}*[0; 1];
    Spoofdetectionerror_mtx_filt(:, n) = abs(filt_spoof_prob(:, n) - (X{n}==2)');
    Spoofdetectionerror_mtx_smooth(:, n) = abs(smooth_spoof_prob(:, n) - (X{n}==2)');
end
mean_Spoofdetectionerror_filt = mean(Spoofdetectionerror_mtx_filt, 2);
mean_Spoofdetectionerror_smooth = mean(Spoofdetectionerror_mtx_smooth, 2);

%% false detection and misdetection across various thresholds

L_thr = length(thr_vec);
FP_overall_filt = zeros(1, L_thr);
FN_overall_filt = zeros(1, L_thr);
FP_overall_smooth = zeros(1, L_thr);
FN_overall_smooth = zeros(1, L_thr);

CM_filt = cell(1, L_thr);
CM_smooth = cell(1, L_thr);
X_mtx = cell2mat(X')'; % T x N
Negatives = sum(X_mtx == 1, 2);
Positives = sum(X_mtx == 2, 2);

for i = 1:L_thr
    %%% Compute the confusion matrix and the F1-score
    filt_estLabel = 1*(filt_spoof_prob < thr_vec(i)) + 2*(filt_spoof_prob >= thr_vec(i));
    smooth_estLabel = 1*(smooth_spoof_prob < thr_vec(i)) + 2*(smooth_spoof_prob >= thr_vec(i));

    CM_filt{i} = zeros(K,K,T);
    CM_smooth{i} = zeros(K,K,T);
    
    for k = 1:K
        for l = 1:K
            CM_filt{i}(k, l, :) = sum(X_mtx == k & filt_estLabel == l, 2);
            CM_smooth{i}(k, l, :) = sum(X_mtx == k & smooth_estLabel == l, 2);
        end
    end
    FP_rate_filt = squeeze(CM_filt{i}(1, 2, :))./Negatives;
    FN_rate_filt = squeeze(CM_filt{i}(2, 1, :))./Positives;
    FP_rate_filt(isnan(FP_rate_filt)) = 0;
    FN_rate_filt(isnan(FN_rate_filt)) = 0;

    FP_overall_filt(i) = mean(FP_rate_filt);
    FN_overall_filt(i) = mean(FN_rate_filt);

    FP_rate_smooth = squeeze(CM_smooth{i}(1, 2, :))./Negatives;
    FN_rate_smooth = squeeze(CM_smooth{i}(2, 1, :))./Positives;
    FP_rate_smooth(isnan(FP_rate_smooth )) = 0;
    FN_rate_smooth(isnan(FN_rate_smooth)) = 0;

    FP_overall_smooth(i) = mean(FP_rate_smooth);
    FN_overall_smooth(i) = mean(FN_rate_smooth);
end

performance_results.Assoc_err = Assoc_err;
performance_results.mean_Spoofdetectionerror_filt = mean_Spoofdetectionerror_filt;
performance_results.mean_Spoofdetectionerror_smooth = mean_Spoofdetectionerror_smooth;
performance_results.FP_overall_filt = FP_overall_filt;
performance_results.FN_overall_filt = FN_overall_filt;
performance_results.FP_overall_smooth = FP_overall_smooth;
performance_results.FN_overall_smooth = FN_overall_smooth;
performance_results.filt_spoof_prob = filt_spoof_prob;
performance_results.smooth_spoof_prob = smooth_spoof_prob;
performance_results.Negatives = Negatives;
performance_results.Positives = Positives;
performance_results.CM_filt = CM_filt;
performance_results.CM_smooth= CM_smooth;


