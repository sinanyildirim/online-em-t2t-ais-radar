Assoc_err_online = performance_results_online.Assoc_err;
mean_Spoofdetectionerror_filt_online = performance_results_online.mean_Spoofdetectionerror_filt;
mean_Spoofdetectionerror_smooth_online = performance_results_online.mean_Spoofdetectionerror_smooth;
FP_overall_filt_online = performance_results_online.FP_overall_filt;
FN_overall_filt_online = performance_results_online.FN_overall_filt;
FP_overall_smooth_online = performance_results_online.FP_overall_smooth;
FN_overall_smooth_online = performance_results_online.FN_overall_smooth;
filt_spoof_prob_online = performance_results_online.filt_spoof_prob;
smooth_spoof_prob_online = performance_results_online.smooth_spoof_prob;

Assoc_err_offline = performance_results_offline.Assoc_err;
mean_Spoofdetectionerror_filt_offline = performance_results_offline.mean_Spoofdetectionerror_filt;
mean_Spoofdetectionerror_smooth_offline = performance_results_offline.mean_Spoofdetectionerror_smooth;
FP_overall_filt_offline = performance_results_offline.FP_overall_filt;
FN_overall_filt_offline = performance_results_offline.FN_overall_filt;
FP_overall_smooth_offline = performance_results_offline.FP_overall_smooth;
FN_overall_smooth_offline = performance_results_offline.FN_overall_smooth;
filt_spoof_prob_offline = performance_results_offline.filt_spoof_prob;
smooth_spoof_prob_offline = performance_results_offline.smooth_spoof_prob;

fc = fc + 1; fig = figure(fc);
fig.Units = 'centimeters';
fig.Position = [2 2 4.4 3];   % [left bottom width height]
plot(thr_vec, FN_overall_filt_online, 'b')
hold on; plot(thr_vec, FP_overall_filt_online, '-k');
hold on; plot(thr_vec, FN_overall_filt_offline, 'r');
hold on; plot(thr_vec, FP_overall_filt_offline, '-.g');
hold off;
legend('false negative rate', 'false positive rate');
xlabel('threshold value', 'Interpreter', 'Latex');
ylabel('error rate', 'Interpreter', 'Latex');
exportgraphics(gca, sprintf('%s/FP_overall_online_offline.pdf', figuresavefolder), 'ContentType', 'vector');
close(fig);
%
fc = fc + 1; fig = figure(fc);
fig.Units = 'centimeters';
fig.Position = [2 2 4.4 3];   % [left bottom width height]
plot(FP_overall_filt_online, 1-FN_overall_filt_online, '.-k', 'markersize', 0.5);
hold on;
plot(FP_overall_filt_offline, 1-FN_overall_filt_offline, '.-b', 'markersize', 0.5);
hold off;
set(gca, 'xlim', [0, 0.1]);
xlabel('FPR', 'Interpreter', 'Latex');
ylabel('TPR', 'Interpreter', 'Latex');
legend('online', 'offline', 'Location', 'southeast');
exportgraphics(gca, sprintf('%s/ROC_online_offline.pdf', figuresavefolder), 'ContentType', 'vector');
close(fig);

%% This generates a plot for each target. Relevant for only K = 2
fc = fc + 1; fig = figure(fc);
fig.Units = 'centimeters';
fig.Position = [2 2 8.8 6];   % [left bottom width height]
ylabels = cell(1, Ls);
for i = 1:Ls
    temp_ind = ind_spoofed_tracks(i);
    plot(X{temp_ind}-1 + 2*i, '-.r');
    hold on;    
    plot(filt_spoof_prob_offline(:, temp_ind) + 2*i, 'k');
    plot(filt_spoof_prob_online(:, temp_ind) + 2*i, 'b');
    % plot(smooth_spoof_prob(:, n) + 2*n, '-.r');
    ylabels{i} = sprintf('no. %d', temp_ind);
end
hold off;
xl = xlabel('$t$', 'Interpreter', 'Latex');
% Get current position
xl.Units = 'normalized';
pos = xl.Position;
pos(2) = pos(2) + 0.1;   % small normalized shift
xl.Position = pos;
set(gca, 'ytick', 2*(1:Ls), 'yticklabel', ylabels, 'ylim', [0, 2*Ls+2], 'xlim', [0, T], 'xtick', [0, T], 'xticklabel', [0, T]);
% legend('filtering spoofing probability', 'smoothing spoofing probability', 'ground truth');
% legend('filtering probability', 'ground truth');
exportgraphics(gca, sprintf('%s/filtering_spoofing_all_targets_online_offline.pdf', figuresavefolder), 'Resolution', 600, 'ContentType', 'vector');
close(fig);

%% C. Spoof detection
fc = fc + 1; fig = figure(fc);
fig.Units = 'centimeters';
fig.Position = [2 2 4.4 3];   % [left bottom width height]        fij = squeeze(Parameter_est.Fs(i, j, :));        
plot(mean_Spoofdetectionerror_filt_offline, 'k'); 
hold on;
plot(mean_Spoofdetectionerror_filt_online, 'b');
hold off;
% hold on; plot(mean_Spoofdetectionerror_smooth); hold off;
xl = xlabel('$t$', 'Interpreter', 'Latex');
% Get current position
xl.Units = 'normalized';
pos = xl.Position;
pos(2) = pos(2) + 0.1;   % small normalized shift
xl.Position = pos;
ylabel('$\mathcal{E}^{s}_{t}$', 'Interpreter', 'Latex');
set(gca, 'xlim', [0, T], 'xtick', [0, T]);
% title('avg. abs. error of filtering probab.', 'Interpreter', 'Latex');
% legend('online', 'offline', 'Location', 'northeast');
exportgraphics(gca, sprintf('%s/error_spoofing_prob_online_offline.pdf', figuresavefolder), 'ContentType', 'vector');
close(fig);