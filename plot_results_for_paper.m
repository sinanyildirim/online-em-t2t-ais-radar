Assoc_err = performance_results.Assoc_err;
mean_Spoofdetectionerror_filt = performance_results.mean_Spoofdetectionerror_filt;
mean_Spoofdetectionerror_smooth = performance_results.mean_Spoofdetectionerror_smooth;
FP_overall_filt = performance_results.FP_overall_filt;
FN_overall_filt = performance_results.FN_overall_filt;
FP_overall_smooth = performance_results.FP_overall_smooth;
FN_overall_smooth = performance_results.FN_overall_smooth;
filt_spoof_prob = performance_results.filt_spoof_prob;
smooth_spoof_prob = performance_results.smooth_spoof_prob;
Negatives = performance_results.Negatives;
Positives = performance_results.Positives;
CM_filt = performance_results.CM_filt;
CM_smooth = performance_results.CM_smooth;


%% Plot spoofed and genuine tracks for spoofed targets
ind_spoofed_tracks = find(spoofed_tracks);
Ls = length(ind_spoofed_tracks);
nc = ceil(Ls/2);
for i = 1:Ls
    temp_ind = ind_spoofed_tracks(i);
    fc = fc + 1; fig = figure(fc);
    fig.Units = 'centimeters';
    fig.Position = [2 2 1.75 2.5];   % [left bottom width height]        
    scatter(path_genuine{temp_ind}(1, :), path_genuine{temp_ind}(3, :), '_'); hold on;
    scatter(path_spoofed{temp_ind}(1, :), path_spoofed{temp_ind}(3, :), '.k'); hold off;
    set(gca, 'xtick', [], 'ytick', []);
    title(sprintf('target no: %d', temp_ind), 'Interpreter', 'Latex');
    if i == Ls
        legend('genuine', 'spoofed');
    end
    exportgraphics(gca, sprintf('%s/spoofedtrack_%d.pdf', figuresavefolder, i));
    close(fig);
end

%% A. Plot the parameter estimates
% The transition probabilities
for i = 1:K
    for j = 1:K
        fc = fc + 1; fig = figure(fc);
        fig.Units = 'centimeters';
        fig.Position = [2 2 4.4 3];   % [left bottom width height]        
        fij = squeeze(Parameter_est.Fs(i, j, :));
        plot(fij);
        hold on;
        plot(F(i, j)*ones(1, length(fij)), 'r');
        hold off;
        title(sprintf('$F(%d, %d)$', i, j), 'interpreter', 'Latex');
        set(gca, 'xlim', [0, length(fij)], 'xtick',[0, length(fij)], 'ylim', [min(fij)-0.1, max(fij)+0.1]);
        if length(fij) == T
            xl = xlabel('$t$', 'Interpreter', 'Latex');
        else
            xl = xlabel('iteration', 'Interpreter', 'Latex');
        end
        xl.Units = 'normalized';
        pos = xl.Position;
        pos(2) = pos(2) + 0.2;   % small normalized shift
        xl.Position = pos;        
        exportgraphics(gca, sprintf('%s/F%d%d.pdf', figuresavefolder, i, j));
        close(fig);
    end
end
%%
% Theta
fc = fc + 1; fig = figure(fc);
fig.Units = 'centimeters';
fig.Position = [2 2 4.4 3];   % [left bottom width height] 
plot(Parameter_est.Thetas);
hold on;
plot(theta*ones(1, length(Parameter_est.Thetas)), 'r');
hold off;
title(sprintf('$\\theta$'), 'interpreter', 'Latex');
set(gca, 'xlim', [0, length(Parameter_est.Thetas)], 'xtick',[0, length(Parameter_est.Thetas)]);
if length(Parameter_est.Thetas) == T
    xl = xlabel('$t$', 'Interpreter', 'Latex');
else
    xl = xlabel('iteration', 'Interpreter', 'Latex');
end
xl.Units = 'normalized';
pos = xl.Position;
pos(2) = pos(2) + 0.2;   % small normalized shift
xl.Position = pos;        
exportgraphics(gca, sprintf('%s/theta.pdf', figuresavefolder));
close(fig);

%%
% error variance
sigma_titles = {{'$\sigma_{1, a}^{2}$', '$\sigma_{1, r}^{2}$'}, ...
    {'$\sigma_{2, a}^{2}$', '$\sigma_{2, r}^{2}$'}};
for i = 1:K
    for j = 1:d
        fc = fc + 1; fig = figure(fc);
        fig.Units = 'centimeters';
        fig.Position = [2 2 4.4 3];   % [left bottom width height] 
        sigmaij = squeeze(Parameter_est.Sigmas(i, j, :));
        plot(sigmaij);
        hold on;
        plot(sigmaNoise_MLE(j, i)*ones(1, length(sigmaij)), 'r');
        hold off;
        title(sigma_titles{j}{i}, 'interpreter', 'Latex');
        if length(Parameter_est.Thetas) == T
            xl = xlabel('$t$', 'Interpreter', 'Latex');
        else
            xl = xlabel('iteration', 'Interpreter', 'Latex');
        end
        % Get current position
        xl.Units = 'normalized';
        pos = xl.Position;
        pos(2) = pos(2) + 0.2;   % small normalized shift
        xl.Position = pos;        
        set(gca, 'xlim', [0, length(sigmaij)], 'xtick',[0, length(sigmaij)]);
        exportgraphics(gca, sprintf('%s/sigma%d%d.pdf', figuresavefolder, j, i));
        close(fig);
    end
end

%% B. Associations
fc = fc + 1; fig = figure(fc);
fig.Units = 'centimeters';
fig.Position = [2 2 4.4 3];   % [left bottom width height]
plot(1:T, Assoc_err);
set(gca, 'ylim', [-2/N, 1], 'xlim', [0, T], 'xtick', [0, T]);
ylabel('$\mathcal{E}^{a}_{t}$', 'Interpreter', 'Latex');
xl = xlabel('$t$', 'Interpreter', 'Latex');
% Get current position
xl.Units = 'normalized';
pos = xl.Position;
pos(2) = pos(2) + 0.2;   % small normalized shift
xl.Position = pos;
% title('Association error', 'Interpreter', 'Latex');
exportgraphics(gca, sprintf('%s/Associationerrror.pdf', figuresavefolder));
close(fig);

%% C. Spoof detection
fc = fc + 1; fig = figure(fc);
fig.Units = 'centimeters';
fig.Position = [2 2 4.4 3];   % [left bottom width height]        fij = squeeze(Parameter_est.Fs(i, j, :));        
plot(mean_Spoofdetectionerror_filt); 
% hold on; plot(mean_Spoofdetectionerror_smooth); hold off;
xl = xlabel('$t$', 'Interpreter', 'Latex');
% Get current position
xl.Units = 'normalized';
pos = xl.Position;
pos(2) = pos(2) + 0.05;   % small normalized shift
xl.Position = pos;
ylabel('$\mathcal{E}^{s}_{t}$', 'Interpreter', 'Latex');
set(gca, 'xlim', [0, T], 'xtick', [0, T]);
% title('avg. abs. error of filtering probab.', 'Interpreter', 'Latex');
% legend('filtering', 'smoothing');
exportgraphics(gca, sprintf('%s/error_spoofing_prob.pdf', figuresavefolder), 'ContentType', 'vector');
close(fig);

%% This generates a plot for each target. Relevant for only K = 2
fc = fc + 1; fig = figure(fc);
fig.Units = 'centimeters';
fig.Position = [2 2 8.8 6];   % [left bottom width height]
ylabels = cell(1, Ls);
for i = 1:Ls
    temp_ind = ind_spoofed_tracks(i);
    plot(filt_spoof_prob(:, temp_ind) + 2*i, 'b');
    hold on;
    % plot(smooth_spoof_prob(:, n) + 2*n, '-.r');
    plot(X{temp_ind}-1 + 2*i, '-.k');
    ylabels{i} = sprintf('no. %d', temp_ind);
end
hold off;
xlabel('$t$', 'Interpreter', 'Latex');
set(gca, 'ytick', 2*(1:Ls), 'yticklabel', ylabels, 'ylim', [0, 2*Ls+2], 'xlim', [0, T], 'xtick', [0, T], 'xticklabel', [0, T]);
% legend('filtering spoofing probability', 'smoothing spoofing probability', 'ground truth');
% legend('filtering probability', 'ground truth');
exportgraphics(gca, sprintf('%s/filtering_spoofing_all_targets.pdf', figuresavefolder), 'Resolution', 600, 'ContentType', 'vector');
close(fig);


%%

fc = fc + 1; fig = figure(fc);
fig.Units = 'centimeters';
fig.Position = [2 2 4.4 3];   % [left bottom width height]        
plot(thr_vec, FN_overall_filt, 'b')
hold on; plot(thr_vec, FP_overall_filt, '-.k');
% hold on; plot(FN_overall_smooth, 'r');
% hold on; plot(FP_overall_smooth, '-.r');
hold off;
legend('false negative rate', 'false positive rate');
xlabel('threshold value', 'Interpreter', 'Latex');
ylabel('error rate', 'Interpreter', 'Latex');
exportgraphics(gca, sprintf('%s/FP_overall.pdf', figuresavefolder), 'ContentType', 'vector');
close(fig);

fc = fc + 1; fig = figure(fc);
fig.Units = 'centimeters';
fig.Position = [2 2 4.4 3];   % [left bottom width height]        
plot(FP_overall_filt, 1-FN_overall_filt, '.-b');
% hold on;
% plot(FP_overall_smooth, FN_overall_smooth, 'r');
% hold off;
xlabel('FPR', 'Interpreter', 'Latex');
ylabel('TPR', 'Interpreter', 'Latex');
exportgraphics(gca, sprintf('%s/ROC.pdf', figuresavefolder), 'ContentType', 'vector');
close(fig);

% legend('filtering', 'smoothing');
%%
thr_plot = 0.5;
i_plot = find(thr_vec == thr_plot);
fc = fc + 1; fig = figure(fc);
fig.Units = 'centimeters';
fig.Position = [2 2 4.4 3];   % [left bottom width height]
plot(Negatives, '.-b'); hold on; plot(squeeze(CM_filt{i_plot}(1, 2, :)), 'k'); hold off;
xl = xlabel('$t$', 'Interpreter', 'Latex'); legend('negatives', 'false detection');
xl.Units = 'normalized';
pos = xl.Position;
pos(2) = pos(2) + 0.05;   % small normalized shift
xl.Position = pos;
set(gca, 'xlim', [0, T], 'xtick', [0, T]);
exportgraphics(gca, sprintf('%s/FP_fixed_thr.pdf', figuresavefolder), 'ContentType', 'vector');
close(fig);

fc = fc + 1; fig = figure(fc);
fig.Units = 'centimeters';
fig.Position = [2 2 4.4 3];   % [left bottom width height]
plot(Positives, '.-b'); hold on; plot(squeeze(CM_filt{i_plot}(2, 1, :)), 'k'); hold off;
xl = xlabel('$t$', 'Interpreter', 'Latex'); legend('positives', 'misdetection');
xl.Units = 'normalized';
pos = xl.Position;
pos(2) = pos(2) + 0.05;   % small normalized shift
xl.Position = pos;
set(gca, 'xlim', [0, T], 'xtick', [0, T]);
exportgraphics(gca, sprintf('%s/FN_fixed_thr.pdf', figuresavefolder), 'ContentType', 'vector');
close(fig);