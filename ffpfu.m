function [alpha_tp1, beta_tp1, log_lkl_tp1] = ffpfu(alpha_t, log_lkl_t, y_tp1, z_tp1, F, theta_vec, Sigma)

% [alpha_tp1, beta_tp1, log_lkl_tp1] = ffpfu(alpha_t, log_lkl_t, y_tp1, z_tp1, F, theta_vec, Sigma)
% 
% Performs forward filtering for the T2TA problem.

beta_tp1 = alpha_t*F; % matrix multiplication

%log de la densité de probabilité g
log_g = sum(-0.5*log(2*pi*Sigma) - 0.5*((y_tp1 - z_tp1 - theta_vec).^2)./Sigma, 2); 

%Calcule le ln du numérateur dans la formule de filtrage
alpha_unnorm_vec = log(beta_tp1) + log_g'; 
%Calcul le ln de la somme située dans le denominateur dans la formule de filtrage
log_inc_lkl = log_sum_exp(alpha_unnorm_vec); 
%Calcule le ln de alpha à partir de la formule de filtrage
log_alpha_tp1 = alpha_unnorm_vec - log_inc_lkl; 
alpha_tp1 = exp(log_alpha_tp1);

%Mise à jour de la log-vraisemblance marginale d'après la formule du rapport
log_lkl_tp1 = log_lkl_t + log_inc_lkl; 