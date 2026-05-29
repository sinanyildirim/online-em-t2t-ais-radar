function y = log_sum_exp(x)

y = log(sum(exp(x - max(x)))) + max(x); 