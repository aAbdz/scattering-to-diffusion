function w_std = util_weighted_std(x,w)

N = length(w);
x_bar_w = sum(x .* w);

num = sum(w.* (x - x_bar_w).^2);
den = (N-1) .* sum(w) ./ N;

w_std = sqrt(num / den);
