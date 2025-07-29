function cD = util_cD(area_i, AD_i, L, dx, beta, t, t0) 

[~, Gamma0, ~, ~, ~] = util_gamma0(area_i, L, dx, beta);

cv = std(area_i)/mean(area_i);

t_idx = t > t0;
x = 1./sqrt(t(t_idx));
y = AD_i(t_idx);
p = polyfit(x,y,1);

D0 = 2;
tortuosity_geo = mean(mean(area_i) ./ area_i);
d_inf = D0 / tortuosity_geo;
    
cD(1) = p(1); % based on diffusion measure
cD(2) = 2 * Gamma0 * sqrt(d_inf/pi); % based on geometrical measures
cD(3) = cv;


