function [Gamma, Gamma0, ind_max, xval_g, yval_g] = util_gamma0(area_i, L, dx, beta) 

dk = 2*pi/L;
kmax = pi/dx;
kx = 0:dk:2*kmax-dk;

alpha = area_i / mean(area_i);
alpha = log(alpha);

alpha_fft = fft(alpha)*dx;
Gamma = abs(alpha_fft).^2 / L ;

ind_max = floor(length(kx)/2);
Gamma = Gamma(1:ind_max);

cumsum_Gamma = cumsum(Gamma(2:end));
p = mean(cumsum_Gamma(end-30:end));

th = p * beta;
e = find(cumsum_Gamma>th, 1, 'first');
qmax = kx(e)/2/pi; 

x = (kx(2:e)/2/pi).^2;
y = Gamma(2:e);    

c = qmax;
fun = @(a, x) a*x.^2 - a*c^2;

a0 = 1;
opts = optimoptions('lsqcurvefit', 'Display','off');
a_fit = lsqcurvefit(fun, a0, x, y', [],[], opts);

xval_g = (kx(2:e)/2/pi);
yval_g = a_fit * xval_g.^2 - a_fit*c^2;
Gamma0 = -a_fit*c^2;

