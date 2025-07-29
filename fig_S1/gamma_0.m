
clear

rng(15)

d_in = '../data/';
d_save = './fig/';

roi = 'cc';
brain_h = 'contra';

load([d_in 'cc_morph.mat']);
h_25_contra = morph.h_25_contra;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand_idx = 100;
rm = 15;

area_i = h_25_contra(rand_idx).area;
area_i = area_i(rm:end-rm);
area_i = util_remove_zeros(area_i); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(area_i);
step_sz = 2;
voxel_size = [0.05, 0.05, 0.05];
dx = voxel_size(3) * step_sz;
L = N*dx;
dk = 2*pi/L;
kmax = pi/dx;
kx = 0:dk:2*kmax-dk;
beta = 1 - 0.07;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Gamma, Gamma0, ind_max, xval, yval] = util_gamma0(area_i, L, dx, beta);
cumsum_Gamma = cumsum(Gamma(2:end));
p = mean(cumsum_Gamma(end-30:end));
th = p * beta;
e = find(cumsum_Gamma>th, 1, 'first');
qmax = kx(e)/2/pi; 

color_list = lines(5);

% Log scale
figure('unit','inch','position',[0 0 3 3])
plot(kx(1:ind_max)/2/pi, Gamma, 'color', color_list(1,:), 'linewidth',1); hold on
plot(kx(2:ind_max)/2/pi, cumsum_Gamma, 'color', color_list(2,:), 'linewidth',1.2); hold on    
plot(xval, yval, '--', 'LineWidth', 1.5, 'color', color_list(3,:));    
xline(qmax, '--k', 'LineWidth', 1.2)
set(gca, 'xscale', 'log', 'yscale', 'log');
grid on; grid minor
xlim([10^-3 5]); ylim([1e-6 50]); pbaspect([1 1 1]);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(gca,'xtick',[1e-2 1e-1 1e0 1e1],'xticklabels',{'$10^{-2}$' '$10^{-1}$' '$10^0$' '$10^1$'})
set(gca,'ytick',[1e-4 1e-2 1e0],'xticklabels',{'$10^{-4}$' '$10^{-2}$' '$10^0$'})
set(gca, 'fontsize',12)
xtickangle(0)
% exportgraphics(gca,[d_save 'gamma0.png'],'Resolution',300)

% Linear scale
figure('unit','inch','position',[0 0 2 2])
plot(kx(1:ind_max)/2/pi, Gamma, 'color', color_list(1,:), 'linewidth',1); hold on
plot(kx(2:ind_max)/2/pi, cumsum_Gamma, 'color', color_list(2,:), 'linewidth',1.2); hold on    
plot(xval, yval, '--', 'LineWidth', 1.5, 'color', color_list(3,:));    
xline(qmax, '--k', 'LineWidth', 1.2)
set(gca, 'xscale', 'log', 'yscale', 'linear');
grid on; grid minor
xlim([0 5]); ylim([1e-7 10]); pbaspect([1 1 1]);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(gca,'xtick',[1e-2 1e-1 1e0 1e1],'xticklabels',{'$10^{-2}$' '$10^{-1}$' '$10^0$' '$10^1$'})
xtickangle(0)
set(gca, 'fontsize',12)  
% exportgraphics(gca,[d_save 'gamma0_lin.png'],'Resolution',300)


%%

clear

d_in = '../data/';
d_save = './fig/';

load([d_in, 'simulation_synth.mat']);
area_all = sim.area;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand_idx = 10;
area_i = area_all(rand_idx,:);
area_i = util_remove_zeros(area_i); 
area_i = area_i(1:2:end); area_i = area_i(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(area_i);
step_sz = 2;
voxel_size = [0.05, 0.05, 0.05];
dx = voxel_size(3) * step_sz;
L = N*dx;
dk = 2*pi/L;
kmax = pi/dx;
kx = 0:dk:2*kmax-dk;
beta = 1 - 0.02;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Gamma, Gamma0, ind_max, xval, yval] = util_gamma0(area_i, L, dx, beta);
cumsum_Gamma = cumsum(Gamma(2:end));
p = mean(cumsum_Gamma(end-30:end));
th = p * beta;
e = find(cumsum_Gamma>th, 1, 'first');
qmax = kx(e)/2/pi; 

color_list = lines(5);

% Log scale
figure('unit','inch','position',[0 0 3 3])
plot(kx(1:ind_max)/2/pi, Gamma, 'color', color_list(4,:), 'linewidth',1); hold on
plot(kx(2:ind_max)/2/pi, cumsum_Gamma, 'color', color_list(2,:), 'linewidth',1.2); hold on    
plot(xval, yval, '--', 'LineWidth', 1.5, 'color', color_list(3,:));    
xline(qmax, '--k', 'LineWidth', 1.2)
set(gca, 'xscale', 'log', 'yscale', 'log');
grid on; grid minor
xlim([10^-3 5]); ylim([1e-6 50]); pbaspect([1 1 1]);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(gca,'xtick',[1e-2 1e-1 1e0 1e1],'xticklabels',{'$10^{-2}$' '$10^{-1}$' '$10^0$' '$10^1$'})
set(gca,'ytick',[1e-4 1e-2 1e0],'xticklabels',{'$10^{-4}$' '$10^{-2}$' '$10^0$'})
set(gca, 'fontsize',12)
xtickangle(0)
% exportgraphics(gca,[d_save 'synth_gamma0.png'],'Resolution',300)


% Linear scale
figure('unit','inch','position',[0 0 2 2])
plot(kx(1:ind_max)/2/pi, Gamma, 'color', color_list(4,:), 'linewidth',1); hold on
plot(kx(2:ind_max)/2/pi, cumsum_Gamma, 'color', color_list(2,:), 'linewidth',1.2); hold on    
plot(xval, yval, '--', 'LineWidth', 1.5, 'color', color_list(3,:));    
xline(qmax, '--k', 'LineWidth', 1.2)
set(gca, 'xscale', 'log', 'yscale', 'linear');
grid on; grid minor
xlim([0 5]); ylim([1e-7 26]); pbaspect([1 1 1]);

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(gca,'xtick',[1e-2 1e-1 1e0 1e1],'xticklabels',{'$10^{-2}$' '$10^{-1}$' '$10^0$' '$10^1$'})
xtickangle(0)
set(gca, 'fontsize',12)  
% exportgraphics(gca,[d_save 'synth_lin.png'],'Resolution',300)
