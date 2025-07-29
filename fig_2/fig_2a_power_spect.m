
% Power spectrum

% Synthetic axon-----------------------------------------------------------

clear; close all; clc

d_in = './data2/';
d_save = './fig/';

fn = 'axon_0.33_47.mat';
fn_in = [d_in fn];
bin_obj = load(fn_in); bin_obj = bin_obj.axon;

color_ = [148 8 125]/255;

% Gamma0
voxel_size = [0.036, 0.036, 0.036];
N = size(bin_obj,3);
dx = voxel_size(3);
L = N*dx;
dk = 2*pi/L;
kmax = pi/dx;
kx = 0:dk:2*kmax-dk;
beta = 1 - 0.004; % From Eq. S1

area_i = util_fast_area_quantification(bin_obj, voxel_size(3));
[Gamma, ~, ind_max, xval, yval] = util_gamma0(area_i, L, dx, beta);

figure('unit','inch','position',[0 0 3 3])
plot(kx(1:ind_max)/2/pi, Gamma, 'color', color_, 'linewidth',1); hold on
plot(xval, yval, '--', 'LineWidth', 1.2, 'Color', color_/1.6);    
set(gca, 'xscale', 'log', 'yscale', 'log'); hold on; grid on

% EM axon------------------------------------------------------------------

fn = 'cc_LM_25_contra_axon_1';
fn_in = [d_in fn];
bin_obj = load(fn_in); bin_obj = bin_obj.axon;

color_ = [48 172 66]/255;

% skeleton parameters
voxel_size = [0.05, 0.05, 0.05];
skl_smooth = 20;
step_sz = 1;
fineness = 1000;
max_r = 7;
grid_res = 0.25;

skel_i = util_skeleton(bin_obj, skl_smooth);    
eq_sampled_skel = util_equidistance_sampling(skel_i, fineness, step_sz);
stat = util_xsection_plane(bin_obj, eq_sampled_skel, max_r, grid_res, voxel_size);
area_i = [stat.Area];
area_i = util_remove_zeros(area_i); area_i = area_i(:);

% Gamma0
N = length(area_i);
dx = voxel_size(3);
L = N * dx;
dk = 2*pi/L;
kmax = pi/dx;
kx = 0:dk:2*kmax-dk;
beta = 1 - 0.07;

[Gamma, ~, ind_max, xval, yval] = util_gamma0(area_i, L, dx, beta);

plot(kx(1:ind_max)/2/pi, Gamma, 'color', color_, 'linewidth',1); hold on
plot(xval, yval, '--', 'LineWidth', 1.2, 'Color', color_/1.6);    
xlim([0 1.5e1]); ylim([0 1.7e1]); pbaspect([1 1 1]);
xticks([1e-2, 1]); yticks(10.^[-3:2:1]); grid on;
set(gca,'fontsize',12)

% EM axon TBI ------------------------------------------------------------------

fn = 'cc_LM_28_ipsi_axon_962';
fn_in = [d_in fn];
bin_obj = load(fn_in); bin_obj = bin_obj.axon;

color_ = [210 82 24]/255;

skel_i = util_skeleton(bin_obj, skl_smooth);    
eq_sampled_skel = util_equidistance_sampling(skel_i, fineness, step_sz);
stat = util_xsection_plane(bin_obj, eq_sampled_skel, max_r, grid_res, voxel_size);
area_i = [stat.Area]';
area_i = util_remove_zeros(area_i);

N = length(area_i);
dx = voxel_size(3);
L = N * dx;
dk = 2*pi/L;
kmax = pi/dx;
kx = 0:dk:2*kmax-dk;

% Gamma0
[Gamma, ~, ind_max, xval, yval] = util_gamma0(area_i, L, dx, beta);

plot(kx(1:ind_max)/2/pi, Gamma, 'color', color_, 'linewidth',1); hold on
plot(xval, yval, '--', 'LineWidth', 1.2, 'Color', color_/1);    
xlim([0 1.5e1]); ylim([0 10]); pbaspect([1 1 1]);
xticks([1e-2, 1]); yticks(10.^[-3:2:1]); grid on;
set(gca,'fontsize',12)

% exportgraphics(gca, [d_save 'power_spectrum.png'],'Resolution',300)