

% Validating the predicted amplitude c_D from geometry against MC simulatioins
% for individual synthetic and EM axons
 
clear

d = '../data/';
load([d, 'simulation_synth.mat']);
d_save = './fig/';

voxel_size = [0.036, 0.036, 0.036];
N = length(sim.area(1,:));
dx = voxel_size(3);
L = N*dx;

% color bar
rgb_points = [
    255, 226, 255;
    255, 21, 255;
    199, 0, 199;
    110, 10, 150;
    52, 0, 54] / 255;

orig_pos = linspace(1, 100, size(rgb_points, 1));
n = 100; new_pos = 1:n;
cmap = [
    interp1(orig_pos, rgb_points(:,1), new_pos);
    interp1(orig_pos, rgb_points(:,2), new_pos);
    interp1(orig_pos, rgb_points(:,3), new_pos)]';


% Error bars reflect errors in estimating the plateau Gamma_0 from \Gamma_\eta(q)
% for a fixed t0 = 10.
t0 = 10;

num_iter = 10;

th = 0.004; perc = th*0.5;
range = linspace(th-perc, th+perc, num_iter);
beta = 1 - range;

cDs_g = zeros(1,3,num_iter);
for ii = 1:size(beta,2)
    beta_i = beta(ii);
    for i = 1:size(sim.area,1)
        area_i = sim.area(i,:); area_i = area_i(:);
        AD_i = sim.axial_D(i,:);
        cDs_g(i,:,ii) = util_cD(area_i, AD_i, L, dx, beta_i, sim.t, t0);
    end
end


% Error bars reflect errors in estimating c_D from MC-simulated D(t) 
% for a fixed beta = 1 - 0.004;
beta = 1 - 0.004;

num_iter = 10;
t0 = 5; t0_end = t0+10;
range = linspace(t0, t0_end, num_iter);
cDs_t = zeros(1,3,num_iter);
for ii = 1:size(range,2)
    t0 = range(ii);
    for i = 1:size(sim.area,1)
        area_i = sim.area(i,:); area_i = area_i(:);
        area_i = util_remove_zeros(area_i);
        AD_i = sim.axial_D(i,:);
        cDs_t(i,:,ii) = util_cD(area_i, AD_i, L, dx, beta, sim.t, t0);
    end
end

% Color bars are defined based on coeff. of variations CV.
lin_cvs = linspace(0,0.5,n);
[~, cl_inx] = min(pdist2(lin_cvs', cDs_g(:,3).^2), [], 1);


figure('unit','inch','position',[0 0 2 2]); hold on
plot([0 3.5], [0 3.5], '--k', 'LineWidth', 1.2);

x = cDs_t(:,1,:); x_av = mean(x,3); x_std = std(x,0,3);
y = cDs_g(:,2,:); y_av = mean(y,3); y_std = std(y,0,3);

errorbar(x_av,y_av,x_std,'horizontal', 'LineStyle','none', 'Color','k'); hold on
errorbar(x_av,y_av,y_std,'vertical', 'LineStyle','none', 'Color','k'); hold on

scatter(x_av, y_av, 15, cmap(cl_inx,:), 'filled', ...
    'MarkerFaceAlpha', 1, ...
    'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerEdgeAlpha', 0.5); 

set(gca,'xtick',[0:0.5:1.5],'ytick',[0:0.5:1.5])
box on; grid on; set(gca,'GridLineStyle','--')
xlim([0 1.5]); ylim([0 1.5]); pbaspect([1 1 1]);

% exportgraphics(gca, [d_save 'cd_synth.png'],'Resolution',300)

% Generate just a colorbar for synth axons
colormap(cmap)
c = colorbar;
clim([0.2 0.7])
c.Ticks = [0.2:0.25:0.7];


%%

clear

d = '../data/';
load([d, 'simulation_EM.mat']);
d_save = './fig/';

voxel_size = [0.05, 0.05, 0.05];
dx = voxel_size(3);

% color bar
rgb_points = [
    200, 255, 200;
    75, 209, 75;
    33, 187, 33;
    20, 122, 20;
    10, 41, 10;
    ] / 255;

orig_pos = linspace(1, 100, size(rgb_points, 1));
n = 100; new_pos = 1:n;
cmap_green = [
    interp1(orig_pos, rgb_points(:,1), new_pos);
    interp1(orig_pos, rgb_points(:,2), new_pos);
    interp1(orig_pos, rgb_points(:,3), new_pos)]';

rgb_points = [
    255, 200, 200;
    209, 75, 75;
    187, 33, 33;
    122, 20, 20;
    41, 10, 10;
    ] / 255;

orig_pos = linspace(1, 100, size(rgb_points, 1));
n = 100; new_pos = 1:n;
cmap_red = [
    interp1(orig_pos, rgb_points(:,1), new_pos);
    interp1(orig_pos, rgb_points(:,2), new_pos);
    interp1(orig_pos, rgb_points(:,3), new_pos)]';


% Error bars reflect errors in estimating the plateau Gamma_0 from \Gamma_\eta(q)
% for a fixed t0 = 10.
t0 = 10;

num_iter = 10;
th = 0.04; perc = th*0.5;
range = linspace(th-perc,th+perc,num_iter);
beta = 1 - range;

cDs_g = zeros(1,3,num_iter);
groups = zeros(1);
for ii = 1:size(beta,2)
    beta_i = beta(ii);
    for i = 1:size(sim.area, 1)
        fn = sim.ids(i).name;
        groups(i) = 1;
        if contains(fn,'LM_25') | contains(fn,'LM_49')
            groups(i) = 0;
        end           
        area_i = sim.area(i,:); area_i = area_i(:);
        area_i = util_remove_zeros(area_i);

        N = length(area_i);
        L = N*dx;
        dk = 2*pi/L;
        kmax = pi/dx;
        kx = 0:dk:2*kmax-dk;
        
        AD_i = sim.axial_D(i,:);
        cDs_g(i,:,ii) = util_cD(area_i, AD_i, L, dx, beta_i, sim.t, t0);             
    end
end


% Error bars reflect errors in estimating c_D from MC-simulated D(t) 
% for a fixed beta = 1 - 0.04;
beta = 1 - 0.04;

num_iter = 10;
t0 = 3; t_end = t0+10;
range = linspace(t0,t_end,num_iter);

cDs_t = zeros(1,3,num_iter);
for ii = 1:size(range,2)
    t0 = range(ii);
    for i = 1:size(sim.area, 1)
        area_i = sim.area(i,:); area_i = area_i(:);
        area_i = util_remove_zeros(area_i);

        N = length(area_i);
        L = N*dx;
        dk = 2*pi/L;
        kmax = pi/dx;
        kx = 0:dk:2*kmax-dk;
        AD_i = sim.axial_D(i,:);
        
        cDs_t(i,:,ii) = util_cD(area_i, AD_i, L, dx, beta, sim.t, t0);
    end
end


% Color bars are defined based on coeff. of variations CV.
lin_cvs = linspace(0,0.5,n);
[~, cl_inx] = min(pdist2(lin_cvs',cDs_g(:,3,1).^2),[],1);

figure('unit','inch','position',[0 0 2 2]); hold on
plot([0 3.5], [0 3.5], '--k', 'LineWidth', 1.2);
x = cDs_t(:,1,:); x_av = mean(x,3); x_std = std(x,0,3);
y = cDs_g(:,2,:); y_av = mean(y,3); y_std = std(y,0,3);

x_best = x(:,:,2);
y_best = y(:,:,6);

errorbar(x_av,y_av,x_std,'horizontal', 'LineStyle','none', 'Color','k'); hold on
errorbar(x_av,y_av,y_std,'vertical', 'LineStyle','none', 'Color','k'); hold on

inx = groups == 0; 
scatter(x_av(inx), y_av(inx), 15, cmap_green(cl_inx(inx),:), 'filled', ...
    'MarkerFaceAlpha', 1, ...
    'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerEdgeAlpha', 0.5); 

inx = groups == 1; 
scatter(x_av(inx), y_av(inx), 15, cmap_red(cl_inx(inx),:), 'filled', ...
    'MarkerFaceAlpha', 1, ...
    'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerEdgeAlpha', 0.5); 

set(gca,'xtick',[0:0.5:2],'ytick',[0:0.5:2])
box on; grid on; set(gca,'GridLineStyle','--')
xlim([0 2]); ylim([0 2]); pbaspect([1 1 1]);

% exportgraphics(gca, [d_save 'cd_real.png'],'Resolution',300)


% Generate just a colorbar for TBIs
figure('unit','inch','position',[0 0 2 2]); hold on
surf(peaks); colormap(cmap_red)
c = colorbar;
clim([0.2 0.7])
c.Ticks = [0.2:0.25:0.7];
% exportgraphics(gca, [d_save 'colorbar_red.png'],'Resolution',300)


% Generate just a colorbar for Shams
figure('unit','inch','position',[0 0 2 2]); hold on
surf(peaks); colormap(cmap_green)
c = colorbar;
clim([0.2 0.7])
c.Ticks = [0.2:0.25:0.7];
% exportgraphics(gca, [d_save 'colorbar_green.png'],'Resolution',300)