
% Axon from sham dataset

clear; clc; close all

d_in = './data/';
d_save = './fig/';

color_ = [48 172 66]/255;
voxel_size = [0.05, 0.05, 0.05];

skl_smooth = 20;
step_sz = 1;
fineness = 1000;
max_r = 7;
grid_res = 0.25;

% Axon's name in the ./data directory
fn = 'cc_LM_25_contra_axon_1';
fn_in = [d_in fn];

bin_obj = load(fn_in); bin_obj = bin_obj.axon;

% Skeleton and xsectional analysis
skel_i = util_skeleton(bin_obj, skl_smooth);    
eq_sampled_skel = util_equidistance_sampling(skel_i, fineness, step_sz);
stat = util_xsection_plane(bin_obj, eq_sampled_skel, max_r, grid_res, voxel_size);

area_i = [stat.Area]';
area_i = util_remove_zeros(area_i);
mean_area_i = mean(area_i);
alpha = (area_i - mean_area_i) / mean_area_i; 

figure
p = patch(isosurface(smooth3(bin_obj,'gaussian',5),0.5));
p.FaceColor = color_; p.EdgeColor = 'none'; 
view(40, 0)
camlight; material dull; lighting gouraud;
box off; axis equal; grid off; axis off
exportgraphics(gca, [d_save fn '.png'],'BackgroundColor','none','Resolution',300)

figure('unit','inch','position',[0 0 3 1.2])
h1 = area(alpha, 'FaceColor', color_); hold on
h1.EdgeAlpha = 0; h1.FaceAlpha = 1;
ylim([-1.2 3.1]); xlim([0 length(area_i)]); pbaspect([2 1 1]);
set(gca,'xtick',[],'ytick',[-1,0,1],'TickDir', 'out')
xtickangle(0); grid off; box off
exportgraphics(gca, [d_save fn '_hist.png'],'Resolution',300)



%% Axon from TBI dataset

clc; close all

color_ = [210 82 24]/255;

fn = 'cc_LM_28_ipsi_axon_962';
fn_in = [d_in fn];
bin_obj = load(fn_in); bin_obj = bin_obj.axon;

skel_i = util_skeleton(bin_obj, skl_smooth);    
eq_sampled_skel = util_equidistance_sampling(skel_i, fineness, step_sz);
stat = util_xsection_plane(bin_obj, eq_sampled_skel, max_r, grid_res, voxel_size);
area_i = [stat.Area]';
area_i = util_remove_zeros(area_i);
mean_area_i = mean(area_i);
alpha = (area_i - mean_area_i) / mean_area_i; 

figure
p = patch(isosurface(smooth3(bin_obj,'gaussian',5),0.5));
p.FaceColor = color_; p.EdgeColor = 'none'; 
camlight; material dull; lighting gouraud;
view(42, 42)
set(gca,'Visible','off', 'Color','k'); axis equal; grid off; axis off
box off; grid off;
exportgraphics(gca, [d_save fn '.png'],'BackgroundColor','none','Resolution',300)

figure('unit','inch','position',[0 0 3 1.2])
h1 = area(alpha, 'FaceColor', color_); hold on
h1.EdgeAlpha = 0; h1.FaceAlpha = 1;
ylim([-1.2 3.1]); xlim([0 length(area_i)]); pbaspect([2 1 1]);
set(gca,'xtick',[],'ytick',[-1,0,1],'TickDir', 'out')
xtickangle(0); grid off; box off
exportgraphics(gca, [d_save fn '_hist.png'],'Resolution',300)


%% Synthetic axon: 500 micron

clc; close all

color_ = [148 8 125]/255;

fn = 'axon_0.33_47.mat';
fn_in = [d_in fn];
bin_obj = load(fn_in); bin_obj = bin_obj.axon;

voxel_size = [0.036, 0.036, 0.036];
area_i = util_fast_area_quantification(bin_obj, voxel_size(3));
mean_area_i = mean(area_i);
alpha = (area_i - mean_area_i) / mean_area_i; 

% Show the first 100 micron of the axon
% 100/0.036
inx = floor(100 / voxel_size(3));

figure
bin_obj_sec = bin_obj(:,:,1:inx);
p = patch(isosurface(smooth3(bin_obj_sec,'gaussian',5),0.5));
p.FaceColor = color_; p.EdgeColor = 'none'; 
camlight; material dull; lighting gouraud;
view(-100, 0)
set(gca,'Visible','off', 'Color','k'); axis equal; grid off; axis off
box off; grid off;
exportgraphics(gca, [d_save fn '.png'],'BackgroundColor','none','Resolution',300)

alpha_sec = alpha(1:inx);
figure('unit','inch','position',[0 0 3 1.2])
h1 = area(alpha_sec, 'FaceColor', color_); hold on
h1.EdgeAlpha = 0; h1.FaceAlpha = 1;
ylim([-1.2 3.1]); xlim([0 length(alpha_sec)]); pbaspect([2 1 1]);
set(gca,'xtick',[],'ytick',[-1,0,1],'TickDir', 'out')
xtickangle(0); grid off; box off
exportgraphics(gca, [d_save fn '_hist_sec1.png'],'Resolution',300)



% Show the last 25 micron of the axon
inx = floor(25 / voxel_size(3));

figure
bin_obj_sec = bin_obj(:,:,end-inx:end);
p = patch(isosurface(smooth3(bin_obj_sec,'gaussian',5),0.5));
p.FaceColor = color_; p.EdgeColor = 'none'; 
camlight; material dull; lighting gouraud;
view(-100, 0)
set(gca,'Visible','off', 'Color','k'); axis equal; grid off; axis off
box off; grid off;
exportgraphics(gca, [d_save fn '_sec2.png'],'BackgroundColor','none','Resolution',300)

alpha_sec = alpha(end-inx:end);
figure('unit','inch','position',[0 0 3 1.2])
h1 = area(alpha_sec, 'FaceColor', color_); hold on
h1.EdgeAlpha = 0; h1.FaceAlpha = 1;
ylim([-1.2 3.1]); xlim([0 length(alpha_sec)]); pbaspect([2 1 1]);
set(gca,'xtick',[],'ytick',[-1,0,1],'yticklabel',[],'TickDir', 'out')
xtickangle(0); grid off; box off
exportgraphics(gca, [d_save fn '_hist_sec2.png'],'Resolution',300)









