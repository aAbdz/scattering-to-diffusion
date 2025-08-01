
% Visualization of full axon length

clear; clc; close all

d_in = './data1/';
d_save = './fig/';

tbi_color = [217 108 58]/255;

% Axon's name in the ../data directory

% CC axons
% cc_LM_24_contra_axon_1000
% cc_LM_24_contra_axon_1014
% cc_LM_24_contra_axon_1029
% cc_LM_24_contra_axon_1042
 
% Cg axons
% cg_LM_24_contra_axon_10
% cg_LM_24_contra_axon_100
% cg_LM_24_contra_axon_101
% cg_LM_24_contra_axon_1001

fn = 'cg_LM_24_contra_axon_10.mat';

fn_in = [d_in fn];
bin_obj = load(fn_in); bin_obj = bin_obj.axon;

figure
p = patch(isosurface(smooth3(bin_obj,'gaussian',5),0.5));
p.FaceColor = tbi_color; p.EdgeColor = 'none'; 
camlight; material dull; lighting gouraud;
view(30, 50)
set(gca,'Visible','off', 'Color','k'); axis equal; grid off; axis off
box off; grid off;

% exportgraphics(gca, [d_save fn '.png'],'BackgroundColor','none','Resolution',300)


%% Visualization of a small axon segment

% The following two axons
% cg_LM_24_contra_axon_10 and 
% cc_LM_24_contra_axon_1014, where used for visualizing small axon segments. 

clear;  clc; close all

d_in = './data1/';
d_save = './fig/';

tbi_color = [217 108 58]/255;
voxel_size = [0.05, 0.05, 0.05]; %in microns

% Axon's name in the ./data directory
fn = 'cg_LM_24_contra_axon_10';

fn_in = [d_in fn];
bin_obj = load(fn_in); bin_obj = bin_obj.axon;

% Axon segment; 
% Axon is clipped along the z-axis between 700-900 
bin_obj_sec = bin_obj(:,:,700:900);


% Calculate area along z-axis
area_i = util_fast_area_quantification(bin_obj_sec, voxel_size(3));

figure
p = patch(isosurface(smooth3(bin_obj_sec,'gaussian',5),0.5));
p.FaceColor = tbi_color; p.EdgeColor = 'none'; 
camlight; material dull; lighting gouraud;
view(-105, 17)
set(gca,'Visible','off', 'Color','k'); axis equal; grid off; axis off
box off; grid off;
% exportgraphics(gca, [d_save fn '_seg.png'],'BackgroundColor','none','Resolution',300)


figure('unit','inch','position',[0 0 1.2 1.2])
h1 = area(area_i, 'FaceColor', tbi_color); hold on
h1.EdgeAlpha = 0; h1.FaceAlpha = 1;
ylim([0 1.2]); pbaspect([2 1 1]);
set(gca,'xtick',[0,100,200], 'xticklabel',[],'ytick',[0,1],'yticklabel',[],'TickDir', 'out')
xtickangle(0); grid off; box off
% exportgraphics(gca, [d_save fn '_seg_area.png'],'Resolution',300)