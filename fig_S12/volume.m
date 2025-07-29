
clear; clc

% [tortuosity, Gamma_0, sinosity, volume, l_gd, l_eq]; 

rng(15)

d_in = '../data/';

roi = 'cg';
brain_h = 'ipsi';

load([d_in roi '_geo_param.mat']);

x1 = getfield(gamma_tortuosity, ['h_25_' brain_h]);
x2 = getfield(gamma_tortuosity, ['h_49_' brain_h]);
x3 = getfield(gamma_tortuosity, ['t_2_'  brain_h]);
x4 = getfield(gamma_tortuosity, ['t_24_' brain_h]);
x5 = getfield(gamma_tortuosity, ['t_28_' brain_h]);

x1 = [x1, 1*ones(length(x1),1)];
x2 = [x2, 2*ones(length(x2),1)];
x3 = [x3, 3*ones(length(x3),1)];
x4 = [x4, 4*ones(length(x4),1)];
x5 = [x5, 5*ones(length(x5),1)];

sham_x = [x1; x2];
[r,~] = find(sham_x<=0); sham_x(r,:) = [];
tbi_x = [x3; x4; x5];
[r,~] = find(tbi_x<=0); tbi_x(r,:) = [];

sham_x = sham_x(randperm(length(sham_x)), :);
tbi_x = tbi_x(randperm(length(tbi_x)), :);
sz_ = min(length(sham_x), length(tbi_x));

sham_smpl_x = sham_x(1:sz_, :);
tbi_smpl_x = tbi_x(1:sz_, :);

sham_v = sham_smpl_x(:,4);
tbi_v = tbi_smpl_x(:,4);

sham_l = sham_smpl_x(:,5);
tbi_l = tbi_smpl_x(:,5);


%% v distribution

d_fig = '/Users/aliabd/Projects/Project_unmAxons/NYU_DL/proj_unmyel/fig_volume/';

sham_cl = [40 163 127]/255;
tbi_cl = [210 82 24]/255;

figure('unit','inch','position',[0 0 2 2])
h1 = histogram(sham_v, 'FaceColor', sham_cl); hold on
h2 = histogram(tbi_v, 'FaceColor', tbi_cl); hold on
h1.Normalization = 'probability';
h1.BinWidth = 4; h1.EdgeAlpha = 0.1; h1.FaceAlpha = 1;
h2.Normalization = 'probability';
h2.BinWidth = 4; h2.EdgeAlpha = 0.1; h2.FaceAlpha = 0.75;
xlim([0 150]); pbaspect([1 1 1]);

set(gca,'xtick',0:50:150,'ytick',[], 'linewidth',0.7)
xtickangle(0); grid off

xline(median(sham_v), '--', 'Color', sham_cl, 'linewidth',1)
xline(median(tbi_v), '--', 'Color', tbi_cl, 'linewidth',1)

effect = util_nonparam_effect_size(sham_v, tbi_v);
effect = round(-effect*100)/100;
title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)

% exportgraphics(gca,[d_fig roi '_' brain_h '_volume.png'], 'Resolution', 300)


%%

figure('unit','inch','position',[0 0 2 2])
h1 = histogram(sham_l, 'FaceColor', sham_cl); hold on
h2 = histogram(tbi_l, 'FaceColor', tbi_cl); hold on
h1.Normalization = 'probability';
h1.BinWidth = 4; h1.EdgeAlpha = 0.1; h1.FaceAlpha = 1;
h2.Normalization = 'probability';
h2.BinWidth = 4; h2.EdgeAlpha = 0.1; h2.FaceAlpha = 0.75;
xlim([0 150]); pbaspect([1 1 1]);

set(gca,'xtick',0:50:200,'ytick',[], 'linewidth',0.7)
xtickangle(0); grid off

xline(median(sham_l), '--', 'Color', sham_cl, 'linewidth',1)
xline(median(tbi_l), '--', 'Color', tbi_cl, 'linewidth',1)

effect = util_nonparam_effect_size(sham_l, tbi_l);
effect = round(-effect*100)/100;
title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)

% exportgraphics(gca,[d_fig roi '_' brain_h '_length.png'], 'Resolution', 300)
