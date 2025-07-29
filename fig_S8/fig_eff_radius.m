
clear; clc

rng(15)

d_in = '../data/';
d_save = './fig/';

sham_cl = [40 163 127]/255;
tbi_cl = [210 82 24]/255;

% Define 'cg' or 'cc'
roi = 'cg';

% Define 'ipsi' or 'contra'
brain_h = 'ipsi';

load([d_in roi '_morph.mat']);

x1 = getfield(morph, ['h_25_' brain_h]);
x2 = getfield(morph, ['h_49_' brain_h]);
x3 = getfield(morph, ['t_2_'  brain_h]);
x4 = getfield(morph, ['t_24_' brain_h]);
x5 = getfield(morph, ['t_28_' brain_h]);

sham_x = [x1, x2];
tbi_x = [x3, x4, x5];

inx = randperm(length(sham_x));
sham_x = sham_x(inx);

inx = randperm(length(tbi_x));
tbi_x = tbi_x(inx);

sz = min(length(sham_x), length(tbi_x));
sham_smpl_x = sham_x(1:sz);
tbi_smpl_x = tbi_x(1:sz);

rm = 15;
sham_v = zeros(sz,1);
sham_r = [];
sham_r_eff = [];
sham_walpha = [];
for ii = 1:sz
    area_i = sham_smpl_x(ii).area;
    area_i = area_i(rm:end-rm);
    area_i = util_remove_zeros(area_i); 
    
    sham_v(ii,1) = sum(area_i);

    r = sqrt(area_i/pi);
    sham_r = [sham_r; r];

    eff_r = mean(r.^6) ./ mean(r.^2);
    sham_r_eff = [sham_r_eff; eff_r];
  
    alpha = area_i ./ mean(area_i);
    walpha = 1. / mean(1 ./ alpha);
    sham_walpha = [sham_walpha, walpha];
end

tbi_v = zeros(sz,1);
tbi_r = [];
tbi_r_eff = [];
tbi_walpha = [];
for ii = 1:sz
    area_i = tbi_smpl_x(ii).area;
    area_i = area_i(rm:end-rm);
    area_i = util_remove_zeros(area_i); 
    
    tbi_v(ii,1) = sum(area_i);

    r = sqrt(area_i/pi);
    tbi_r = [tbi_r; r];

    eff_r = mean(r.^6) ./ mean(r.^2);
    tbi_r_eff = [tbi_r_eff; eff_r];
  
    alpha = area_i ./ mean(area_i);
    walpha = 1. / mean(1 ./ alpha);
    tbi_walpha = [tbi_walpha, walpha];
end

sham_w = sham_v./sum(sham_v);
tbi_w = tbi_v./sum(tbi_v);

%% r distribution

figure('unit','inch','position',[0 0 2 2])
h1 = histogram(sham_r, 'FaceColor', sham_cl); hold on
h2 = histogram(tbi_r, 'FaceColor', tbi_cl); hold on
h1.Normalization = 'probability';
h1.BinWidth = 0.02; h1.EdgeAlpha = 0.1; h1.FaceAlpha = 1;
h2.Normalization = 'probability';
h2.BinWidth = 0.02; h2.EdgeAlpha = 0.1; h2.FaceAlpha = 0.75;
xlim([0 1.5]); pbaspect([2 1 1]);

set(gca,'xtick',0:0.5:1.5,'ytick',[], 'linewidth',0.7)
xtickangle(0); grid off

xline(median(sham_r), '--', 'Color', sham_cl, 'linewidth',1)
xline(median(tbi_r), '--', 'Color', tbi_cl, 'linewidth',1)

effect = util_nonparam_effect_size(sham_r, tbi_r);
effect = round(-effect*100)/100;
title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)

% exportgraphics(gca,[d_save roi '_' brain_h '_radius.png'], 'Resolution', 300)


%% w * r_eff^4

ws_rp = ceil(1e3 * (sham_w / max(sham_w)));
sham_effR_rp = repelem(sham_r_eff, ws_rp); 

wt_rp = ceil(1e3 * (tbi_w / max(tbi_w)));
tbi_effR_rp = repelem(tbi_r_eff, wt_rp);

figure('unit','inch','position',[0 0 2 2])
h1 = histogram(sham_effR_rp, 'FaceColor', sham_cl); hold on
h2 = histogram(tbi_effR_rp, 'FaceColor', tbi_cl); hold on
h1.Normalization = 'probability';
h1.BinWidth = 0.001; h1.EdgeAlpha = 0.1; h1.FaceAlpha = 1;
h2.Normalization = 'probability';
h2.BinWidth = 0.001; h2.EdgeAlpha = 0.1; h2.FaceAlpha = 0.75;
xlim([0 0.15]); pbaspect([2 1 1]);

set(gca,'xtick',0:0.05:0.15,'ytick',[], 'linewidth',0.7)
box on; h = gca; h.YAxis.Visible = 'on';
xtickangle(0); grid off

xline(median(sham_effR_rp), '--', 'Color', sham_cl, 'linewidth',1)
xline(median(tbi_effR_rp), '--', 'Color', tbi_cl, 'linewidth',1)

effect = util_nonparam_effect_size(sham_effR_rp, tbi_effR_rp);
effect = round(-effect*100)/100;
title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)

% exportgraphics(gca,[d_save roi '_' brain_h '_eff_radius.png'], 'Resolution', 300)


%% w / < 1/alpha >

sham_walpha_rp = repelem(sham_walpha, ws_rp); 
tbi_walpha_rp = repelem(tbi_walpha, wt_rp);

figure('unit','inch','position',[0 0 2 2])
h1 = histogram(sham_walpha_rp, 'FaceColor', sham_cl); hold on
h2 = histogram(tbi_walpha_rp, 'FaceColor', tbi_cl);
h1.Normalization = 'probability';
h1.BinWidth = 0.02; h1.EdgeAlpha = 0.1; h1.FaceAlpha = 1;
h2.Normalization = 'probability';
h2.BinWidth = 0.02; h2.EdgeAlpha = 0.1; h2.FaceAlpha = 0.75;
xlim([0 1]); pbaspect([2 1 1]);

box on; h = gca; h.YAxis.Visible = 'on';
set(gca,'xtick',0:0.5:1.5,'ytick',[], 'linewidth',.7)
xtickangle(0); grid off

xline(median(sham_walpha_rp), '--', 'Color', sham_cl, 'linewidth',1)
xline(median(tbi_walpha_rp), '--', 'Color', tbi_cl, 'linewidth',1)

effect = util_nonparam_effect_size(sham_walpha_rp, tbi_walpha_rp);
effect = round(-effect*100)/100;
title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)

% exportgraphics(gca,[d_save roi '_' brain_h '_cv2.png'], 'Resolution', 300)

