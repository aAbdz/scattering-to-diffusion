
clear

rng(15)

d_in = '../data/';
d_save = './fig/';

sham_cl = [40 163 127]/255;
tbi_cl = [210 82 24]/255;


% Define 'cg' or 'cc'
roi = 'cg';

% Define 'ipsi' or 'contra'
brain_h = 'ipsi';

% column1 = tortuosity % column2 = Gamma0
% column3 = sinuosity  % column4 = axon volume
% column5 = geodesic   % column6 = Euclidean length
load([d_in roi '_geo_param.mat']);

x1 = getfield(gamma_tortuosity, ['h_25_' brain_h]);
x2 = getfield(gamma_tortuosity, ['h_49_' brain_h]);
x3 = getfield(gamma_tortuosity, ['t_2_'  brain_h]);
x4 = getfield(gamma_tortuosity, ['t_24_' brain_h]);
x5 = getfield(gamma_tortuosity, ['t_28_' brain_h]);

sham_x = [x1; x2];
[r,~] = find(sham_x<=0); sham_x(r,:) = [];

tbi_x = [x3; x4; x5];
[r,~] = find(tbi_x<=0); tbi_x(r,:) = [];

inx = randperm(length(sham_x));
sham_x = sham_x(inx, :);

inx = randperm(length(tbi_x));
tbi_x = tbi_x(inx, :);

sz_ = min(length(sham_x), length(tbi_x));
sham_smpl_x = sham_x(1:sz_, :);
tbi_smpl_x = tbi_x(1:sz_, :);

% volume of every axon
sham_v = sum(sham_smpl_x(:,4));
tbi_v = sum(tbi_smpl_x(:,4));

% sinuosity 
sham_s = sham_smpl_x(:,3);
tbi_s = tbi_smpl_x(:,3);


%%
% sinuosity distribution

figure('unit','inch','position',[0 0 2 2])
h1 = histogram(sham_s, 'FaceColor', sham_cl); hold on
h2 = histogram(tbi_s, 'FaceColor', tbi_cl); hold on
h1.Normalization = 'probability';
h1.BinWidth = 0.002; h1.EdgeAlpha = 0.1; h1.FaceAlpha = 1;
h2.Normalization = 'probability';
h2.BinWidth = 0.002; h2.EdgeAlpha = 0.1; h2.FaceAlpha = 0.75;
xlim([1 1.1]); pbaspect([2 1 1]);

set(gca,'xtick',0:0.05:1.1,'ytick',[], 'linewidth',0.7)
xtickangle(0); grid off

xline(median(sham_s), '--', 'Color', sham_cl, 'linewidth',1)
xline(median(tbi_s), '--', 'Color', tbi_cl, 'linewidth',1)

effect = util_nonparam_effect_size(sham_s, tbi_s);
effect = round(-effect*100)/100;
title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)

% exportgraphics(gca,[d_save roi '_' brain_h '_sinuosity.png'], 'Resolution', 300)


