
clear

rng(15)
d_in = './data/';
d_save = './fig/';

sham_color = [27 184 99]/255;
tbi_color = [213 91 36]/255;

sham_cl1 = [32, 219, 75]/255;
sham_cl2 = [40, 127, 81]/255;
tbi_cl1  = [210 82 24]/255;
tbi_cl2  = [251 3 27]/255;
tbi_cl3  = [255 140 0]/255;

all_color = [sham_cl1; sham_cl2; tbi_cl1; tbi_cl2; tbi_cl3];

bbox_cl = [50,50,120]/255;
proj_cl = [50,50,120]/255;
svm_cl  = [0,230,240]/255;


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

x1 = [x1, 1*ones(length(x1),1)];
x2 = [x2, 2*ones(length(x2),1)];
x3 = [x3, 3*ones(length(x3),1)];
x4 = [x4, 4*ones(length(x4),1)];
x5 = [x5, 5*ones(length(x5),1)];

% pool sham data
sham_x = [x1; x2];
[r,~] = find(sham_x<=0); sham_x(r,:) = [];

% pool tbi data
tbi_x = [x3; x4; x5];
[r,~] = find(tbi_x<=0); tbi_x(r,:) = [];

% sample to the size of minimum of sham and tbi sizes
sham_x = sham_x(randperm(length(sham_x)), :);
tbi_x = tbi_x(randperm(length(tbi_x)), :);
sz_ = min(length(sham_x), length(tbi_x));
sham_smpl_x = (sham_x(1:sz_,:));
tbi_smpl_x  = (tbi_x(1:sz_,:));

smpl_X = [sham_smpl_x; tbi_smpl_x];
smpl_Y = [ones(sz_, 1); 2*ones(sz_, 1)];


%% Histogram of tortuosity

figure('unit','inch','position',[0 0 2 2])
h1 = histogram(sham_smpl_x(:,1)); hold on
h2 = histogram(tbi_smpl_x(:,1));

h1.Normalization = 'probability'; h1.BinWidth = 0.05; 
h1.EdgeColor = sham_color; h1.FaceColor = sham_color; h1.FaceAlpha = 0.85;

h2.Normalization = 'probability'; h2.BinWidth = 0.05; 
h2.EdgeColor = tbi_color; h2.FaceColor = tbi_color; h2.FaceAlpha = 0.85;

xlim([1 3]); pbaspect([2 1 1]);
box off; h = gca; h.YAxis.Visible = 'off';
set(gca,'xtick',[],'ytick',[],'linewidth', .7)

xline(median(sham_smpl_x(:,1)), '--', 'LineWidth', 1, 'Color', max(sham_color-0.15,0))
xline(median(tbi_smpl_x(:,1)),  '--', 'LineWidth', 1, 'Color', max(tbi_color-0.15,0) )

effect = util_nonparam_effect_size(sham_smpl_x(:,1), tbi_smpl_x(:,1));
effect = round(-effect*100)/100;

title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)

% exportgraphics(gca,[d_fig roi '_' brain_h '_1d_hist_tortuosity.png'],'Resolution',300)


%% Histogram of Gamma0


figure('unit','inch','position',[0 0 2 2])
h1 = histogram(sham_smpl_x(:,2)); hold on
h2 = histogram(tbi_smpl_x(:,2));

h1.Normalization = 'probability'; h1.BinWidth = 0.05;
h1.EdgeColor = sham_color; h1.FaceColor = sham_color; h1.FaceAlpha = 0.85;

h2.Normalization = 'probability'; h2.BinWidth = 0.05; 
h2.EdgeColor = tbi_color; h2.FaceColor = tbi_color; h2.FaceAlpha = 0.85;

xlim([0 2]); pbaspect([2 1 1]);
box off; h = gca; h.YAxis.Visible = 'off';

%%%%%% x-dir is reversed %%%%%%
set(gca,'xtick',[],'ytick',[], 'linewidth',.7, 'xdir', 'reverse')

xline(median(sham_smpl_x(:,2)), '--', 'LineWidth', 1, 'Color', max(sham_color-0.15,0))
xline(median(tbi_smpl_x(:,2)), '--', 'LineWidth', 1, 'Color', max(tbi_color-0.15,0))

effect = util_nonparam_effect_size( exp(sham_smpl_x(:,2)), exp(tbi_smpl_x(:,2))  );
effect = round(-effect*100)/100;

title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)

box off; h = gca; h.YAxis.Visible = 'off'; 
% exportgraphics(gca,[d_fig roi '_' brain_h '_1d_hist_g0.png'],'Resolution',300)


%% Fit linear SVM to sham and tbi data

mdl = fitcsvm(smpl_X(:,1:2), smpl_Y, ...
    'KernelFunction', 'linear', ...
    'PolynomialOrder', [], 'KernelScale', 1, 'BoxConstraint', 1, ...
    'Standardize', false, 'ClassNames', [1, 2]);

mn = min(smpl_X);
mx = max(smpl_X);

beta = mdl.Beta;
bias = mdl.Bias;

slope = -beta(1)/beta(2);
intercept = -bias/beta(2);
slope_perp = -1 / slope;


%%

pause(0.5)

figure('unit','inch','position',[0 0 4 4])
set(groot,'defaultAxesTickLabelInterpreter','default');  

inx = randperm(length(sham_smpl_x));
num_parts = 20;
set_size = floor(numel(inx) / num_parts);
sets = mat2cell(inx', [ones(1, num_parts)*set_size, mod(numel(inx), num_parts)]);

% plot tortuosity and Gamma0 points
for i = 1:length(sets)
    ss = sets{i}; 
    scatter((tbi_smpl_x(ss,1)), (tbi_smpl_x(ss,2)), 15, all_color(tbi_smpl_x(ss,7),:), 'filled', ...
        'MarkerFaceAlpha', 0.8, ...
        'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerEdgeAlpha', 0.3); hold on

    scatter((sham_smpl_x(ss,1)), (sham_smpl_x(ss,2)), 15, all_color(sham_smpl_x(ss,7),:), 'filled', ...
        'MarkerFaceAlpha', 0.8, ...
        'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerEdgeAlpha', 0.3); hold on    
    drawnow
end

% plot SVM hyperplane
xplot = linspace(0, mx(1)+2, 200);
yplot = slope*xplot + intercept;
plot((xplot), (yplot), '--', 'color', svm_cl, 'LineWidth', 2)

% plot the perpendicular line to the SVM hyperplane

% perp_v berings the perpendicular line into the plot view
% you may need to set it manually
perp_v = 4;
yplot_perp = slope_perp*xplot + intercept + perp_v;
plot(xplot, yplot_perp, '--', 'color', proj_cl, 'LineWidth', 2)

set(gca, 'xtick', [1:0.5:3.5], 'ytick',[0:0.5:2.5], 'TickDir','in', 'TickLength', [0,0], ...
    'xaxisLocation','bottom', 'yaxisLocation','left', 'linewidth',1.1, 'fontsize', 17)
grid off; box on; axis equal
xlim([1 3]); ylim([0 2]); pbaspect([1 1 1]);

% exportgraphics(gca,[d_fig roi '_' brain_h '_tortuosity_g0.png'], 'Resolution',300)


%%

% The optimal linear combination zG of the morphological parameters
% is derived from a trained support vector machine

proj_sham = (([1; slope_perp]*[1; slope_perp]') ./ ([1; slope_perp]' * [1; slope_perp]) * sham_smpl_x(:,1:2)')';
proj_tbi = (([1; slope_perp]*[1; slope_perp]') ./ ([1; slope_perp]' * [1; slope_perp]) * tbi_smpl_x(:,1:2)')';

figure('unit','inch','position',[0 0 3 3])
h1 = histogram(proj_sham(:,1)); hold on
h2 = histogram(proj_tbi(:,1));

h1.Normalization = 'probability'; h1.BinWidth = 0.01; 
h1.EdgeColor = sham_color; h1.FaceColor = sham_color; h1.FaceAlpha = 0.85;

h2.BinWidth = 0.01; h2.Normalization = 'probability';
h2.EdgeColor = tbi_color; h2.FaceColor = tbi_color; h2.FaceAlpha = 0.85;

xlim([0.6 1.2]); pbaspect([1 1 1]);
set(gca,'xtick',[0.6 0.9 1.2], 'ytick',[], 'xaxisLocation','bottom', ...
    'xColor', proj_cl, 'YColor', proj_cl, 'linewidth',1.1, 'fontsize', 8)

xline(median(proj_sham(:,1)), '--', 'LineWidth', 1, 'Color', sham_color)
xline(median(proj_tbi(:,1)), '--', 'LineWidth', 1, 'Color', max(tbi_color-0.15,0))

effect = util_nonparam_effect_size( exp(proj_sham(:,1)), exp(proj_tbi(:,1)) );
effect = round(-effect*100)/100;
title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)

% exportgraphics(gca,[d_fig roi '_' brain_h '_hist_projection.png'],'Resolution',300)
