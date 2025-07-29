
clear

rng(15)

d_in = '../data/';
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
sham_smpl_x = sham_x(1:sz_, :);
tbi_smpl_x = tbi_x(1:sz_, :);

% total volume of all shams
sham_v = sum(sham_smpl_x(:,4));

% total volume of all TBIs
tbi_v = sum(tbi_smpl_x(:,4));

% volume weight of axons
sham_w = sham_smpl_x(:,4) ./ sum(sham_v);
tbi_w = tbi_smpl_x(:,4) ./ sum(tbi_v);

% sinuosity 
sham_s = sham_smpl_x(:,3);
tbi_s = tbi_smpl_x(:,3);


%% converts sham and tbi groups

D0 = 2;

% calculate d_inf based on torruosity
sham_dinf = D0 ./ sham_smpl_x(:,1);
tbi_dinf = D0 ./  tbi_smpl_x(:,1);

% calculate c_D based on Gamma0
sham_cd = 2 .* sham_smpl_x(:,2) .* sqrt(sham_dinf ./ pi);
tbi_cd = 2 .* tbi_smpl_x(:,2) .* sqrt(tbi_dinf ./ pi);

% apply sinuosity 
sham_ = [sham_dinf, sham_cd] ./ (sham_s.^2);
tbi_ = [tbi_dinf, tbi_cd] ./ (tbi_s.^2);

ws_rp = ceil(1e3 * (sham_w / min(sham_w)));
sham_smpl_x_rp = repelem(sham_(:,1), ws_rp); 
sham_smpl_y_rp = repelem(sham_(:,2), ws_rp); 

wt_rp = ceil(1e3 * (tbi_w / min(tbi_w)));
tbi_smpl_x_rp = repelem(tbi_(:,1), wt_rp);
tbi_smpl_y_rp = repelem(tbi_(:,2), wt_rp);

smpl_X = [sham_; tbi_];
smpl_W = [sham_w; tbi_w];
smpl_Y = [ones(sz_, 1); 2*ones(sz_, 1)];


%% convert for individual animals

s_x1 = sham_smpl_x(sham_smpl_x(:,7)==1, :);
s_x2 = sham_smpl_x(sham_smpl_x(:,7)==2, :);
s_x3 = tbi_smpl_x(tbi_smpl_x(:,7)==3, :);
s_x4 = tbi_smpl_x(tbi_smpl_x(:,7)==4, :);
s_x5 = tbi_smpl_x(tbi_smpl_x(:,7)==5, :);

% calculate d_inf based on torruosity
dinf1 = D0 ./ s_x1(:,1);
dinf2 = D0 ./ s_x2(:,1);
dinf3 = D0 ./ s_x3(:,1);
dinf4 = D0 ./ s_x4(:,1);
dinf5 = D0 ./ s_x5(:,1);

% calculate c_D based on Gamma0
cd1 = 2 .* s_x1(:,2) .* sqrt(dinf1 ./ pi);
cd2 = 2 .* s_x2(:,2) .* sqrt(dinf2 ./ pi);
cd3 = 2 .* s_x3(:,2) .* sqrt(dinf3 ./ pi);
cd4 = 2 .* s_x4(:,2) .* sqrt(dinf4 ./ pi);
cd5 = 2 .* s_x5(:,2) .* sqrt(dinf5 ./ pi);


%% apply sinuosity 

dinf1 = dinf1 ./ (s_x1(:,3).^2);
dinf2 = dinf2 ./ (s_x2(:,3).^2);
dinf3 = dinf3 ./ (s_x3(:,3).^2);
dinf4 = dinf4 ./ (s_x4(:,3).^2);
dinf5 = dinf5 ./ (s_x5(:,3).^2);

cd1 = cd1 ./ (s_x1(:,3).^2);
cd2 = cd2 ./ (s_x2(:,3).^2);
cd3 = cd3 ./ (s_x3(:,3).^2);
cd4 = cd4 ./ (s_x4(:,3).^2);
cd5 = cd5 ./ (s_x5(:,3).^2);

%% Volume weighting

v1 = sum(s_x1(:,4));
v2 = sum(s_x2(:,4));
v3 = sum(s_x3(:,4));
v4 = sum(s_x4(:,4));
v5 = sum(s_x5(:,4));

w1 = s_x1(:,4)./v1;
w2 = s_x2(:,4)./v2;
w3 = s_x3(:,4)./v3;
w4 = s_x4(:,4)./v4;
w5 = s_x5(:,4)./v5;

%%

w1_rp = ceil(1e3 * (w1 / max(w1)));
w2_rp = ceil(1e3 * (w2 / max(w2)));
w3_rp = ceil(1e3 * (w3 / max(w3)));
w4_rp = ceil(1e3 * (w4 / max(w4)));
w5_rp = ceil(1e3 * (w5 / max(w5)));

dinf1_rp = repelem(dinf1, w1_rp); 
if isempty(w2_rp)
    dinf2_rp = [];
else
    dinf2_rp = repelem(dinf2, w2_rp); 
end
dinf3_rp = repelem(dinf3, w3_rp); 
dinf4_rp = repelem(dinf4, w4_rp); 
dinf5_rp = repelem(dinf5, w5_rp); 

cd1_rp = repelem(cd1, w1_rp);
if isempty(w2_rp)
    cd2_rp = [];
else
    cd2_rp = repelem(cd2, w2_rp); 
end
cd3_rp = repelem(cd3, w3_rp); 
cd4_rp = repelem(cd4, w4_rp); 
cd5_rp = repelem(cd5, w5_rp); 


%% Histogram of d_inf

figure('unit','inch','position',[0 0 2 2])
h1 = histogram([dinf1_rp; dinf2_rp]); hold on
h2 = histogram([dinf3_rp; dinf4_rp; dinf5_rp]);

h1.Normalization = 'probability'; h1.BinWidth = 0.04; 
h1.EdgeColor = sham_color; h1.FaceColor = sham_color; h1.FaceAlpha = 0.85;

h2.Normalization = 'probability'; h2.BinWidth = 0.04; 
h2.EdgeColor = tbi_color; h2.FaceColor = tbi_color; h2.FaceAlpha = 0.85;

xlim([0.5 2]); pbaspect([2 1 1]);
box off; h = gca; h.YAxis.Visible = 'off';
set(gca,'xtick',[],'ytick',[],'linewidth', .7)

xline(median(sham_smpl_x_rp), '--', 'LineWidth', 1, 'Color', max(sham_color-0.15,0))
xline(median(tbi_smpl_x_rp), '--', 'LineWidth', 1, 'Color', max(tbi_color-0.15,0))

effect = util_nonparam_effect_size( sham_smpl_x_rp, tbi_smpl_x_rp  );
effect = round(-effect*100)/100;

title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)

% exportgraphics(gca,[d_fig roi '_' brain_h '_1d_hist_dinf.png'],'Resolution',300)


%% Histogram of c_D

figure('unit','inch','position',[0 0 2 2])
h1 = histogram([cd1_rp; cd2_rp]); hold on
h2 = histogram([cd3_rp; cd4_rp; cd5_rp]);

h1.Normalization = 'probability'; h1.BinWidth = 0.08; 
h1.EdgeColor = sham_color; h1.FaceColor = sham_color; h1.FaceAlpha = 0.85;

h2.Normalization = 'probability'; h2.BinWidth = 0.08; h2.FaceAlpha = 0.85;
h2.EdgeColor = tbi_color; h2.EdgeColor = tbi_color;

xlim([0 3]); pbaspect([2 1 1]);
box off; h = gca; h.YAxis.Visible = 'off';
set(gca,'xtick',[],'ytick',[],'linewidth', .7,'xdir', 'reverse')

xline(median(sham_smpl_y_rp), '--', 'LineWidth', 1, 'Color', max(sham_color-0.15,0))
xline(median(tbi_smpl_y_rp), '--', 'LineWidth', 1, 'Color', max(tbi_color-0.15,0))

effect = util_nonparam_effect_size( sham_smpl_y_rp, tbi_smpl_y_rp  );
effect = round(-effect*100)/100;

title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)

% exportgraphics(gca,[d_fig roi '_' brain_h '_1d_hist_cd.png'],'Resolution',300)

pause(0.5)

%% Fit SVM

mdl = fitcsvm(smpl_X, smpl_Y, 'KernelFunction', 'linear', ...
    'PolynomialOrder', [], 'KernelScale', 1, ...
    'BoxConstraint', 1, 'Weights', smpl_W, ...
    'Standardize', false, 'ClassNames', [1, 2]);

mn = min(smpl_X);
mx = max(smpl_X);

[x1, x2] = meshgrid(linspace(mn(1), mx(1), 200), ...
                    linspace(mn(2), mx(2), 200));
xGrid = [x1(:),x2(:)];
[~, scores] = predict(mdl, xGrid);

beta = mdl.Beta;
bias = mdl.Bias;

slope = -beta(1)/beta(2);
intercept = -bias/beta(2);
slope_perp = -1 / slope;


%% plot data with decision boundary

figure('unit','inch','position',[0 0 2 2])

inx = randperm(length(sham_smpl_x));
num_parts = 10;
set_size = floor(numel(inx) / num_parts);
sets = mat2cell(inx', [ones(1, num_parts)*set_size, mod(numel(inx), num_parts)]);

for i = 1:length(sets)
    ss = sets{i};  
    scatter(tbi_(ss,1), tbi_(ss,2), 10000*tbi_w(ss), all_color(tbi_smpl_x(ss,7),:), 'filled', ...
        'MarkerFaceAlpha', 0.8, ...
        'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerEdgeAlpha', 0.3); hold on

    scatter(sham_(ss,1), sham_(ss,2), 10000*sham_w(ss), all_color(sham_smpl_x(ss,7),:), 'filled', ...
        'MarkerFaceAlpha', 0.8, ...
        'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerEdgeAlpha', 0.3); hold on    
    drawnow
end
box on;
set(gca, 'xtick',0:0.5:2, 'ytick',0:0.5:9)
xlim([0 2]); ylim([0 2]); pbaspect([2 2 1]);

xplot = linspace(mn(1)-1, mx(1)+0.5, 200);
yplot = slope*xplot + intercept;
plot(xplot, yplot, '--', 'color', [0,230,240]/255, 'LineWidth', 1.2)

% perp_v berings the perpendicular line into the plot view
% you may need to set it manually
perp_v = -3.5;
yplot_perp = slope_perp*xplot + intercept + perp_v;
plot(xplot, yplot_perp, '--', 'color', bbox_cl, 'LineWidth', 1.2)
set(groot,'defaultAxesTickLabelInterpreter','default');  
set(gca, 'TickDir','in', 'TickLength', [0,0], 'xaxisLocation','bottom', 'yaxisLocation','left', 'linewidth',0.7, 'fontsize', 9)
box on;
xtickangle(0);

% exportgraphics(gca,[d_fig roi '_' brain_h '_2d.png'],'Resolution',300)


%%

% The optimal linear combination zD of the diffusional parameters
% is derived from a trained support vector machine

sham1 = [dinf1 cd1];
sham2 = [dinf2 cd2];
tbi1 = [dinf3 cd3];
tbi2 = [dinf4 cd4];
tbi3 = [dinf5 cd5];

proj1 = (([1; slope_perp]*[1; slope_perp]') ./ ([1; slope_perp]' * [1; slope_perp]) * sham1')';
proj2 = (([1; slope_perp]*[1; slope_perp]') ./ ([1; slope_perp]' * [1; slope_perp]) * sham2')';
proj3 = (([1; slope_perp]*[1; slope_perp]') ./ ([1; slope_perp]' * [1; slope_perp]) * tbi1')';
proj4 = (([1; slope_perp]*[1; slope_perp]') ./ ([1; slope_perp]' * [1; slope_perp]) * tbi2')';
proj5 = (([1; slope_perp]*[1; slope_perp]') ./ ([1; slope_perp]' * [1; slope_perp]) * tbi3')';

w1_rp = ceil(1e3 * (w1 / min(w1)));
w2_rp = ceil(1e3 * (w2 / min(w2)));
w3_rp = ceil(1e3 * (w3 / min(w3)));
w4_rp = ceil(1e3 * (w4 / min(w4)));
w5_rp = ceil(1e3 * (w5 / min(w5)));

proj1_rp = repelem(proj1(:,1), w1_rp); 
if isempty(w2_rp)
    proj2_rp = [];
else
    proj2_rp = repelem(proj2(:,1), w2_rp); 
end
proj3_rp = repelem(proj3(:,1), w3_rp); 
proj4_rp = repelem(proj4(:,1), w4_rp); 
proj5_rp = repelem(proj5(:,1), w5_rp); 

figure('unit','inch','position',[0 0 3 3])
if isempty(w2_rp)
    h1 = histogram(proj1_rp(:,1)); hold on
else
    h1 = histogram([proj1_rp(:,1); proj2_rp(:,1)]); hold on
end
h2 = histogram([proj3_rp(:,1); proj4_rp(:,1); proj5_rp(:,1)]);

h1.Normalization = 'probability'; h1.BinWidth = 0.04; 
h1.EdgeColor = sham_color; h1.FaceColor = sham_color; h1.FaceAlpha = 0.85;

h2.Normalization = 'probability'; h2.BinWidth = 0.04; 
h2.EdgeColor = tbi_color; h2.FaceColor = tbi_color; h2.FaceAlpha = 0.85;

xlim([0.5 2.2]); pbaspect([1 1 1]);

set(gca,'xtick',[0.5 1.35 2.2], 'ytick',[], 'xaxisLocation','bottom', ...
    'xColor', proj_cl, 'YColor', proj_cl, 'linewidth',1.1, 'fontsize', 10)

proj_sham = (([1; slope_perp]*[1; slope_perp]') ./ ([1; slope_perp]' * [1; slope_perp]) * sham_')';
proj_tbi = (([1; slope_perp]*[1; slope_perp]') ./ ([1; slope_perp]' * [1; slope_perp]) * tbi_')';

proj_sham_rp = repelem(proj_sham(:,1), ws_rp); 
proj_tbi_rp = repelem(proj_tbi(:,1), wt_rp);

xline(median(proj_sham_rp), '--', 'LineWidth', 1, 'Color', max(sham_color-0.15,0))
xline(median(proj_tbi_rp), '--', 'LineWidth', 1, 'Color', max(tbi_color-0.15,0))

effect = util_nonparam_effect_size( proj_sham_rp, proj_tbi_rp  );
effect = round(-effect*100)/100;

title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)

% exportgraphics(gca,[d_fig roi '_' brain_h '_diff_projection.png'],'Resolution',300)











