

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
sham_smpl_x = sham_x(1:sz_, :);
tbi_smpl_x = tbi_x(1:sz_, :);

smpl_X = [sham_smpl_x; tbi_smpl_x];
smpl_Y = [ones(sz_, 1); 2*ones(sz_, 1)];


%%

s_x1 = sham_smpl_x(sham_smpl_x(:,7)==1, :);
s_x2 = sham_smpl_x(sham_smpl_x(:,7)==2, :);
s_x3 = tbi_smpl_x(tbi_smpl_x(:,7)==3, :);
s_x4 = tbi_smpl_x(tbi_smpl_x(:,7)==4, :);
s_x5 = tbi_smpl_x(tbi_smpl_x(:,7)==5, :);

%%

mdl = fitcsvm(smpl_X(:,1:2), smpl_Y, 'KernelFunction', 'linear', ...
    'PolynomialOrder', [], 'KernelScale', 1, 'BoxConstraint', 1, ...
    'Standardize', false, 'ClassNames', [1, 2]);

mn = min(smpl_X);
mx = max(smpl_X);

beta = mdl.Beta;
bias = mdl.Bias;

xplot_save = linspace(mn(1), mx(1), 200);
yplot_save = (-beta(1)*xplot_save - bias) / beta(2);


%%

[r,~] = find(x1<=0); x1(r,:) = [];
[r,~] = find(x2<=0); x2(r,:) = []; 
[r,~] = find(x3<=0); x3(r,:) = [];
[r,~] = find(x4<=0); x4(r,:) = [];
[r,~] = find(x5<=0); x5(r,:) = [];

v1 = sum(x1(:,4));
v2 = sum(x2(:,4));
v3 = sum(x3(:,4));
v4 = sum(x4(:,4));
v5 = sum(x5(:,4));

w1 = x1(:,4)./v1;
w2 = x2(:,4)./v2;
w3 = x3(:,4)./v3;
w4 = x4(:,4)./v4;
w5 = x5(:,4)./v5;


%% converts

D0 = 2;

dinf1 = D0 ./ x1(:,1);
dinf2 = D0 ./ x2(:,1);
dinf3 = D0 ./ x3(:,1);
dinf4 = D0 ./ x4(:,1);
dinf5 = D0 ./ x5(:,1);

cd1 = 2 .* x1(:,2) .* sqrt(dinf1 ./ pi);
cd2 = 2 .* x2(:,2) .* sqrt(dinf2 ./ pi);
cd3 = 2 .* x3(:,2) .* sqrt(dinf3 ./ pi);
cd4 = 2 .* x4(:,2) .* sqrt(dinf4 ./ pi);
cd5 = 2 .* x5(:,2) .* sqrt(dinf5 ./ pi);

%% apply sinuosity 

dinf1 = dinf1 ./ (x1(:,3).^2);
dinf2 = dinf2 ./ (x2(:,3).^2);
dinf3 = dinf3 ./ (x3(:,3).^2);
dinf4 = dinf4 ./ (x4(:,3).^2);
dinf5 = dinf5 ./ (x5(:,3).^2);

cd1 = cd1 ./ (x1(:,3).^2);
cd2 = cd2 ./ (x2(:,3).^2);
cd3 = cd3 ./ (x3(:,3).^2);
cd4 = cd4 ./ (x4(:,3).^2);
cd5 = cd5 ./ (x5(:,3).^2);


%%

% d_inf weighted
sham_cl = [40 163 127]/255;
tbi_cl = [210 82 24]/255;

w_dinf1 = sum(dinf1 .* w1);
std_dinf1 = util_weighted_std(dinf1, w1);

w_dinf2 = sum(dinf2 .* w2);
std_dinf2 = util_weighted_std(dinf2, w2);

w_dinf3 = sum(dinf3 .* w3);
std_dinf3 = util_weighted_std(dinf3, w3);

w_dinf4 = sum(dinf4 .* w4);
std_dinf4 = util_weighted_std(dinf4, w4);

w_dinf5 = sum(dinf5 .* w5);
std_dinf5 = util_weighted_std(dinf5, w5);


% cd weighted
w_cd1 = sum(cd1 .* w1);
std_cd1 = util_weighted_std(cd1, w1);

w_cd2 = sum(cd2 .* w2);
std_cd2 = util_weighted_std(cd2, w2);

w_cd3 = sum(cd3 .* w3);
std_cd3 = util_weighted_std(cd3, w3);

w_cd4 = sum(cd4 .* w4);
std_cd4 = util_weighted_std(cd4, w4);

w_cd5 = sum(cd5 .* w5);
std_cd5 = util_weighted_std(cd5, w5);


%%

tortuosity1 = D0 / w_dinf1;
tortuosity2 = D0 / w_dinf2;
if tortuosity2==inf; tortuosity2=[]; end
tortuosity3 = D0 / w_dinf3;
tortuosity4 = D0 / w_dinf4;
tortuosity5 = D0 / w_dinf5;

err_tort1 = std_dinf1 / (w_dinf1^2);
err_tort2 = std_dinf2 / (w_dinf2^2);
if isnan(err_tort2); err_tort2=[]; end
err_tort3 = std_dinf3 / (w_dinf3^2);
err_tort4 = std_dinf4 / (w_dinf4^2);
err_tort5 = std_dinf5 / (w_dinf5^2);

gamma0_1 = w_cd1 ./ (2*sqrt(w_dinf1 ./ pi));
gamma0_2 = w_cd2 ./ (2*sqrt(w_dinf2 ./ pi));
if isnan(gamma0_2); gamma0_2=[]; end
gamma0_3 = w_cd3 ./ (2*sqrt(w_dinf3 ./ pi));
gamma0_4 = w_cd4 ./ (2*sqrt(w_dinf4 ./ pi));
gamma0_5 = w_cd5 ./ (2*sqrt(w_dinf5 ./ pi));

err_g1 = (std_cd1 / 4 / (sqrt(w_dinf1/pi))) + ((w_cd1 * std_dinf1)/ 8 / sqrt((w_dinf1^3)/pi) );
err_g2 = (std_cd2 / 4 / (sqrt(w_dinf2/pi))) + ((w_cd2 * std_dinf2)/ 8 / sqrt((w_dinf2^3)/pi) );
if isnan(err_g2); err_g2=[]; end
err_g3 = (std_cd3 / 4 / (sqrt(w_dinf3/pi))) + ((w_cd3 * std_dinf3)/ 8 / sqrt((w_dinf3^3)/pi) );
err_g4 = (std_cd4 / 4 / (sqrt(w_dinf4/pi))) + ((w_cd4 * std_dinf4)/ 8 / sqrt((w_dinf4^3)/pi) );
err_g5 = (std_cd5 / 4 / (sqrt(w_dinf5/pi))) + ((w_cd5 * std_dinf5)/ 8 / sqrt((w_dinf5^3)/pi) );


%%

figure('unit','inch','position',[0 0 2 2]) 
errorbar(tortuosity1, gamma0_1, ...
    err_g1, err_g1, err_tort1, err_tort1, ...
    'o','MarkerSize',7, 'Color', sham_cl1, ...
    'MarkerFaceColor',sham_cl1, 'LineWidth',1); hold on

errorbar(tortuosity2, gamma0_2, ...
    err_g2, err_g2, err_tort2, err_tort2, ...
    'o','MarkerSize',7, 'Color', sham_cl2, ...
    'MarkerFaceColor',sham_cl2, 'LineWidth',1); hold on


errorbar(tortuosity3, gamma0_3, ...
    err_g3, err_g3, err_tort3, err_tort3, ...
    'o','MarkerSize',7, 'Color', tbi_cl1, ...
    'MarkerFaceColor',tbi_cl1, 'LineWidth',1); hold on

errorbar(tortuosity4, gamma0_4, ...
    err_g4, err_g4, err_tort4, err_tort4, ...
    'o','MarkerSize',7, 'Color', tbi_cl2, ...
    'MarkerFaceColor',tbi_cl2, 'LineWidth',1); hold on

errorbar(tortuosity5, gamma0_5, ...
    err_g5, err_g5, err_tort5, err_tort5, ...
    'o','MarkerSize',7, 'Color', tbi_cl3, ...
    'MarkerFaceColor',tbi_cl3, 'LineWidth',1); hold on

xlim([1 2]); ylim([0 1]);
set(gca, 'xtick', [1:0.5:3], 'ytick', [0:0.5:2], 'TickDir','in', 'TickLength', [0,0], ...
    'xaxisLocation','bottom', 'yaxisLocation','left', 'linewidth',.7, 'fontsize', 9)

grid off; box on
plot(xplot_save, yplot_save, '--', 'color', [0,230,240]/255, 'LineWidth', 1.2)

% exportgraphics(gca,[d_save roi '_' brain_h '_inverse_model.png'],'Resolution',300)

