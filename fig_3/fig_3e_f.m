
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


%% sham and tbi groups

D0 = 2;

% calculate d_inf based on torruosity
sham_dinf = D0 ./ sham_smpl_x(:,1);
tbi_dinf = D0 ./  tbi_smpl_x(:,1);

% calculate c_D based on Gamma0
sham_cd = 2 .* sham_smpl_x(:,2) .* sqrt(sham_dinf ./ pi);
tbi_cd = 2 .* tbi_smpl_x(:,2) .* sqrt(tbi_dinf ./ pi);

% apply sinuosity 
sham_smpl_x = [sham_dinf, sham_cd] ./ (sham_s.^2);
tbi_smpl_x = [tbi_dinf, tbi_cd] ./ (tbi_s.^2);

smpl_X = [sham_smpl_x; tbi_smpl_x];
smpl_W = [sham_w; tbi_w];
smpl_Y = [ones(sz_, 1); 2*ones(sz_, 1)];


%%

mdl = fitcsvm(smpl_X, smpl_Y, 'KernelFunction', 'linear', ...
    'PolynomialOrder', [], 'KernelScale', 1, ...
    'BoxConstraint', 1, 'Weights', smpl_W, ...
    'Standardize', false, 'ClassNames', [1, 2]);

mn = min(smpl_X);
mx = max(smpl_X);

beta = mdl.Beta;
bias = mdl.Bias;

slope = -beta(1)/beta(2);
intercept = -bias/beta(2);

xplot_save = linspace(mn(1), mx(1)+0.5, 200);
yplot_save = slope*xplot_save + intercept;


%% Draw per animal stats over the seperating line of the fitted model

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
wmu_d1 = sum(dinf1 .* w1);
wstd_d1 = util_weighted_std(dinf1, w1);

wmu_d2 = sum(dinf2 .* w2);
wstd_d2 = util_weighted_std(dinf2, w2);
if wmu_d2==0; wmu_d2=[]; wstd_d2=[]; end


wmu_d3 = sum(dinf3 .* w3);
wstd_d3 = util_weighted_std(dinf3, w3);

wmu_d4 = sum(dinf4 .* w4);
wstd_d4 = util_weighted_std(dinf4, w4);

wmu_d5 = sum(dinf5 .* w5);
wstd_d5 = util_weighted_std(dinf5, w5);


% cd weighted
wmu_cd1 = sum(cd1 .* w1);
wstd_cd1 = util_weighted_std(cd1, w1);

wmu_cd2 = sum(cd2 .* w2);
wstd_cd2 = util_weighted_std(cd2, w2);
if wmu_cd2==0; wmu_cd2=[]; wstd_cd2=[]; end

wmu_cd3 = sum(cd3 .* w3);
wstd_cd3 = util_weighted_std(cd3, w3);

wmu_cd4 = sum(cd4 .* w4);
wstd_cd4 = util_weighted_std(cd4, w4);

wmu_cd5 = sum(cd5 .* w5);
wstd_cd5 = util_weighted_std(cd5, w5);


%%

figure('unit','inch','position',[0 0 2 2]); hold on
errorbar(wmu_d1, wmu_cd1, ...
    wstd_cd1, wstd_cd1,wstd_d1, wstd_d1, ...
    'o','MarkerSize',8, 'Color', sham_cl1, ...
    'MarkerFaceColor',sham_cl1, 'LineWidth',1)

errorbar(wmu_d2, wmu_cd2, ...
    wstd_cd2,wstd_cd2,wstd_d2,wstd_d2, ...
    'o','MarkerSize',8, 'Color', sham_cl2, ...
    'MarkerFaceColor',sham_cl2, 'LineWidth',1)

errorbar(wmu_d3, wmu_cd3, ...
    wstd_cd3, wstd_cd3, wstd_d3, wstd_d3, ...
    'o','MarkerSize',8, 'Color', tbi_cl1, ...
    'MarkerFaceColor',tbi_cl1, 'LineWidth',1)

errorbar(wmu_d4, wmu_cd4, ...
    wstd_cd4, wstd_cd4, wstd_d4, wstd_d4, ...
    'o','MarkerSize',8, 'Color', tbi_cl2, ...
    'MarkerFaceColor',tbi_cl2, 'LineWidth',1)

errorbar(wmu_d5, wmu_cd5, ...
    wstd_cd5, wstd_cd5, wstd_d5, wstd_d5, ...
    'o','MarkerSize',8, 'Color', tbi_cl3, ...
    'MarkerFaceColor',tbi_cl3, 'LineWidth',1)

set(gca, 'xtick',0:0.5:2, 'ytick',0:0.5:2)
xlim([1 2]); ylim([0 1]); pbaspect([2 2 1]);

h = gca; h.XAxis.Visible = 'on';

set(gca, 'TickDir','in', 'TickLength', [0,0], 'xaxisLocation','bottom', 'yaxisLocation','left', 'linewidth',0.7, 'fontsize', 9)
grid off; box on; xtickangle(0);
plot(xplot_save, yplot_save, '--', 'color', [0,230,240]/255, 'LineWidth', 1.2)

% exportgraphics(gca,[d_save roi '_' brain_h '_voxel.png'],'Resolution',300)


%%

figure('unit','inch','position',[0 0 2 2]); hold on

% animal 1
xplot = 1./sqrt([8, 10, 20, 70, 100, 1000000]);
yplot = wmu_cd1 * xplot + wmu_d1;

std_y = sqrt((xplot.^2) * wstd_cd1^2 + wstd_d1^2);

yplot_norm = yplot / D0;
std_y_norm = std_y / D0;

x_fill = [xplot, fliplr(xplot)];
y_fill = [yplot_norm + std_y_norm, fliplr(yplot_norm - std_y_norm)];

hold on;
fill(x_fill, y_fill, sham_cl1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(xplot, yplot_norm, 'Color', sham_cl1, 'LineWidth', 1.2);

% animal 2
if ~isempty(wmu_cd2)
    yplot = wmu_cd2 * xplot + wmu_d2;

    std_y = sqrt((xplot.^2) * wstd_cd2^2 + wstd_d2^2);
    
    yplot_norm = yplot / D0;
    std_y_norm = std_y / D0;
    
    x_fill = [xplot, fliplr(xplot)];
    y_fill = [yplot_norm + std_y_norm, fliplr(yplot_norm - std_y_norm)];
    
    hold on;
    fill(x_fill, y_fill, sham_cl2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(xplot, yplot_norm, 'Color', sham_cl2, 'LineWidth', 1.2);
end

% animal 3
yplot = wmu_cd3 * xplot + wmu_d3;
std_y = sqrt((xplot.^2) * wstd_cd3^2 + wstd_d3^2);

yplot_norm = yplot / D0;
std_y_norm = std_y / D0;

x_fill = [xplot, fliplr(xplot)];
y_fill = [yplot_norm + std_y_norm, fliplr(yplot_norm - std_y_norm)];
hold on;
fill(x_fill, y_fill, tbi_cl1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(xplot, yplot_norm, 'Color', tbi_cl1, 'LineWidth', 1.2);


% animal 4
yplot = wmu_cd4 * xplot + wmu_d4;
std_y = sqrt((xplot.^2) * wstd_cd4^2 + wstd_d4^2);

yplot_norm = yplot / D0;
std_y_norm = std_y / D0;

x_fill = [xplot, fliplr(xplot)];
y_fill = [yplot_norm + std_y_norm, fliplr(yplot_norm - std_y_norm)];
hold on;
fill(x_fill, y_fill, tbi_cl2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(xplot, yplot_norm, 'Color', tbi_cl2, 'LineWidth', 1.2);

% animal 5
yplot = wmu_cd5 * xplot + wmu_d5;
std_y = sqrt((xplot.^2) * wstd_cd5^2 + wstd_d5^2);
yplot_norm = yplot / D0;
std_y_norm = std_y / D0;
x_fill = [xplot, fliplr(xplot)];
y_fill = [yplot_norm + std_y_norm, fliplr(yplot_norm - std_y_norm)];
hold on;
fill(x_fill, y_fill, tbi_cl3, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(xplot, yplot_norm, 'Color', tbi_cl3, 'LineWidth', 1.2);

xlim([0 0.35]); ylim([0.4 1.12]); pbaspect([2 2 1]);
set(gca,'TickDir','in', 'TickLength', [0,0], 'xtick',[0,0.175,0.35],'ytick',0:0.2:1.2, 'linewidth',0.7, 'fontsize', 9)

box on; grid off; xtickangle(0)
xline(0.1,'--','Alpha',0.5)
xline(0.3162,'--','Alpha',0.5)
box on; grid off;
xtickangle(0)

% exportgraphics(gca,[d_save 'sim_diff_linear_err.png'],'Resolution',300)



