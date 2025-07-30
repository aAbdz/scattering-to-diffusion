
clear; clc; close all

rng(15)

d_in = '../data/';
d_save = './fig/';

sham_cl = [27 184 99]/255;
tbi_cl = [213 91 36]/255;
proj_cl = [50,50,120]/255;

animal_condition = {'TBI1','TBI2','Sham1','TBI3','Sham2'};
YY = [1, 1, 2, 1, 2];

samples = {'exvivo_250110_AA-1','exvivo_250112_AA-3','exvivo_250113_AA_4','exvivo_250117_AA-5','exvivo_250118_AA_6'};
colors = {[208,88,40]/255, [249,41,43]/255, [23,217,86]/255, [252,144,43]/255, [39,126,83]/255};

region = 'ipsi';
anatomy_labels = {'Spl-CC', 'Bdy-CC', 'Cg'};

diff_t = [7, 15, 20, 30, 40];
D0 = 2; 

load([d_in 'exvivo_AD_' region '.mat']);

param = struct('sample', {}, 'anatomy', {}, 'region', {}, 'tort', [], 'G0', [], 'D_inf', [], 'cD', []);
entry = struct;
for kk = 1:length(anatomy_labels)
    anat = anatomy_labels{kk};
    for ii = 1:length(samples)
        sample = samples{ii};
        
        idx = strcmp({AD.sample}, sample) & ...
              strcmp({AD.anatomy}, anat) & ...
              strcmp({AD.region}, region);

        data = [AD(idx).value];

        tort_ = []; G0_ = []; cD_ = []; D_inf_ = [];
        for jj = 1:size(data,1)
            p = polyfit(1./sqrt(diff_t), data(jj,:), 1); 

            D_inf_(jj) = p(2);
            cD_(jj) = p(1);

            f = 2 * sqrt(D_inf_(jj) / pi);

            tort_(jj) = D0 / D_inf_(jj);
            G0_(jj) = cD_(jj) / f;
        end

        entry.sample = sample;
        entry.anatomy = anat;
        entry.region = region;

        entry.tort = tort_;
        entry.G0 = G0_;
        entry.D_inf = D_inf_;
        entry.cD = cD_;
        param(end+1) = entry;
    end
end


%% G0 and lambda

slope_perps = NaN*ones(15,1);
cc = 1;
for kk = 1:length(anatomy_labels)
    anat = anatomy_labels{kk};
       
    smpl_X1 = []; smpl_X2 = []; smpl_Y = [];

    figure('unit','inch','position',[0 0 4 4]); hold on
    for ii = 1:length(samples)
        cl = colors{ii};
        sample = samples{ii};
       
        idx = strcmp({param.sample}, sample) & ...
              strcmp({param.anatomy}, anat) & ...
              strcmp({param.region}, region);

        cond = [param(idx).cD]>0;

        tort = [param(idx).tort]; tort = tort(cond);
        G0    = [param(idx).G0]; G0 = G0(cond);

        nPoints = length(tort);
     
        smpl_X1 = [smpl_X1 tort];
        smpl_X2 = [smpl_X2 G0];
        smpl_Y = [smpl_Y; repmat(YY(ii), nPoints, 1)];

        errorbar(mean(tort), mean(G0), std(G0), std(G0), std(tort), std(tort), '.', 'Color', cl, 'MarkerSize', 35, 'LineWidth', 1.5); hold on
        
        scatter(tort, G0, 'SizeData', 12, ...
            'MarkerFaceColor', cl, 'MarkerFaceAlpha', 0.7, ...
            'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerEdgeAlpha', 0.3); hold on
    end
    grid off; box on; axis equal
    set(gca, 'xtick', 2.5:0.5:4, 'ytick', 0:0.5:1.5, 'TickDir','in', 'TickLength', [0,0], ...
        'xaxisLocation','bottom', 'yaxisLocation','left', 'linewidth',1.1, 'fontsize', 15)

    xlim([1.5 3.5]); ylim([0 2]);

    smpl_X = [smpl_X1; smpl_X2]';
    mn = min(smpl_X);
    mx = max(smpl_X);

    idx1 = find(smpl_Y == 1);
    idx2 = find(smpl_Y == 2);
    
    % Determine the minimum class count
    n_min = min(length(idx1), length(idx2));
    
    % Randomly sample from both classes
    idx1_sel = randsample(idx1, n_min);
    idx2_sel = randsample(idx2, n_min);
    
    % Combine selected samples
    sel_idx = [idx1_sel; idx2_sel];
    smpl_X = smpl_X(sel_idx, :);
    smpl_Y = smpl_Y(sel_idx);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mdl = fitcsvm(smpl_X, smpl_Y, ...
        'KernelFunction', 'linear', ...
        'KernelScale', 1, ...
        'BoxConstraint', 3, ...
        'Standardize', true, ...
        'ClassNames', [1, 2]);
   
    [x1, x2] = meshgrid(linspace(mn(1), mx(1), 200), ...
                        linspace(mn(2), mx(2), 200));
    xGrid = [x1(:),x2(:)];
    [~, scores] = predict(mdl, xGrid);
    
    C = contour(x1, x2, reshape(scores(:,2),size(x1)), [-10 0 10], '--', 'color', [0,230,240]/255, 'LineWidth',2);

    idx = 1;
    while idx < size(C,2)
        level = C(1,idx);
        n_points = C(2,idx);
        if level == 0
            % Found zero contour
            x0 = C(1, idx+1:idx+n_points);
            y0 = C(2, idx+1:idx+n_points);
            break
        else
            % Skip to next block
            idx = idx + n_points + 1;
        end
    end
    p = polyfit(x0, y0, 1);
    slope_perp = -1 / p(1);

    slope_perps(cc) = slope_perp; cc = cc+1;
    
    x_mid = mean(x0);
    y_mid = mean(y0);
    bias_perp = y_mid - slope_perp * x_mid;
    
    x_vals = linspace(0, max(x0)+5, 100);
    y_vals_perp = slope_perp * x_vals + bias_perp - 0;
    
    plot(x_vals, y_vals_perp, '--', 'color', [50,50,120]/255, 'LineWidth', 1.5)

    % exportgraphics(gcf, [d_save, 'g0-lambda-', anat '-' region, '.png'], 'Resolution', 300)
    title([anat '-' region])    
end
            

%% lambda

for kk = 1:length(anatomy_labels)
    anat = anatomy_labels{kk};
      
    smpl_X1 = []; smpl_X2 = []; smpl_Y = [];
    for ii = 1:length(samples)
        cl = colors{ii};
        sample = samples{ii};
       
        idx = strcmp({param.sample}, sample) & ...
              strcmp({param.anatomy}, anat) & ...
              strcmp({param.region}, region);

        cond = [param(idx).cD]>0;

        tort = [param(idx).tort]; tort = tort(cond);
        smpl_X1 = [smpl_X1 tort];
        
        nPoints = length(tort);
        smpl_Y = [smpl_Y; repmat(YY(ii), nPoints, 1)];
    end

    smpl_X = [smpl_X1; smpl_X2]';

    idx1 = find(smpl_Y == 1);
    idx2 = find(smpl_Y == 2);
    
    % Determine the minimum class count
    n_min = min(length(idx1), length(idx2));
    
    % Randomly sample from both classes
    idx1_sel = randsample(idx1, n_min);
    idx2_sel = randsample(idx2, n_min);
    
    % Combine selected samples
    sel_idx = [idx1_sel; idx2_sel];
    smpl_X = smpl_X(sel_idx, :);
    smpl_Y = smpl_Y(sel_idx);

    smpl_sham = [smpl_X(smpl_Y==2)];
    smpl_tbi  = [smpl_X(smpl_Y==1)];

    figure('unit','inch','position',[0 0 3 3]); hold on
    h1 = histogram(smpl_sham); hold on
    h2 = histogram(smpl_tbi);
    
    h1.Normalization = 'probability'; h1.BinWidth = 0.04; 
    h1.EdgeColor = sham_cl; h1.FaceColor = sham_cl; h1.FaceAlpha = 0.85;
    
    h2.Normalization = 'probability'; h2.BinWidth = 0.04; 
    h2.EdgeColor = tbi_cl; h2.FaceColor = tbi_cl; h2.FaceAlpha = 0.85;

    pbaspect([2 1 1]);
    xlim([1.5 3.5]);
    yl = ylim; ylim([yl(1), yl(2) * 1.05]); 

    set(gca,'xtick',[],'ytick',[], 'xaxisLocation','bottom', 'linewidth',.7)
    box off; h = gca; h.YAxis.Visible = 'off'; 
    
    xline(median(smpl_sham), '--', 'LineWidth', 1, 'Color', max(sham_cl-0.15,0))
    xline(median(smpl_tbi), '--', 'LineWidth', 1, 'Color', max(tbi_cl-0.15,0))
    
    effect = util_nonparam_effect_size(smpl_sham, smpl_tbi);
    effect = round(-effect*100)/100;

    % title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)
    % exportgraphics(gcf, [d_save, 'lambda-hist-', anat '-' region, '.png'], 'Resolution', 300)
    
    title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f') ' \\ ' [anat '-' region]] ,'interpreter', 'latex','fontsize',10)
end


%% G0

for kk = 1:length(anatomy_labels)
    anat = anatomy_labels{kk};
        
    smpl_X1 = []; smpl_X2 = []; smpl_Y = [];
    for ii = 1:length(samples)
        cl = colors{ii};
        sample = samples{ii};
       
        idx = strcmp({param.sample}, sample) & ...
              strcmp({param.anatomy}, anat) & ...
              strcmp({param.region}, region);

        cond = [param(idx).cD]>0;

        G0 = [param(idx).G0]; G0 = G0(cond);
        smpl_X1 = [smpl_X1 G0];
        
        nPoints = length(G0);
        smpl_Y = [smpl_Y; repmat(YY(ii), nPoints, 1)]; 
        
    end

    smpl_X = [smpl_X1; smpl_X2]';

    idx1 = find(smpl_Y == 1);
    idx2 = find(smpl_Y == 2);
    
    % Determine the minimum class count
    n_min = min(length(idx1), length(idx2));
    
    % Randomly sample from both classes
    idx1_sel = randsample(idx1, n_min);
    idx2_sel = randsample(idx2, n_min);
    
    % Combine selected samples
    sel_idx = [idx1_sel; idx2_sel];
    smpl_X = smpl_X(sel_idx, :);
    smpl_Y = smpl_Y(sel_idx);

    smpl_sham = [smpl_X(smpl_Y==2)];
    smpl_tbi  = [smpl_X(smpl_Y==1)];

    figure('unit','inch','position',[0 0 3 3]); hold on
    h1 = histogram(smpl_sham); hold on
    h2 = histogram(smpl_tbi);
    
    h1.Normalization = 'probability'; h1.BinWidth = 0.06; 
    h1.EdgeColor = sham_cl; h1.FaceColor = sham_cl; h1.FaceAlpha = 0.85;
    
    h2.Normalization = 'probability'; h2.BinWidth = 0.06; 
    h2.EdgeColor = tbi_cl; h2.FaceColor = tbi_cl; h2.FaceAlpha = 0.85;

    pbaspect([2 1 1]);
    xlim([0 2]);
    yl = ylim; ylim([yl(1), yl(2) * 1.05]); 

    set(gca,'xtick',[],'ytick',[], 'xaxisLocation','bottom', 'linewidth',.7, 'XDir', 'reverse')
    box off; h = gca; h.YAxis.Visible = 'off'; 
    
    xline(median(smpl_sham), '--', 'LineWidth', 1, 'Color', max(sham_cl-0.15,0))
    xline(median(smpl_tbi), '--', 'LineWidth', 1, 'Color', max(tbi_cl-0.15,0))
    
    effect = util_nonparam_effect_size(smpl_sham, smpl_tbi);
    effect = round(-effect*100)/100;

    % title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)
    % exportgraphics(gcf, [d_save, 'g0-hist-', anat '-' region, '.png'], 'Resolution', 300)

    title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f') ' \\ ' [anat '-' region]] ,'interpreter', 'latex','fontsize',10)
end


%% Zg

cc = 1;
for kk = 1:length(anatomy_labels)
    anat = anatomy_labels{kk};
  
    smpl_X1 = []; smpl_X2 = []; smpl_Y = [];
    slope_perp = slope_perps(cc); cc = cc+1;
  
    for ii = 1:length(samples)
        cl = colors{ii};
        sample = samples{ii};
       
        idx = strcmp({param.sample}, sample) & ...
              strcmp({param.anatomy}, anat) & ...
              strcmp({param.region}, region);

        cond = [param(idx).cD]>0;

        tort = [param(idx).tort]; tort = tort(cond);
        G0    = [param(idx).G0]; G0 = G0(cond);

        nPoints = length(tort);
     
        smpl_X1 = [smpl_X1 tort];
        smpl_X2 = [smpl_X2 G0];
        x = [smpl_X1; smpl_X2];

        nPoints = length(G0);
        smpl_Y = [smpl_Y; repmat(YY(ii), nPoints, 1)];
    end
    
    smpl_X = [smpl_X1; smpl_X2]';

    idx1 = find(smpl_Y == 1);
    idx2 = find(smpl_Y == 2);
    
    % Determine the minimum class count
    n_min = min(length(idx1), length(idx2));
    
    % Randomly sample from both classes
    idx1_sel = randsample(idx1, n_min);
    idx2_sel = randsample(idx2, n_min);
    
    % Combine selected samples
    sel_idx = [idx1_sel; idx2_sel];
    smpl_X = smpl_X(sel_idx, :);
    smpl_Y = smpl_Y(sel_idx);

    smpl_sham = [smpl_X(smpl_Y==2,:)];
    smpl_tbi  = [smpl_X(smpl_Y==1,:)];

    proj_sham = (([1; slope_perp]*[1; slope_perp]') ./ ([1; slope_perp]' * [1; slope_perp]) * smpl_sham')';
    proj_tbi = (([1; slope_perp]*[1; slope_perp]') ./ ([1; slope_perp]' * [1; slope_perp]) * smpl_tbi')';

    figure('unit','inch','position',[0 0 2 2]); hold on
    h1 = histogram(proj_sham(:,1)); hold on
    h2 = histogram(proj_tbi(:,1));
    
    h1.Normalization = 'probability'; h1.BinWidth = 0.06; 
    h1.EdgeColor = sham_cl; h1.FaceColor = sham_cl; h1.FaceAlpha = 0.85;
    
    h2.Normalization = 'probability'; h2.BinWidth = 0.06; 
    h2.EdgeColor = tbi_cl; h2.FaceColor = tbi_cl; h2.FaceAlpha = 0.85;

    set(gca,'xtick',[1 2.5 4], 'ytick',[], 'xaxisLocation','bottom', ...
        'xColor', [50,50,120]/255, 'YColor', [50,50,120]/255, 'linewidth',1.1, 'fontsize', 10)

    xlim([1 4]);
    yl = ylim; ylim([yl(1), yl(2) * 1.05]); 

    xline(median(proj_sham(:,1)), '--', 'LineWidth', 1, 'Color', max(sham_cl-0.15,0))
    xline(median(proj_tbi(:,1)), '--', 'LineWidth', 1, 'Color', max(tbi_cl-0.15,0))
    
    effect = util_nonparam_effect_size(proj_sham(:,1), proj_tbi(:,1));
    effect = round(-effect*100)/100;
    box on

    % title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)
    % exportgraphics(gcf, [d_save, 'zg-hist-', anat '-' region, '.png'], 'Resolution', 300)
    title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f') ' \\ ' [anat '-' region]] ,'interpreter', 'latex','fontsize',10)    
end
