
clear; clc; close all

rng(15)

d_in = '../data/';
d_save = './fig/';

animal_condition = {'TBI1','TBI2','Sham1','TBI3','Sham2'};
samples = {'exvivo_250110_AA-1','exvivo_250112_AA-3','exvivo_250113_AA_4','exvivo_250117_AA-5','exvivo_250118_AA_6'};
colors = {[208,88,40]/255, [249,41,43]/255, [23,217,86]/255, [252,144,43]/255, [39,126,83]/255};

% ipsi or contra
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


%% D(t) vs t

for kk = 1:length(anatomy_labels)
    anat = anatomy_labels{kk};
    
    figure('unit','inch','position',[0 0 3 4]); hold on
    for ii = 1:length(samples)
        cl = colors{ii};
        sample = samples{ii};

        idx = strcmp({AD.sample}, sample) & ...
              strcmp({AD.anatomy}, anat) & ...
              strcmp({AD.region}, region);
        data = [AD(idx).value];

        idx = strcmp({param.sample}, sample) & ...
              strcmp({param.anatomy}, anat) & ...
              strcmp({param.region}, region);

        % cD should be >= 0 as a condition
        cond = [param(idx).cD];

        data_c = data(cond>=0,:);
       
        mu = mean(data_c);
        err = std(data_c);

        upper = mu + err;
        lower = mu - err;

        fill([diff_t, fliplr(diff_t)], ...
             [upper, fliplr(lower)],cl,...
             'FaceAlpha', 0.2, 'EdgeColor', 'none'); 

        plot(diff_t, mu, '^-','MarkerSize', 7, 'MarkerFaceColor', cl, 'LineWidth', .7 ,'Color', cl); hold on
    end

    xlim([6 41]); pbaspect([1 2 1])
    set(gca, 'xtick', diff_t, 'ytick', 0.4:0.1:1.5, 'TickDir','in', 'TickLength', [0,0], ...
        'xaxisLocation','bottom', 'yaxisLocation','left', 'linewidth',.7, 'fontsize', 12)
    grid off; box on; xtickangle(0)
    % exportgraphics(gcf, [d_save, 'ad_dti-', anat '-' region, '.png'], 'Resolution', 300)
    title([anat '-' region])
end


%% D(t) vs 1/sqrt(t)


for kk = 1:length(anatomy_labels)
    anat = anatomy_labels{kk};
    
    figure('unit','inch','position',[0 0 3 4]); hold on
    for ii = 1:length(samples)
        cl = colors{ii};
        sample = samples{ii};

        idx = strcmp({AD.sample}, sample) & ...
              strcmp({AD.anatomy}, anat) & ...
              strcmp({AD.region}, region);
        data = [AD(idx).value];

        idx = strcmp({param.sample}, sample) & ...
              strcmp({param.anatomy}, anat) & ...
              strcmp({param.region}, region);

        % cD should be >= 0 as a condition
        cond = [param(idx).cD];

        data_c = data(cond>=0,:);
           
        mu = mean(data_c);
        err = std(data_c);

        upper = mu + err;
        lower = mu - err;
    
        x_ = 1 ./ sqrt(diff_t);
    
        fill([x_, fliplr(x_)], ...
             [upper, fliplr(lower)],cl,...
             'FaceAlpha', 0.2, 'EdgeColor', 'none'); 

        plot(x_, mu, '^','MarkerSize', 7, 'MarkerFaceColor', cl, 'LineWidth', .7 ,'Color', cl); hold on
    
        x_all = repmat(x_, size(data_c,1), 1); x_flat = x_all(:);
        y_flat = data_c(:);
        p = polyfit(x_flat, y_flat, 1);
    
        x_fit = linspace(0.001, max(x_flat), 100);
        y_fit = polyval(p, x_fit);
        plot(x_fit, y_fit, '--', 'color', cl, 'LineWidth', 1.2)
    end
    xlim([0.1 0.38]);
    set(gca, 'xtick', [0.16 0.18 0.22 0.26 0.38], 'ytick', 0.5:0.1:1.3, 'TickDir','in', 'TickLength', [0,0], ...
        'xaxisLocation','bottom', 'yaxisLocation','left', 'linewidth',.7, 'fontsize', 12)
    xtickangle(90); grid off; box on
    % exportgraphics(gcf, [d_save, 'ad_dti-1sqrtt-', anat '-' region, '.png'], 'Resolution', 300)
    title([anat '-' region])
end


