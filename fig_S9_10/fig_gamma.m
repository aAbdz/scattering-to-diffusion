
clear; clc
rng(15)

d_in = '../data/';
d_save = './fig/';

sham_cl = [40 163 127]/255;
tbi_cl = [210 82 24]/255;

roi = 'cg';
brain_h = 'ipsi';

load([d_in roi '_morph.mat']);

step_sz = 2;
voxel_size = [0.05, 0.05, 0.05];
dx = voxel_size(3) * step_sz;
rm = 15;

sham_gamma_dec = zeros(1,3);
sham_bead_pos = cell(1);

tbi_gamma_dec = zeros(1,3);
tbi_bead_pos = cell(1);

count_s = 1;
count_t = 1;

data_name = fieldnames(morph);
for jj = 1:length(data_name)
    if contains(data_name{jj}, brain_h)
        data_ = getfield(morph, data_name{jj});
        sz = size(data_, 2);
        for ii = 1:sz
    
            area_i = data_(ii).area;
            skel_i = data_(ii).skl;
            if sum(area_i) == 0
                continue
            end
            area_i = area_i(rm:end-rm);
            skel_i = skel_i(rm:end-rm,:);

            area_i = util_remove_zeros(area_i); 

            [pos, pos_stats, bead_vol] = util_gamma0_decomposition (area_i', skel_i, step_sz, voxel_size(3));    

            if contains(data_name{jj}, 'h_25') || contains(data_name{jj}, 'h_49') 
                sham_gamma_dec(count_s,1:2) = pos_stats;
                sham_gamma_dec(count_s,3) = bead_vol;

                sham_bead_pos{count_s} = pos;
                count_s = count_s+1;
            else
                tbi_gamma_dec(count_t,1:2) = pos_stats;
                tbi_gamma_dec(count_t,3) = bead_vol; 

                tbi_bead_pos{count_t} = pos;
                count_t = count_t+1;
            end
        end
    end
end


%%

t = isnan(sham_gamma_dec);
[r,~] = find(t(:,1)); 
sham_gamma_dec(r,:) = []; 
sham_bead_pos(r) = [];

t = isnan(tbi_gamma_dec);
[r,~] = find(t(:,1)); 
tbi_gamma_dec(r,:) = []; 
tbi_bead_pos(r) = [];

inx = randperm(size(sham_gamma_dec,1));
sham_gamma_dec = sham_gamma_dec(inx, :);
sham_bead_pos = sham_bead_pos(inx);

inx = randperm(size(tbi_gamma_dec,1));
tbi_gamma_dec = tbi_gamma_dec(inx, :);
tbi_bead_pos = tbi_bead_pos(inx);

sz_ = min(size(sham_gamma_dec,1), size(tbi_gamma_dec,1));

sham_gamma_dec = sham_gamma_dec(1:sz_, :);
tbi_gamma_dec = tbi_gamma_dec(1:sz_, :);

sham_bead_pos = sham_bead_pos(1:sz_);
tbi_bead_pos = tbi_bead_pos(1:sz_);


%% sigma^2 / bar{a}

shamx = sham_gamma_dec(:,2).^2 ./ sham_gamma_dec(:,1);
shamy = (sham_gamma_dec(:,3) ./ sham_gamma_dec(:,1)).^2;

tbix = tbi_gamma_dec(:,2).^2 ./ tbi_gamma_dec(:,1);
tbiy = (tbi_gamma_dec(:,3) ./ sham_gamma_dec(:,1)).^2;

figure('unit','inch','position',[0 0 1.2 1.2])
h1 = histogram(shamx, 'FaceColor', sham_cl); hold on
h2 = histogram(tbix, 'FaceColor', tbi_cl);
h1.Normalization = 'probability';
h1.BinWidth = 0.1; h1.EdgeAlpha = 0.1; h1.FaceAlpha = 1;
h2.Normalization = 'probability';
h2.BinWidth = 0.1; h2.EdgeAlpha = 0.1; h2.FaceAlpha = 0.75;
xlim([0 5]);
set(gca,'xtick',[],'ytick',[],'linewidth',0.5)
box off; h = gca; h.YAxis.Visible = 'off';
xtickangle(0); grid off

xline(median(shamx), '--', 'LineWidth', 1, 'Color', sham_cl-0.15)
xline(median(tbix), '--', 'LineWidth', 1, 'Color', max(tbi_cl-0.15,0))

effect = util_nonparam_effect_size(shamx, tbix);
effect = round(-effect*100)/100;

title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)

% exportgraphics(gca,[d_save roi '_' brain_h '_1d_hist_stats.png'],'Resolution',300)


%% phi^2

figure('unit','inch','position',[0 0 1.2 1.2])
h1 = histogram(shamy, 'FaceColor', sham_cl); hold on
h2 = histogram(tbiy, 'FaceColor', tbi_cl);
h1.Normalization = 'probability';
h1.BinWidth = 1.5; h1.EdgeAlpha = 0.1; h1.FaceAlpha = 1;
h2.Normalization = 'probability';
h2.BinWidth = 1.5; h2.EdgeAlpha = 0.1; h2.FaceAlpha = 0.75;

xlim([0 25]); pbaspect([1 1 1]);
set(gca,'xtick',[],'ytick',[],'linewidth',0.5, 'xdir','reverse')
box off; h = gca; h.YAxis.Visible = 'off';
xtickangle(0); grid off

xline(median(shamy), '--', 'LineWidth',1, 'Color', sham_cl-0.15)
xline(median(tbiy), '--', 'LineWidth', 1, 'Color', max(tbi_cl-0.15,0))

effect = util_nonparam_effect_size(shamy, tbiy);
effect = round(-effect*100)/100;

title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',10)

% exportgraphics(gca,[d_save roi '_' brain_h '_1d_hist_vFrac.png'],'Resolution',300)

pause(0.5)

%%

set(groot,'defaultAxesTickLabelInterpreter','default');  

inx = randperm(length(shamx));
num_parts = 10;
set_size = floor(numel(inx) / num_parts);
sets = mat2cell(inx', [ones(1, num_parts)*set_size, mod(numel(inx), num_parts)]);
figure('unit','inch','position',[0 0 2 2])
for i = 1:length(sets)
    ss = sets{i};  
    scatter(tbix(ss), tbiy(ss), 'SizeData', 15, ...
        'MarkerFaceColor', tbi_cl, 'MarkerFaceAlpha', 0.8, ...
        'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerEdgeAlpha', 0.3); hold on

    scatter(shamx(ss), shamy(ss), 'SizeData', 15, ...
        'MarkerFaceColor', sham_cl, 'MarkerFaceAlpha', 0.8, ...
        'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerEdgeAlpha', 0.3); hold on    
    drawnow
end
xlim([0 5]); ylim([0 25]); 
box on; set(gca, 'xtick', 0:1:5, 'ytick',0:5:60)

set(gca, 'TickDir','in', 'TickLength', [0,0], 'xaxisLocation', 'bottom', 'yaxisLocation','left', 'linewidth',0.7, 'fontsize', 10)
grid off; box on;
xtickangle(0);

% exportgraphics(gca,[d_save roi '_' brain_h '_2d.png'],'Resolution',300)


%% Fig S10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Gamma_0 of bead positions

if contains(roi,'cg')
    L = 40;
else
    L = 70;
end
N = L/dx;
dk = 2*pi/L;
kmax = pi/dx;
kx = 0:dk:2*kmax-dk;
Gammas_sham = zeros(length(sham_bead_pos), N);
for i = 1:length(sham_bead_pos)
    pos = sham_bead_pos{i};
    pos = round(pos/dx);
    m_dist = sham_gamma_dec(i,1);
    if isnan(m_dist)
        continue
    end
    
    dd = pos(end)+round(m_dist/dx);
    if pos(end) > N
        z_pos = zeros(1,dd);
        z_pos(pos) = 1/dx;

        d1 = randi([1 dd-N],1);
        z_pos = z_pos(d1:d1+N-1);

        pos = find(z_pos)*dx;
    else
        pos(pos>N) = [];
        z_pos = zeros(1, N);
        z_pos(pos) = 1 * (1/dx);

        pos = find(z_pos)*dx;
    end
    sham_abar(i) = m_dist;
    pos_fft = fft(z_pos)*dx;
    Gammas_sham(i,:) = abs(pos_fft).^2 / L;   
end
t = isnan(sham_abar);
r = find(t); 
Gammas_sham(r,:) = []; 
sham_abar(r) = [];


Gammas_tbi = zeros(length(tbi_bead_pos), N);
for i = 1:length(tbi_bead_pos)

    pos = tbi_bead_pos{i};
    pos = round(pos/dx);
    m_dist = tbi_gamma_dec(i,1);
    if isnan(m_dist)
        continue
    end
    
    dd = pos(end)+round(m_dist/dx);
    if pos(end) > N
        z_pos = zeros(1,dd);
        z_pos(pos) = 1/dx;

        d1 = randi([1 dd-N],1);
        z_pos = z_pos(d1:d1+N-1);

        pos = find(z_pos)*dx;
    else
        pos(pos>N) = [];
        z_pos = zeros(1, N);
        z_pos(pos) = 1 * (1/dx);

        pos = find(z_pos)*dx;
    end
    tbi_abar(i) = m_dist;
    pos_fft = fft(z_pos)*dx;
    Gammas_tbi(i,:) = abs(pos_fft).^2 / L;   

end
t = isnan(tbi_abar);
r = find(t); 
Gammas_tbi(r,:) = []; 
tbi_abar(r) = [];

ind_max = floor(length(kx)/2);
figure('unit','inch','position',[0 0 2 2])
plot(kx(1:ind_max)/2/pi, mean((sham_abar').*Gammas_sham(:,1:ind_max), 1) , 'linewidth', 1.2, 'color', sham_cl); hold on
plot(kx(1:ind_max)/2/pi, mean((tbi_abar').*Gammas_tbi(:,1:ind_max), 1) , 'linewidth', 1.2, 'color', tbi_cl); hold on
set(gca, 'xscale', 'log', 'yscale', 'log'); hold on; grid on
xlim([0 1e1]); ylim([0.15 1.1]); pbaspect([1 1 1]);
set(gca,'xtick',[0.01, 0.1, 1, 10],'ytick',[0.2:0.2:1], 'xaxisLocation','bottom', 'linewidth', 0.7, 'FontSize',8)

% exportgraphics(gca,[d_save roi '_' brain_h '_gamma1d.png'],'Resolution',300)


%% \bar{a}

shamx = sham_gamma_dec(:,1);
tbix = tbi_gamma_dec(:,1);

figure('unit','inch','position',[0 0 1.2 1.2])
h1 = histogram(shamx, 'FaceColor', sham_cl); hold on
h2 = histogram(tbix, 'FaceColor', tbi_cl);
h1.Normalization = 'probability';
h1.BinWidth = 0.2; h1.EdgeAlpha = 0.1; h1.FaceAlpha = 1;
h2.Normalization = 'probability';
h2.BinWidth = 0.2; h2.EdgeAlpha = 0.1; h2.FaceAlpha = 0.75;

xlim([0 10]); pbaspect([1 1 1]);
set(gca,'xtick',[0:5:10],'ytick',[],'linewidth',0.5)
box on; h = gca; h.YAxis.Visible = 'on';
xtickangle(0); grid off

xline(median(shamx), '--', 'LineWidth',1, 'Color', sham_cl-0.15)
xline(median(tbix), '--', 'LineWidth', 1, 'Color', max(tbi_cl-0.15,0))

effect = util_nonparam_effect_size(shamx, tbix);
effect = round(-effect*100)/100;
title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',8)
% exportgraphics(gca,[d_save roi '_' brain_h '_a-bar_hist1d.png'],'Resolution',300)



%% \sigma / \bar{a}

shamx = sham_gamma_dec(:,2) ./ sham_gamma_dec(:,1);
tbix = tbi_gamma_dec(:,2) ./ tbi_gamma_dec(:,1);

figure('unit','inch','position',[0 0 1.2 1.2])
h1 = histogram(shamx, 'FaceColor', sham_cl); hold on
h2 = histogram(tbix, 'FaceColor', tbi_cl);
h1.Normalization = 'probability';
h1.BinWidth = 0.02; h1.EdgeAlpha = 0.1; h1.FaceAlpha = 1;
h2.Normalization = 'probability';
h2.BinWidth = 0.02; h2.EdgeAlpha = 0.1; h2.FaceAlpha = 0.75;

xlim([0 1]); pbaspect([1 1 1]);
set(gca,'xtick',[0:0.5:1],'ytick',[],'linewidth',0.5)
box on; h = gca; h.YAxis.Visible = 'on';
xtickangle(0); grid off

xline(median(shamx), '--', 'LineWidth',1, 'Color', sham_cl-0.15)
xline(median(tbix), '--', 'LineWidth', 1, 'Color', max(tbi_cl-0.15,0))

effect = util_nonparam_effect_size(shamx, tbix);
effect = round(-effect*100)/100;
title(['$d_{\mathrm{eff}} = \, $', num2str(effect, '%.2f')] ,'interpreter', 'latex','fontsize',8)
% exportgraphics(gca,[d_save roi '_' brain_h '_astd-norm_hist1d.png'],'Resolution',300)




