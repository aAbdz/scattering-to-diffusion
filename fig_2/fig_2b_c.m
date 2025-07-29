
clear

d = './data/';
d_save = './fig/';

synth_color = [148 8 125]/255;
sham_color = [40 163 127]/255;
tbi_color = [210 82 24]/255;

% Synthetic axons simulations
load([d, 'simulation_synth.mat']);

volumes = sum(sim.area, 2);
volume_weights = volumes ./ sum(volumes);

AD_av_synth = sum(sim.axial_D .* volume_weights, 1);

%%

% EM axons simulations
load([d, 'simulation_EM.mat']);

% Animal 25 and 49 are sham-operated: group 0
% Animal 2, 24 and 28 are TBI: group 1

ids = sim.ids;
names = {ids.name};
groups = ~contains(names, {'25', '49'});

% Volume weighting sham signal
inx = groups == 0;
volumes = sum(sim.area(inx,:), 2);
volume_weights = volumes ./ sum(volumes);
AD_av_sham = sum(sim.axial_D(inx,:) .* volume_weights, 1);


% Volume weighting tbi signal
inx = groups == 1;
volumes = sum(sim.area(inx,:), 2);
volume_weights = volumes ./ sum(volumes);
AD_av_tbi = sum(sim.axial_D(inx,:) .* volume_weights, 1);

%%

% Intrinsic diffusivity is 2 um^2/ms 
D0 = 2; 

% Linear fit to 1/sqrt(t) for t0>10 ms
t = sim.t;
idx = t>10; 

x = 1./sqrt(t(idx));

y = AD_av_synth(idx);
p_synth = polyfit(x,y,1);

y = AD_av_sham(idx);
p_sham = polyfit(x,y,1);

y = AD_av_tbi(idx);
p_tbi = polyfit(x,y,1);


% Interpolate for t = [1 500] ms 
idx = t>1;
x = 1./sqrt(t(idx));
xval = [0 max(x)];

y_synth = AD_av_synth(idx);
yval_synth = p_synth(1) * xval+p_synth(2);

y_sham = AD_av_sham(idx);
yval_sham = p_sham(1) * xval+p_sham(2);

y_tbi = AD_av_tbi(idx);
yval_tbi = p_tbi(1) * xval+p_tbi(2);


%% Plot Axial D(t) of synth, sham, and TBI axons vs. t

figure('unit','inch','position',[0 0 1.7 2.2]); hold on 
plot(t(2:end),AD_av_synth(2:end)/D0,'color',synth_color,'linewidth',1.2);
plot(t(2:end),AD_av_sham(2:end)/D0,'color',sham_color,'linewidth',1.2); hold on
plot(t(2:end),AD_av_tbi(2:end)/D0,'color',tbi_color,'linewidth',1.2); hold on

xline(10, '--')
xline(100, '--')

xlim([0 500]); ylim([0.75 0.95]); pbaspect([1 1.7 1]);
set(gca,'xtick',[0:250:500],'ytick',0.7:0.1:0.9)
box on; grid off;
xtickangle(0)

exportgraphics(gca,[d_save 'diff_linear.png'],'Resolution',300)

%% Plot Axial D(t) of synth, sham, and TBI axons vs. 1/sqrt(t)

figure('unit','inch','position',[0 0 1.7 2.2]); hold on 
plot(x,y_synth/D0,'color',synth_color,'linewidth',1.2);
plot(xval,yval_synth/D0,'--','color',synth_color,'linewidth',1);

plot(x,y_sham/D0,'color',sham_color,'linewidth',1.2); hold on
plot(xval,yval_sham/D0,'--','color',sham_color,'linewidth',1);

plot(x,y_tbi/D0,'color',tbi_color,'linewidth',1.2); hold on
plot(xval,yval_tbi/D0,'--','color',tbi_color,'linewidth',1);

xline(0.1, '--')
xline(0.3162, '--')

xlim([0 0.5]); ylim([0.75 0.95]); pbaspect([1 1.7 1]);
set(gca,'xtick',[0:0.25:0.5],'ytick',0.7:0.1:0.9)
box on; grid off;
xtickangle(0)

exportgraphics(gca,[d_save 'diff_sqrt.png'],'Resolution',300)