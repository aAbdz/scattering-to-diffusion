

% The exact tortuosity limit is validated for both synthetic and EM axons. 

clear

% synthetic axons
d = './data/';
load([d, 'simulation_synth.mat']);
d_save = './fig/';

D0 = 2;
t = sim.t;

% color bar
rgb_points = [
    255, 226, 255;
    255, 21, 255;
    199, 0, 199;
    110, 10, 150;
    52, 0, 54] / 255;

orig_pos = linspace(1, 100, size(rgb_points, 1));
n = 100; new_pos = 1:n;
cmap = [
    interp1(orig_pos, rgb_points(:,1), new_pos);
    interp1(orig_pos, rgb_points(:,2), new_pos);
    interp1(orig_pos, rgb_points(:,3), new_pos)]';

% Error bar reflect errors in estimating D_inf for a varying t0.
num_iter = 10;
t0 = 5; t0_end = t0+10;
t0_range = linspace(t0,t0_end,num_iter);
tortuosity = zeros(1,3,num_iter);
for ii = 1:length(t0_range)
    t0 = t0_range(ii);
    for i = 1:size(sim.area,1)    
        area_i = sim.area(i,:); area_i = area_i(:);
        area_i = util_remove_zeros(area_i);        
        cv = std(area_i)/mean(area_i);
    
        % tortuosity on geometry
        tortuosity_geo = mean(mean(area_i) ./ area_i);
    

        % tortuosity on diffusion
        AD_i = sim.axial_D(i,:);
        t_idx = t > t0;
        x = 1./sqrt(t(t_idx));
        y = AD_i(t_idx);
        p = polyfit(x,y,1);
        d_inf = p(2);
        
        tortuosity_diff = D0 / d_inf;
    
        tortuosity(i,1,ii) = tortuosity_diff;
        tortuosity(i,2,ii) = tortuosity_geo;
        tortuosity(i,3,ii) = cv;
    end
end


% Color bars are defined based on coeff. of variations CV.
lin_cvs = linspace(0,0.5,n);
[~, cl_inx] = min(pdist2(lin_cvs',tortuosity(:,3).^2),[],1);


x = tortuosity(:,1,:); x_av = mean(x,3); x_std = std(x,0,3);
y = tortuosity(:,2,:); y_av = mean(y,3);

figure('unit','inch','position',[0 0 2 2]); hold on
plot([0 3.5], [0 3.5], '--k', 'LineWidth', 1.2);
errorbar(x_av,y_av,x_std,'horizontal', 'LineStyle','none', 'Color','k'); hold on

scatter(x_av, y_av, 15, cmap(cl_inx,:), 'filled', ...
    'MarkerFaceAlpha', 1, ...
    'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerEdgeAlpha', 0.5); 

xlim([1 1.6]); ylim([1 1.6]); pbaspect([1 1 1]);
set(gca,'xtick',[1:0.2:2.5],'ytick',[1:0.2:2.5])
box on; grid on; set(gca,'GridLineStyle','--')

exportgraphics(gca, [d_save 'tortuosity_synth.png'],'Resolution',300)


%%

clear

% EM axons
d = './data/';
load([d, 'simulation_EM.mat']);
d_save = './fig/';

t = sim.t;
D0 = 2;
voxel_size = [0.05, 0.05, 0.05];

rgb_points = [
    200, 255, 200;
    75, 209, 75;
    33, 187, 33;
    20, 122, 20;
    10, 41, 10;
    ] / 255;

orig_pos = linspace(1, 100, size(rgb_points, 1));
n = 100; new_pos = 1:n;
cmap_green = [
    interp1(orig_pos, rgb_points(:,1), new_pos);
    interp1(orig_pos, rgb_points(:,2), new_pos);
    interp1(orig_pos, rgb_points(:,3), new_pos)]';


rgb_points = [
    255, 200, 200;
    209, 75, 75;
    187, 33, 33;
    122, 20, 20;
    41, 10, 10;
    ] / 255;

orig_pos = linspace(1, 100, size(rgb_points, 1));
n = 100; new_pos = 1:n;
cmap_red = [
    interp1(orig_pos, rgb_points(:,1), new_pos);
    interp1(orig_pos, rgb_points(:,2), new_pos);
    interp1(orig_pos, rgb_points(:,3), new_pos)]';

% Error bar reflect errors in estimating D_inf for a varying t0.
num_iter = 10;
t0 = 3; t0_end = t0+10;
t0_range = linspace(t0,t0_end,num_iter);
tortuosity = zeros(1,3,num_iter);
groups = zeros(1);
for ii = 1:length(t0_range)
    t0 = t0_range(ii);
    for i = 1:size(sim.area,1)
        fn = sim.ids(i).name;
        groups(i) = 1;
        if contains(fn,'LM_25') | contains(fn,'LM_49')
            groups(i) = 0;
        end

        area_i = sim.area(i,:); area_i = area_i(:);
        area_i = util_remove_zeros(area_i);
        cv = std(area_i)/mean(area_i);
    
        rr = sqrt(area_i/pi);
        corr_fact = sqrt(2/3) * voxel_size(3);
        rr = rr - corr_fact;
        idx = rr>2*voxel_size(3);
        rr = rr(idx);
        area_i = pi*(rr.^2);
    
        tortuosity_geo = mean(mean(area_i) ./ area_i);
    
        AD_i = sim.axial_D(i,:);
        t_idx = t > t0;
        x = 1./sqrt(t(t_idx));
        y = AD_i(t_idx);
        p = polyfit(x,y,1);
        d_inf = p(2);
        
        tortuosity_diff = D0 / d_inf;
    
        tortuosity(i,1,ii) = tortuosity_diff;
        tortuosity(i,2,ii) = tortuosity_geo;
        tortuosity(i,3,ii) = cv;
    end
end

lin_cvs = linspace(0,0.5,n);
[~, cl_inx] = min(pdist2(lin_cvs',tortuosity(:,3).^2),[],1);

x = tortuosity(:,1,:); x_av = mean(x,3); x_std = std(x,0,3);
y = tortuosity(:,2,:); y_av = mean(y,3);

figure('unit','inch','position',[0 0 2 2]); hold on
plot([0 3], [0 3], '--k', 'LineWidth', 1.2);
errorbar(x_av,y_av,x_std,'horizontal', 'LineStyle','none', 'Color','k'); hold on

inx = groups == 0; 
scatter(x_av(inx), y_av(inx), 15, cmap_green(cl_inx(inx),:), 'filled', ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerEdgeAlpha', 0.5); 

inx = groups == 1; 
scatter(x_av(inx), y_av(inx), 15, cmap_red(cl_inx(inx),:), 'filled', ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeColor', [0.1 0.1 0.1], 'MarkerEdgeAlpha', 0.5); 

xlim([1 2.5]); ylim([1 2.5]); pbaspect([1 1 1]);
set(gca,'xtick',[1:0.5:3],'ytick',[1:0.5:3])
box on; grid on; set(gca,'GridLineStyle','--')


exportgraphics(gca, [d_save 'tortuosity_em.png'],'Resolution',300)



