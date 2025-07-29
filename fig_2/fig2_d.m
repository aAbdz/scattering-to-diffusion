
% Coarse-graining over the increasing diffusion length ell(t) makes an axon
% appear increasingly more uniform.

% To illustrate the effect, an axon segment is Gaussian-filtered
% with the standard deviation ell(t) / sqrt{2}
% for ell(t) = 0, 5, 10, 20 um.

% We artificially stitched the 1d axonal Area signal multiple times 
% by concatenating each axon with its flipped version, thereby generating longer synthetic axons.

clear;  clc

d_in = './data/cg_LM_49_contra_axon_1.mat';
d_save = './figs/';

colors = [150 150 170; 
          110 110 150;
          55  55  100;
          20  20  80] / 255;

D0 = 2; % um^2/ms
dims = 1;
voxel_size = 0.05;

% EM axon
bin_obj = load(d_in).axon;
bin_obj = imdilate(bin_obj, true(3,3,7));
area_i = util_fast_area_quantification(bin_obj, voxel_size); 
volume = sum(bin_obj(:));

% 4 time points
ts = [0, 6.88, 27.53, 110.13];

area_i = [area_i; flip(area_i); area_i; flip(area_i); area_i; flip(area_i)];

% To increase randomness, and break the priodicity,
% we removed a segment from the beginning and end of the Area after stitching.
area_i = area_i(360:end-360);
area_i = util_remove_zeros(area_i);

N = length(area_i);
dx = voxel_size;
L = N*dx;
dq = 2*pi/L;
kmax = pi/dx;
q = (0:dq:2*kmax-dq)';
h_inx = floor(length(q)/2);

% D_inf is estimated based on geometry
D_inf = D0 ./ mean(mean(area_i)./area_i);
ell = round(sqrt(2*D_inf*dims*ts));
beta = 0.93;

figure('unit','inch','position',[0 0 2.7 2.7]); hold on
for i = 1:length(ell)
    
    w = ell(i) / 2;  
    aq = fft(area_i)*dx;
    filtq = fftshift(exp(-((q-kmax).^2) * (w^2)));
    filt_aq = aq .* filtq;
    filt_area = real(ifft(filt_aq) * (1/dx));

    [Gamma, Gamma_0, ~, xval, yval] = util_gamma0(filt_area, L, dx, beta);

    plot(q(1:h_inx)/2/pi, Gamma, 'linewidth',1.2, 'Color', colors(i,:)); hold on
    plot(xval, yval, '--', 'LineWidth', 1, 'Color', colors(i,:)); 
    
    xlim([1.9*10^-3 10])
    ylim([10^-6 3])


    set(gca, 'xscale', 'log', 'yscale', 'log'); hold on;

    xticks(10.^[-2:1]); yticks(10.^[-6, -3:1:1]); grid off; box on
    xline(10.^[-2:1], ':'); yline(10.^[-6, -3:1:1], ':')

    fn = [d_save 'synth_gamma' num2str(i) '.png'];
    % exportgraphics(gca,fn,'Resolution',300)

end


%%

% To illustrate the effect of coarse-graining on 3d axonal geometry,
% an axon segment is Gaussian-filtered
% with the standard deviation ell(t) / sqrt{2} for ell(t) = 0, 5, 10, 20 um.
% Because the axon segment is a 3d geometry, the transverse dimensions
% are also Gaussian filtered with a samll increasing width.

clc

dxy = [1, 7, 13, 17]; 
ell = sqrt(2*D_inf*dims*ts) * 1/voxel_size;
dz = ell / sqrt(2);

% As the Gaussian filter width increases in all dimensions,
% the filtered object becomes more diluted, requiring a lower threshold. 
ths = [0.45, 0.2, 0.11, 0.07]; 

% EM axon will be Gaussian filtered to visualize the effect of coarse graining
% pad with zero in-plane 
% pad with replication along the z-axis
pad_v = 400;
obj_padded = single(padarray(bin_obj, [50,50,0]));
obj_padded = single(padarray(obj_padded, [0,0,pad_v], 'replicate'));

for i = 1:length(ell)
    dx = dxy(i); dx = max(1, dx);
    dz_i = dz(i); dz_i = max(1, dz_i);
    th = ths(i);

    obj_gauss = imgaussfilt3(obj_padded, [dx,dx,dz_i]);
    obj_gauss = obj_gauss(:,:,pad_v+1:end-pad_v); % remove padding. 
    bin_gauss = obj_gauss>th;
   
    inx1 = 500; inx2 = 1000;  % Cutting 25 um of the axon for visualization
    figure('unit','inch','position',[0 0 3 3])
    p = patch(isosurface(bin_gauss(:,:,inx1:inx2),0.5));
    p.FaceColor = colors(i,:); p.EdgeColor = 'none'; p.FaceAlpha = 1;
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    grid off; axis off; axis equal; material dull; lighting gouraud; camlight; camlight
    view(180, 0) 

    fn = [d_save num2str(ts(i)) '.png'];
    % exportgraphics(gca,fn,'Resolution','300')
end
