function [loc_maxs, pos_stats, bead_vol] = util_gamma0_decomposition(area_i, skel_i, step_sz, vox)

N = length(area_i);
L = sum(sqrt(sum(diff(skel_i.*vox).^2, 2)));
dx = step_sz * vox; %nm

alpha = area_i / mean(area_i);
eta = log(alpha);
filt_eta = smoothdata(eta, 'gaussian', 20);

filt_area = smoothdata(area_i, 'gaussian', 20);
filt_area(filt_area<min(area_i)) = min(area_i);
a0 = min(filt_area);


min_d = 2/dx;   % minimum of 2um dist between beads
min_w = 0.5/dx; % minimum of 0.5um bead width
[~, loc_mins] = findpeaks(-filt_area, 'MinPeakDistance', min_d, 'MinPeakWidth', min_w); 
loc_mins = [1, loc_mins, N];

% Compute relative area drop in each interval
% Filter out intervals where drop is too small >>>> not a bead
paired_ind = [loc_mins(1:end-1)', loc_mins(2:end)'];
for i = 1:size(paired_ind,1)
    a = filt_area(paired_ind(i,1):paired_ind(i,2));
    a1 = a - a0;
    v1s(i) = sum(a1) / sum(a);
end
paired_ind(v1s < 0.2, :) = [];


loc_maxs = [];
for i = 1:size(paired_ind,1)
    a = filt_area(paired_ind(i,1):paired_ind(i,2));
    [~,ind] = max(a);
    loc_maxs(i) = paired_ind(i,1) + ind - 1;
end
loc_maxs = loc_maxs*dx;

abar = mean(diff(loc_maxs));
astd = std(diff(loc_maxs));

e0 = min(filt_eta);
bead_v = sum(filt_eta - e0) * dx;
av_bead_v = bead_v / abar;

pos_stats = [abar, astd]; 
bead_vol = av_bead_v;    

if length(loc_maxs)<3
    pos_stats = [NaN NaN]; 
    bead_vol = NaN;
    loc_maxs = [];
end

