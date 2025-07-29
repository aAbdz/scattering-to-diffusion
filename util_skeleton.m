function skel = util_skeleton(axon, mov_av)
z = find(sum(sum(axon, 1), 2));
skel = zeros(numel(z), 3);
sz = size(axon);
for jj = z'
    im_2d = axon(:,:,jj);
    [R, C] = ind2sub([sz(1), sz(2)], find(im_2d));
    skel(jj, :) = [mean(R), mean(C), jj];
end
skel = movmean(skel, mov_av);