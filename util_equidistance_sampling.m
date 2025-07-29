function eq_sampled_skel = util_equidistance_sampling(skel, fineness, step_sz)

z = size(skel,1);
n = fineness * z;
q = linspace(1,z,n);

xq = interp1(1:z, skel(:,1), q, 'spline');
yq = interp1(1:z, skel(:,2), q, 'spline');
zq = interp1(1:z, skel(:,3), q, 'spline');
intp_skel = [xq;yq;zq]';

d = sqrt(sum(diff(intp_skel).^2, 2));
cum_d = cumsum(d);
eq_d = (0:step_sz:max(cum_d))';

inx = zeros(length(eq_d),1);
for i = 1:length(eq_d)
    p = eq_d(i);
    [~,inx(i)] = min(abs(cum_d - p));
end
eq_sampled_skel = intp_skel(inx, :);


