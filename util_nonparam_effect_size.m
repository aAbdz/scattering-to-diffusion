function c_effect = util_nonparam_effect_size(x,y)

nx = length(x);
ny = length(y);

mx = 1.4826 * median(abs(x - median(x)));
my = 1.4826 * median(abs(y - median(y)));

pmad = ((((nx-1) * mx^2) + ((ny-1) * my^2)) ./ (nx+ny-2))^(1/2);

c_effect = (median(x) - median(y)) / pmad;
