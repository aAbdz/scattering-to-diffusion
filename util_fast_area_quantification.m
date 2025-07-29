function areas = util_fast_area_quantification(obj, res)
z_sum = squeeze(sum(sum(obj, 1), 2));
areas = z_sum * (res^2);
