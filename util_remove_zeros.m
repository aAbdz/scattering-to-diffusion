function area_i = util_remove_zeros(area_i)
l = numel(area_i);
th = max(0.1*mean(area_i), 0.01);

if area_i(1)<=th
    area_i(1) = th;
end
if area_i(l)<=th
    area_i(l) = th;
end

z_inx = find(area_i<=th);
z_inx(z_inx==1) = [];
z_inx(z_inx==l) = [];
if ~isempty(z_inx)
    for ii = 1:length(z_inx)
        z = z_inx(ii);
        next_nonzero = find(area_i(z+1:end)>th,1,'first');
        next_val = area_i(next_nonzero+z);
        area_i(z) = mean([area_i(z-1),next_val]);
    end
end

av = movmean(area_i, 17);
z_inx = find(abs(area_i-av)/mean(area_i)>0.3);
area_i(z_inx) = av(z_inx);

% z_inx = find(abs(diff(area_i)/mean(area_i))>0.3)+1;
% while ~isempty(z_inx)
%     z = z_inx(1);
%     area_i(z) = area_i(z-1);  
%     z_inx = find(abs(diff(area_i))>0.3)+1;
% end
