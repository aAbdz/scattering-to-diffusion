function stats = util_xsection_plane(bw,skel,max_diam,grid_res,vox)

bw = double(bw);
sz_bw = size(bw);
sz_skl = size(skel,1);
inv_res = 1 ./ grid_res;

stats = struct('Point',[],'Area',[]);

unit_tang = normr(diff(skel));
unit_tang = [unit_tang(1,:); unit_tang];

max_xy = max_diam / 2 / vox(1);
[x,y] = ndgrid(-max_xy:grid_res:max_xy, -max_xy:grid_res:max_xy); z = zeros(size(x)); 
xyz = [x(:),y(:),z(:)];

dot_p = unit_tang(:,3);
phi = pi - acos(dot_p);
phi_deg = phi * (180/pi);
rot_axis = cross(unit_tang, repmat([0 0 1], [sz_skl,1]), 2);

for j = 1:sz_skl
    if (phi_deg(j) == 0) || (phi_deg(j)==180)
        XYZnew = xyz';
    else
        [XYZnew, ~, ~] = AxelRot(xyz', phi_deg(j), rot_axis(j,:), [0 0 0]);
    end
    xx = XYZnew(1,:); yy = XYZnew(2,:); zz = XYZnew(3,:);

    point = skel(j,:);
    xx = xx + point(1); 
    yy = yy + point(2); 
    zz = zz + point(3);

    mx = floor(min(max(zz), sz_bw(3)));
    mn = floor(max(min(zz), 1));
    if mn==mx; mx = mn+1; end
    
    cross_section = interp3(bw(:,:,mn:mx),yy,xx,zz-mn+1);
    bw_xsection = cross_section >= 0.5;  
    area = sum(bw_xsection);

    stats(j).Point = point;
    stats(j).Area = vox(1)^2 * area / inv_res^2;
end


