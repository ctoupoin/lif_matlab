function d_pts = compute_distance_points(c,rows)

num_pts = length(c{1});

ind_corners = false(1,num_pts);
ind_corners([1 rows num_pts-rows+1 num_pts]) = true;

ind_borders = false(1,num_pts);
ind_borders([1:rows rows+1:rows:num_pts-rows rows:rows:num_pts-rows+1 num_pts-rows+1:num_pts]) = true;
ind_borders(ind_corners) = false;

ind_centers = ~ind_borders & ~ind_corners;

coord_cam0 = c{1};

pwd_cam = pdist2([coord_cam0(:,1) coord_cam0(:,2)],...
    [coord_cam0(:,1) coord_cam0(:,2)],'euclidean');

dist_corners_cam0 = sort(pwd_cam0(:,ind_corners));
dist_corners_cam0 = dist_corners_cam0(1:3,:);
dist_corners_cam0 = dist_corners_cam0(dist_corners_cam0 ~= 0);
dist_corners_cam0 = unique(dist_corners_cam0);

dist_borders_cam0 = sort(pwd_cam0(:,ind_borders));
dist_borders_cam0 = dist_borders_cam0(1:4,:);
dist_borders_cam0 = dist_borders_cam0(dist_borders_cam0 ~= 0);
dist_borders_cam0 = unique(dist_borders_cam0);

dist_centers_cam0 = sort(pwd_cam0(:,ind_centers));
dist_centers_cam0 = dist_centers_cam0(1:5,:);
dist_centers_cam0 = dist_centers_cam0(dist_centers_cam0 ~= 0);
dist_centers_cam0 = unique(dist_centers_cam0);

d_pts = vertcat(dist_corners_cam0,dist_borders_cam0,dist_centers_cam0);

end