function [d_px,O,inclination] = get_target_infos_v2(calib_pattern,c,r)

num_img = length(calib_pattern);
d_px = zeros(1,num_img);
inclination = zeros(1,num_img);
O = zeros(num_img,2);

ScrSize=get(0,'MonitorPositions');
figure('Position',ScrSize),
disp('Please click on all 4 corners points, starting from the upper left and going clockwise')
for ii = 1:num_img
    title([num2str(ii) '/' num2str(num_img)])
    imshow(calib_pattern{ii},[])
    [x,y] = ginput(4);
    
    num_detected_pts = length(c{ii});
    
    dist1 = zeros(1,num_detected_pts);
    dist2 = zeros(1,num_detected_pts);
    dist3 = zeros(1,num_detected_pts);
    dist4 = zeros(1,num_detected_pts);
    for jj = 1:num_detected_pts
        dist1(jj) =  sqrt((c{ii}(jj,1) - x(1)).^2 + (c{ii}(jj,2) - y(1)).^2);
        dist2(jj) =  sqrt((c{ii}(jj,1) - x(2)).^2 + (c{ii}(jj,2) - y(2)).^2);
        dist3(jj) =  sqrt((c{ii}(jj,1) - x(3)).^2 + (c{ii}(jj,2) - y(3)).^2);
        dist4(jj) =  sqrt((c{ii}(jj,1) - x(4)).^2 + (c{ii}(jj,2) - y(4)).^2);
    end
    [~,ind_pt1] = min(dist1);
    [~,ind_pt2] = min(dist2);
    [~,ind_pt3] = min(dist3);
    [~,ind_pt4] = min(dist4);
    
    viscircles(c{ii}(ind_pt1,:),r{ii}(ind_pt1,:),'Color','b');
    viscircles(c{ii}(ind_pt2,:),r{ii}(ind_pt2,:),'Color','r');
    viscircles(c{ii}(ind_pt3,:),r{ii}(ind_pt3,:),'Color','g');
    viscircles(c{ii}(ind_pt4,:),r{ii}(ind_pt4,:),'Color','y');
    
    pause
    
    d_px_top = sqrt( (c{ii}(ind_pt1,1) - c{ii}(ind_pt2,1))^2 + (c{ii}(ind_pt1,2) - c{ii}(ind_pt2,2))^2 );
    d_px_bottom = sqrt( (c{ii}(ind_pt4,1) - c{ii}(ind_pt3,1))^2 + (c{ii}(ind_pt4,2) - c{ii}(ind_pt3,2))^2 );
    
    inclination_top = atan2(c{ii}(ind_pt1,2) - c{ii}(ind_pt2,2),c{ii}(ind_pt1,1) - c{ii}(ind_pt2,1));
%     inclination_bottom = atan2(c{ii}(ind_pt4,2) - c{ii}(ind_pt3,2),c{ii}(ind_pt4,1) - c{ii}(ind_pt3,1));
    
    d_px(ii) = mean([d_px_top, d_px_bottom]);
    inclination(ii) = inclination_top;
    
    O(ii,:) = c{ii}(ind_pt1,:);
end
close

end