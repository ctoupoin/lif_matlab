function [d_px,O,inclination] = get_target_infos(calib_pattern,c,r)

num_img = length(calib_pattern);
d_px = zeros(1,num_img);
inclination = zeros(1,num_img);
O = zeros(num_img,2);
% O1 = zeros(num_img,2);

% ScrSize=get(0,'MonitorPositions');
figure('Position',[100 100 100 1000]),
disp('Please click on the upper left point, then on the upper right one.')
for ii = 1:num_img
    title([num2str(ii) '/' num2str(num_img)])
    imshow(-calib_pattern{ii},[])
    [x,y] = ginput(3);
    
    num_detected_pts = length(c{ii});
    
    dist1 = zeros(1,num_detected_pts);
    dist2 = zeros(1,num_detected_pts);
    dist3 = zeros(1,num_detected_pts);
    for jj = 1:num_detected_pts
        dist1(jj) =  sqrt((c{ii}(jj,1) - x(1)).^2 + (c{ii}(jj,2) - y(1)).^2);
        dist2(jj) =  sqrt((c{ii}(jj,1) - x(2)).^2 + (c{ii}(jj,2) - y(2)).^2);
        dist3(jj) =  sqrt((c{ii}(jj,1) - x(3)).^2 + (c{ii}(jj,2) - y(3)).^2);
    end
    [~,ind_pt1] = min(dist1);
    [~,ind_pt2] = min(dist2);
    [~,ind_pt3] = min(dist3);
    
    viscircles(c{ii}(ind_pt1,:),r{ii}(ind_pt1,:),'Color','b');
    viscircles(c{ii}(ind_pt2,:),r{ii}(ind_pt2,:),'Color','r');
    viscircles(c{ii}(ind_pt3,:),r{ii}(ind_pt3,:),'Color','g');
    
    pause
    
    d_px(ii) = sqrt( (c{ii}(ind_pt1,1) - c{ii}(ind_pt2,1))^2 + (c{ii}(ind_pt1,2) - c{ii}(ind_pt2,2))^2 );
    inclination(ii) = atan2(c{ii}(ind_pt1,2) - c{ii}(ind_pt2,2),c{ii}(ind_pt1,1) - c{ii}(ind_pt2,1));
    O(ii,:) = c{ii}(ind_pt3,:);
%     O1(ii,:) = c{ii}(ind_pt3,:);
%     O(ii,:) = [mean(c{ii}(:,1)) mean(c{ii}(:,2))];
end
close

end