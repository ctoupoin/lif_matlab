function [bkgrd1,bkgrd2] = background_creation_v2(paths,params_cam1,params_cam2,mask_common_roi,tform)

% Changelog :
%   - v2 : included the if loop which checks whether or not the averaged
%   background images exist


if exist(paths.bkgrd_no_laser,'dir')
    
    if exist(paths.bkgrd_no_laser_averaged_cam1,'file') || exist(paths.bkgrd_no_laser_averaged_cam2,'file')
        bkgrd1 = double(imread(paths.bkgrd_no_laser_averaged_cam1));
        bkgrd2 = double(imread(paths.bkgrd_no_laser_averaged_cam2));
    else
        imgs_bkgrd = dir(paths.bkgrd_no_laser);
        format = imgs_bkgrd(3).name(end-2:end);
        imgs_bkgrd = dir([paths.bkgrd_no_laser '*.' format]);
        imgs_names = vertcat(imgs_bkgrd.name);
        
        ind_cam1 = imgs_names(:,31) == '1';
        ind_cam1 = ind_cam1';
        ind_cam2 = ~ind_cam1;
        
        ind_cam1 = find(ind_cam1);
        ind_cam2 = find(ind_cam2);
        
        num_img_bkgrd = length(ind_cam1);
        
        length_name = length([imgs_bkgrd(end).folder '/' imgs_bkgrd(end).name]);
        
        imgs_bkgrd_name1 = repmat('a',length(ind_cam1),length_name);
        imgs_bkgrd_name2 = repmat('a',length(ind_cam2),length_name);
        for ii = ind_cam1
            imgs_bkgrd_name1(ii,:) = [imgs_bkgrd(ii).folder '/' imgs_bkgrd(ii).name];
        end
        count2 = 0;
        for ii = ind_cam2
            count2 = count2 + 1;
            imgs_bkgrd_name2(count2,:) = [imgs_bkgrd(ii).folder '/' imgs_bkgrd(ii).name];
        end
        
        count_loop = 0;
        for ii = 1:min(20,num_img_bkgrd)
            count_loop = count_loop + 1;
            
            [roi1,roi2] = geometric_transforms(imgs_bkgrd_name1,imgs_bkgrd_name2,params_cam1,...
                params_cam2,mask_common_roi,tform,ii);
            
            if ii == 1
                bkgrd1 = roi1;
                bkgrd2 = roi2;
            else
                bkgrd1 = count_loop ./ (count_loop + 1) .* bkgrd1 + 1 ./ (count_loop + 1) .* roi1;
                bkgrd2 = count_loop ./ (count_loop + 1) .* bkgrd2 + 1 ./ (count_loop + 1) .* roi2;
            end
        end

        imwrite(uint16(bkgrd1),paths.bkgrd_no_laser_averaged_cam1);
        imwrite(uint16(bkgrd2),paths.bkgrd_no_laser_averaged_cam2);
    end
    
    
    if all(size(bkgrd1) ~= [mask_common_roi(2)-mask_common_roi(1)+1 mask_common_roi(4)-mask_common_roi(3)+1])
        delete(paths.bkgrd_no_laser_averaged_cam1);
        delete(paths.bkgrd_no_laser_averaged_cam2);
        
        
        imgs_bkgrd = dir(paths.bkgrd_no_laser);
        format = imgs_bkgrd(3).name(end-2:end);
        imgs_bkgrd = dir([paths.bkgrd_no_laser '*.' format]);
        imgs_names = vertcat(imgs_bkgrd.name);
        
        ind_cam1 = imgs_names(:,31) == '1';
        ind_cam1 = ind_cam1';
        ind_cam2 = ~ind_cam1;
        
        ind_cam1 = find(ind_cam1);
        ind_cam2 = find(ind_cam2);
        
        num_img_bkgrd = length(ind_cam1);
        
        length_name = length([imgs_bkgrd(end).folder '/' imgs_bkgrd(end).name]);
        
        imgs_bkgrd_name1 = repmat('a',length(ind_cam1),length_name);
        imgs_bkgrd_name2 = repmat('a',length(ind_cam2),length_name);
        for ii = ind_cam1
            imgs_bkgrd_name1(ii,:) = [imgs_bkgrd(ii).folder '/' imgs_bkgrd(ii).name];
        end
        count2 = 0;
        for ii = ind_cam2
            count2 = count2 + 1;
            imgs_bkgrd_name2(count2,:) = [imgs_bkgrd(ii).folder '/' imgs_bkgrd(ii).name];
        end
        
        count_loop = 0;
        for ii = 1:min(20,num_img_bkgrd)
            count_loop = count_loop + 1;
            
            [roi1,roi2] = geometric_transforms(imgs_bkgrd_name1,imgs_bkgrd_name2,params_cam1,...
                params_cam2,mask_common_roi,tform,ii);
            
            if ii == 1
                bkgrd1 = roi1;
                bkgrd2 = roi2;
            else
                bkgrd1 = count_loop ./ (count_loop + 1) .* bkgrd1 + 1 ./ (count_loop + 1) .* roi1;
                bkgrd2 = count_loop ./ (count_loop + 1) .* bkgrd2 + 1 ./ (count_loop + 1) .* roi2;
            end
        end
        
        imwrite(uint16(bkgrd1),paths.bkgrd_no_laser_averaged_cam1);
        imwrite(uint16(bkgrd2),paths.bkgrd_no_laser_averaged_cam2);
    end
else
    error('No background images with no laser')
end