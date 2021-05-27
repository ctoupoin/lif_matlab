function c_synth = synthetic_grid_creation(inclination,d_mm,gamma_cam,O,num_img,columns,rows)

c_synth = cell(1,num_img);

for ii = 1:num_img
    count = 0;
    d_pts_px0 = d_mm .* gamma_cam(ii);
    rotation_angle = pi + inclination(ii);
    
    R = [cos(rotation_angle), sin(rotation_angle);...
    -sin(rotation_angle), cos(rotation_angle)];
    
    for jj = 1:columns
        for kk = 1:rows
            count = count + 1;
            
            c_synth{ii}(count,1) = (jj - 1) * d_pts_px0;
            c_synth{ii}(count,2) = (kk - 1) * d_pts_px0;
        end
    end
    num_pts = length(c_synth{ii});
    
    c_synth{ii} = c_synth{ii} * R;
    c_synth{ii}(:,1) = c_synth{ii}(:,1) + O(ii,1) - c_synth{ii}(round(num_pts/2),1); 
    c_synth{ii}(:,2) = c_synth{ii}(:,2) + O(ii,2) - c_synth{ii}(round(num_pts/2),2);
end

end