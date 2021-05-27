function [c_matched,r_matched] = points_matching(c,r,c_synth,c_matched,r_matched,num_img)

for ii = 1:num_img
    c_check = zeros(1,length(c_synth{ii}));
    for jj = 1:length(c_synth{ii})
        
        d_match = zeros(1,length(c_synth{ii}));
        for kk = 1:length(c{ii})
            d_match(kk) = sqrt( (c{ii}(kk,1) - c_synth{ii}(jj,1))^2 + (c{ii}(kk,2) - c_synth{ii}(jj,2))^2 );
        end
        
        [~,ind_match] = min(d_match);
        while c_check(ind_match) == 1
            d_match(ind_match) = inf;
            [~,ind_match] = min(d_match);
        end
        c_check(ind_match) = 1;
        
        c_matched{ii}(jj,1) = c{ii}(ind_match,1);
        c_matched{ii}(jj,2) = c{ii}(ind_match,2);
        
        r_matched{ii}(jj) = r{ii}(ind_match,1);
    end
end

end