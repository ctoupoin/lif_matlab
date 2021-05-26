function roi_c = intensity_corrections_v8(roi,transfer_function,threshold,filter_size)

[size_x, size_y] = size(roi);
padding = padarray(roi,[size_x,size_y],0,'both');


% ROI
fft_roi = fft2(padding);

transformed_roi = fftshift(fft_roi .* transfer_function);

amp_im = log(abs(transformed_roi));
max_amp = max(max(amp_im));
threshold_amp = threshold .* max_amp;

spikes = amp_im > threshold_amp;
spikes(round(1.4*size_x:1.6*size_x),:) = 0;
spikes(:,round(1.4*size_y:1.6*size_y)) = 0;

[r,c] = find(spikes);

% figure,imshow(spikes,[])

for ii = 1:length(r)
    r_span = max(r(ii)-filter_size,1):min(r(ii)+filter_size,3*size_x);
    c_span = max(c(ii)-filter_size,1):min(c(ii)+filter_size,3*size_y);
    
    transformed_roi(r(ii),c(ii)) = mean(mean(transformed_roi(r_span,c_span)));
end
transformed_roi = ifftshift(transformed_roi);

roi_c = ifft2( transformed_roi);
roi_c = real(roi_c(size_x+1:2*size_x,size_y+1:2*size_y));

end