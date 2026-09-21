function [ bitmap, nuclei_mask, overlap_mask ] = maskArray2BitMap( masks )

[Y, X, Z] = size(masks);

%if(Z<32)
    
bitmap = zeros(Y, X, 'single');
nuclei_mask = zeros(Y, X, 'uint8');
overlap_mask = zeros(Y, X, 'uint8');
for z=1:Z
    logical_mask = logical(masks(:,:,z));
    int_mask = single(logical_mask) .* (2^(z-1));
    bitmap = bitmap + int_mask;
end

masksum = sum(double(masks), 3);
overlap_id = find(masksum >=2);%%%%%%%%%%%%%%%%%%%OVERLAP
nuclei_id = find(masksum >=1);%%%%%%%%%%%%%%%%%%NUCLEI
nuclei_mask(nuclei_id) = 255;
overlap_mask(overlap_id) = 255;
%else
%    disp('ROI > 32...')
%end
end