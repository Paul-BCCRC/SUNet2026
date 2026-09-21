function [ masks ] = bitMap2MaskArray( bitmap )
[Y, X] = size(bitmap);
max_val = max(max(bitmap));
Z=0;

while(max_val >= (2^Z))
    Z=Z+1;
end

masks = zeros(Y, X, Z, 'uint8');

    for z=Z:-1:1
        int_mask = bitand(uint64(bitmap), ones(Y,X,'uint64')*(2^(z-1)));
        logical_mask = logical(int_mask);
        masks(:,:,z) = uint8(logical_mask);
    end

end