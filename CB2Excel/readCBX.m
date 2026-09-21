function [ features ] = readCBX( path, filename, ext )

CB_path=[path,filename,'.',ext];
cb_file=fopen(CB_path,'r');

cb_data=fread(cb_file, inf, 'uint8');

BYTE_TOTAL = numel(cb_data);
dummy = bytes2type(cb_data( 1:2 ), 'int16');
numFeats = double(bytes2type(cb_data( 3:4 ), 'int16'));

BYTE_FEATS = BYTE_TOTAL - (numFeats*16) - 4;

numCells = floor((BYTE_FEATS/4)/numFeats);

features = cell([numCells+1, numFeats]);

%Get feature labels
for n=1:numFeats
    nn = (n-1)*16;
    features{1,n} = native2unicode(cb_data( nn+5:nn+19 ), 'US-ASCII')';
    features{1,n}  = features{1,n}(find(~isspace(features{1,n} )));%Remove spaces
end

%Get feature values per cell
StartByte = numFeats*16+4+1;
for m=1:numCells
    for n=1:numFeats
        nn = StartByte + (m-1)*(numFeats*4) + (n-1)*4;
        features{m+1, n} = bytes2type(cb_data(nn:nn+3), 'single');
    end
end

    function [returntype] = bytes2type(bytes, type)
        bytes = uint8(bytes);
        returntype = typecast( fliplr(bytes), type);%%flip bytes to change endianness
    end


fclose('all');
end