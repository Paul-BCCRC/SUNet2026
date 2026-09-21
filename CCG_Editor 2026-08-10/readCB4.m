function [ features ] = readCB4( path, filename )

CB4_path=[path,filename,'.cb4'];
cb4_file=fopen(CB4_path,'r');

cb4_data=fread(cb4_file, inf, 'uint8');

BYTE_TOTAL = numel(cb4_data);
dummy = bytes2type(cb4_data( 1:2 ), 'int16');
numFeats = double(bytes2type(cb4_data( 3:4 ), 'int16'));

BYTE_FEATS = BYTE_TOTAL - (numFeats*16) - 4;

numCells = floor((BYTE_FEATS/4)/numFeats);

features = cell([numFeats, numCells+1]);

%Get feature labels
for n=1:numFeats
    nn = (n-1)*16;
    features{n,1} = native2unicode(cb4_data( nn+5:nn+19 ), 'US-ASCII');
    features{n,1}  = features{n,1}(find(~isspace(features{n,1} )));%Remove spaces
end

%Get feature values per cell
StartByte = numFeats*16+4+1;
for m=1:numCells
    for n=1:numFeats
        nn = StartByte + (m-1)*(numFeats*4) + (n-1)*4;
        features{n,m+1} = bytes2type(cb4_data(nn:nn+3), 'single');
    end
end

    function [returntype] = bytes2type(bytes, type)
        bytes = uint8(bytes);
        returntype = typecast( fliplr(bytes), type);%%flip bytes to change endianness
    end


fclose('all');
end