function writeCBX(features, path, filename, ext)

CB_path=[path,'\',filename,'.', ext];
cb_file=fopen(CB_path,'w');

[numCell, numfeat] = size(features);
numCell = numCell-1;

fwrite(cb_file,0,'int16');
fwrite(cb_file,numfeat,'int16');

%Write out feature labels
for i=1:numfeat
    
    flabel = '               ,';
    slabel = numel( features{1, i} );
    if(slabel <= 15)
        flabel(1:slabel) = features{1, i};
        for n=1:numel(flabel)
            fwrite(cb_file, flabel(n), 'uint8');
        end
    else
        disp(['WARNING: Feature label ', features{i, 1}, ' is too long with ', num2str(slabel), ' characters!'])
    end

end

%Write out feature values
for m=1:numCell
    for n=1:numfeat
        fwrite(cb_file, features{m+1, n}, 'single');
    end
end

fclose('all');

disp('Done.')
end