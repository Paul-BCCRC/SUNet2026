function writeFeatureFile(features, path, filename)
%'cb4'
CB4_path=[path,'\',filename];
cb4_file=fopen(CB4_path,'w');

[numfeat, numCell] = size(features);
numCell = numCell-1;

fwrite(cb4_file,0,'int16');
fwrite(cb4_file,numfeat,'int16');

%Write out feature labels
for i=1:numfeat
    
    flabel = '               ,';
    slabel = numel( features{i, 1} );
    if(slabel <= 15)
        flabel(1:slabel) = features{i, 1};
        for n=1:numel(flabel)
            fwrite(cb4_file, flabel(n), 'uint8');
        end
    else
        disp(['WARNING: Feature label ', features{i, 1}, ' is too long with ', num2str(slabel), ' characters!'])
    end

end

%Write out feature values
for m=1:numCell
    for n=1:numfeat
        fwrite(cb4_file, features{n, m+1}, 'single');
    end
end

fclose('all');

disp('Done.')
end