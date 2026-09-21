function writeROI( path, filename, hheader, basalpoints, externalpoints, nucleuspoints )

basalpoints = uint32(basalpoints);
externalpoints = uint32(externalpoints);
nucleuspoints = uint32(nucleuspoints);

ROI_path=[path,filename,'.roi'];
roi_file=fopen(ROI_path,'w+');

Nb_basal    = uint32( size(basalpoints, 1)      );
Nb_external  = uint32( size(externalpoints, 1)   );
Nb_roi      = uint32( Nb_basal + Nb_external     );
Nb_nuclei   = uint32( size(nucleuspoints, 1)    );

for i=1:numel(hheader)
    header{i} = num2str(hheader(i));
end

%Write Header
fwrite(roi_file, header{1}, 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, header{2}, 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, header{3}, 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, header{4}, 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, header{5}, 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, header{6}, 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, header{7}, 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, header{8}, 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, header{9}, 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, header{10}, 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, header{11}, 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, header{12}, 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, header{13}, 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, header{14}, 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, header{15}, 'char');
fwrite(roi_file, char(32), 'ubit8');

fwrite(roi_file, num2str(Nb_basal), 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, num2str(Nb_external), 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, num2str(Nb_roi), 'char');
fwrite(roi_file, char(32), 'ubit8');
fwrite(roi_file, num2str(Nb_nuclei), 'char');
fwrite(roi_file, char(13), 'ubit8');
fwrite(roi_file, char(10), 'ubit8');

if(Nb_basal > 0)
    for n=1:Nb_basal
        fwrite(roi_file, num2str(basalpoints(n, 1)), 'char');
        fwrite(roi_file, char(32), 'ubit8');
        fwrite(roi_file, num2str(basalpoints(n, 2)), 'char');
        fwrite(roi_file, char(13), 'ubit8');
        fwrite(roi_file, char(10), 'ubit8');
    end
end

if(Nb_external > 0)
    for m=1:Nb_external
        fwrite(roi_file, num2str(externalpoints(m, 1)), 'char');
        fwrite(roi_file, char(32), 'ubit8');
        fwrite(roi_file, num2str(externalpoints(m, 2)), 'char');
        fwrite(roi_file, char(13), 'ubit8');
        fwrite(roi_file, char(10), 'ubit8');
    end
end

if(Nb_nuclei > 0)
    for o=1:Nb_nuclei
        fwrite(roi_file, num2str(nucleuspoints(o, 1)), 'char');
        fwrite(roi_file, char(32), 'ubit8');
        fwrite(roi_file, num2str(nucleuspoints(o, 2)), 'char');
        fwrite(roi_file, char(32), 'ubit8');
        fwrite(roi_file, num2str(nucleuspoints(o, 3)), 'char');
        fwrite(roi_file, char(13), 'ubit8');
        fwrite(roi_file, char(10), 'ubit8');
    end
end

fclose('all');
end