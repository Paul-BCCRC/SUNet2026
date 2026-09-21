function [ basal_roi, external_roi, nuclei_centers, header, P, PP ] = readROI( path, filename )

ROI_path=[path,filename];
roi_file=fopen(ROI_path,'r');

roi_data=fread(roi_file, inf, 'uint8');
N = numel(roi_data);

tstr = '';
p=1;
for n=1:N
    
    if(roi_data(n) == 13) % 13 = CR Carriage Return
        Points{p} = tstr;
        tstr = '';
        p=p+1;
    elseif(roi_data(n) == 10) % 10 = LF Line Feed
        
    else
        tstr = [tstr, char(roi_data(n))];
    end
    
    
end

Points = Points';

M = numel(Points);
%pcount =1;
for m=1:M
    C = strsplit(Points{m}, ' ');
    PP{m} = C;
%     if(numel(C)==2)
%         P(pcount,1) = str2num(C{1});
%         P(pcount,2) = str2num(C{2});
%         pcount=pcount+1;
%     end
end
% P(pcount ,:) = P(1,:);
PP = PP';

for i=1:(numel(PP{1})-1)
    header(i) = str2num(PP{1}{i});
end
n_basal = header(16);
n_external = header(17);
n_total = header(18);
try
    n_nuclei = header(19);
catch
    n_nuclei = 0;
end

basal_roi=[];
for n=1:n_basal
    basal_roi(n, 1) = str2num( PP{1+n}{1} );
    basal_roi(n, 2) = str2num( PP{1+n}{2} );
end

external_roi=[];
for n=1:n_external
    external_roi(n, 1) = str2num( PP{1+n_basal+n}{1} );
    external_roi(n, 2) = str2num( PP{1+n_basal+n}{2} );
end

nuclei_centers=[];
for n=1:n_nuclei
    nuclei_centers(n, 1) = str2num( PP{1+n_total+n}{1} );
    nuclei_centers(n, 2) = str2num( PP{1+n_total+n}{2} );
    nuclei_centers(n, 3) = str2num( PP{1+n_total+n}{3} );
end


if(n_basal > 0 && n_external > 0)
    P = zeros(n_total+1, 2);
    P(1:n_basal, :) = basal_roi;
    P((n_basal+1):n_total , :) = flip(external_roi);
    P(n_total+1, :) = basal_roi(1,:);
    
elseif(n_basal > 0 && n_external <= 0)
    
    P = basal_roi;
    P(end+1, :) = basal_roi(1, :);
    
elseif(n_basal <= 0 && n_external > 0)
    
    P = external_roi;
    P(end+1, :) = external_roi(1, :);

else
    P = 0;
end

fclose('all')

end