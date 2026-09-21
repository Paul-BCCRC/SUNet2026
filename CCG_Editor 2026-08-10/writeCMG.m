function [  ] = writeCMG( header, Images, Masks, path, filename )
if( ~strcmp(path(end), '\') )
    path = [path,'\'];
end
cmg_time = tic;

disp('-------------------------------------------------')
disp([path, filename, '.ccg'])

CMG_path=[path,'\',filename,'.cmg'];
cmg_file=fopen(CMG_path,'w');
I = numel(Images);


prc_dis=0;
prc_ =  floor(double(1)*100/double(I));
fprintf(1,'Writing CMG...%3d%%\n',prc_);

for i=1:I
    
    prc_ =  floor(double(i)*100/double(I));
    if(prc_ > prc_dis)
        fprintf(1,'\b\b\b\b%3.0f%%',prc_);
        prc_dis = prc_;
    end

    NbColorMap=size(Images{i},3);
    NbBitMap=size(Masks{i},3);
    
    I_RGB=Images{i};
	I_Bitmap=Masks{i};
	I_Bitmap(I_Bitmap==255)=1;

    
%     if (NbBitMap==3)
%         Background_Bitmap=uint8(~I_Bitmap);
%         BackgroundFloatRGB=I_RGB.*I_Bitmap;
%         BackgroundFloatR=BackgroundFloatRGB(:,:,header.RGB_Order{i}(1)+1);
%         BackgroundFloatG=BackgroundFloatRGB(:,:,header.RGB_Order{i}(2)+1);
%         BackgroundFloatB=BackgroundFloatRGB(:,:,header.RGB_Order{i}(3)+1);
%         MeanR=mean(BackgroundFloatR(BackgroundFloatR~=0));
%         MeanG=mean(BackgroundFloatG(BackgroundFloatG~=0));
%         MeanB=mean(BackgroundFloatB(BackgroundFloatB~=0));
%         BackgroundFloat=header.RedFaction(i).*MeanR+header.GreenFaction(i).*MeanG+header.BlueFaction(i).*MeanB;
%     else
%         Background_Bitmap=uint8(~I_Bitmap);
%         BackgroundFloatGrayscale=I_RGB.*I_Bitmap;
%         BackgroundFloat=mean(BackgroundFloatGrayscale(BackgroundFloatGrayscale~=0));
%     end
    
    fwrite(cmg_file,'c','ubit8');
    fwrite(cmg_file,header.Mode(i),'ubit8');
    fwrite(cmg_file,NbColorMap,'ubit8');
    fwrite(cmg_file,header.Class(i),'ubit32');
    fwrite(cmg_file,header.Screenx(i),'ubit32');
    fwrite(cmg_file,header.Screeny(i),'ubit32');
    fwrite(cmg_file,header.Stagex(i),'ubit64');
    fwrite(cmg_file,header.Stagey(i),'ubit64');
    fwrite(cmg_file,header.Stagez(i),'ubit64');
    fwrite(cmg_file,header.Resolution(i),'single');
    fwrite(cmg_file,header.LowThreshold(i),'ubit16');
    fwrite(cmg_file,header.MidThreshold(i),'ubit16');
    fwrite(cmg_file,header.Group(i),'ubit8');
    fwrite(cmg_file,size(I_RGB,2),'ubit32');%width of the image
    fwrite(cmg_file,size(I_RGB,1),'ubit32');%height of the image
    fwrite(cmg_file,header.Accession(i),'ubit32');
    fwrite(cmg_file,header.Iod(i),'single');
    fwrite(cmg_file,header.Fluor(i),'ubit8');
    fwrite(cmg_file,header.Diagnosis(i),'ubit16');
    fwrite(cmg_file,header.RedFaction(i),'single');
    fwrite(cmg_file,header.GreenFaction(i),'single');
    fwrite(cmg_file,header.BlueFaction(i),'single');
    fwrite(cmg_file,header.Index(i),'ubit32');
    fwrite(cmg_file,header.Objective(i),'ubit32');
    fwrite(cmg_file,header.Calibrated(i),'ubit8');
    fwrite(cmg_file,header.StackX_int(i),'ubit32');
    fwrite(cmg_file,header.StackY_int(i),'ubit32');
    fwrite(cmg_file,NbBitMap,'ubit8');
    fwrite(cmg_file,header.CassettePosition(i),'ubit8');
    fwrite(cmg_file,header.vorx(i),'ubit32');
    fwrite(cmg_file,header.vory(i),'ubit32');
    fwrite(cmg_file,header.BestFocusFrame(i),'ubit8');
    fwrite(cmg_file,header.BackgroundFloat(i),'single');
    fwrite(cmg_file,header.Layer(i),'ubit16');
    % fwrite(cmg_file,header.Layer{i}(1),'ubit8');
    % fwrite(cmg_file,header.Layer{i}(2),'ubit8');
    fwrite(cmg_file,0,'ubit8');
    fwrite(cmg_file,0,'ubit8');
    fwrite(cmg_file,0,'ubit8');
    fwrite(cmg_file,0,'ubit8');
    fwrite(cmg_file,0,'ubit8');
    fwrite(cmg_file,0,'ubit8');
    fwrite(cmg_file,0,'ubit8');
    fwrite(cmg_file,0,'ubit8');
    fwrite(cmg_file,0,'ubit8');
    fwrite(cmg_file,0,'ubit8');
    fwrite(cmg_file,header.NumFeatures(i),'ubit8');
    fwrite(cmg_file,header.RGB_Order{i}(1),'ubit8');
    fwrite(cmg_file,header.RGB_Order{i}(2),'ubit8');
    fwrite(cmg_file,header.RGB_Order{i}(3),'ubit8');
    fwrite(cmg_file,'$','ubit8');% right the terminator

    for j=1:NbColorMap
        for k=1:size(I_RGB,1)
            %for l=1:size(I_RGB,2)
                fwrite(cmg_file,I_RGB(k,:,j),'ubit8');
            %end
        end
    end
 
    for j=1:NbBitMap
        for k=1:size(I_Bitmap,1)
            %for l=1:size(I_Bitmap,2)
                fwrite(cmg_file,I_Bitmap(k,:,j),'ubit8');
            %end
        end
    end

end


fclose('all');

end_time = toc(cmg_time);
fprintf(1,'\b\b\b\b%3.0f%%\n',100);
%fprintf(1,'\b\b\b\b%5.3f%',end_time);
%disp('s')
disp(['Done...',num2str(end_time),'s'])
disp('-------------------------------------------------')
end