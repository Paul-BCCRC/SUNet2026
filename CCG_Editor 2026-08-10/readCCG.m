function [ Images, Masks, header, ImageLabel, MaskLabel, ilist ] = readCCG( path, filename )
if( ~strcmp(path(end), '\') )
    path = [path,'\'];
end
ccg_time = tic;

disp('-------------------------------------------------')
disp([path, filename, '.ccg'])
%disp('Loading CCG...')

%Set path to CCG
CCG_path=[path,filename,'.ccg'];

%Open CCG for reading
ccg_file=fopen(CCG_path,'r');

%Read ALL bytes from CCG file
ccg_data=fread(ccg_file, inf, 'uint8');

%Number of bytes
BYTE_TOTAL = numel(ccg_data);
%Current byte
BYTE_CURRENT = 1;
%Current Cell
icount = 1;

%-----------------------------Read CCG Header------------------------------
%Get number of cell objects
N = bytes2type(ccg_data( (BYTE_CURRENT):(BYTE_CURRENT+7) ), 'uint64');
BYTE_CURRENT = BYTE_CURRENT + 8;

%Get number of image planes
P = bytes2type(ccg_data( (BYTE_CURRENT):(BYTE_CURRENT+1) ), 'uint16');
BYTE_CURRENT = BYTE_CURRENT + 2;

%Space to store image labels
ImageLabel = char(ones(P, 16, 'uint8')*32);

%Read plane labels
for p=1:P
    labelbytes = bytes2type(ccg_data( (BYTE_CURRENT):(BYTE_CURRENT+15) ), 'uint8');
    BYTE_CURRENT = BYTE_CURRENT + 16;
    ImageLabel(p, :) = char(labelbytes)';
end

%Read vorx for each object
for n=1:N
    h_vorx(n) = bytes2type(ccg_data( (BYTE_CURRENT):(BYTE_CURRENT+3) ), 'uint32');
    BYTE_CURRENT = BYTE_CURRENT + 4;
end

%Read vory for each object
for n=1:N
    h_vory(n) = bytes2type(ccg_data( (BYTE_CURRENT):(BYTE_CURRENT+3) ), 'uint32');
    BYTE_CURRENT = BYTE_CURRENT + 4;
end

%Read memory locations for each object
for n=1:N
    memoryloc(n) = bytes2type(ccg_data( (BYTE_CURRENT):(BYTE_CURRENT+7) ), 'uint64');
    BYTE_CURRENT = BYTE_CURRENT + 8;
end
%--------------------------------------------------------------------------

prc_dis=0;
prc_ =  floor(double(BYTE_CURRENT)*100/double(BYTE_TOTAL));
fprintf(1,'Parsing CCG bytes...%3d%%\n',prc_);
%While not at end of file
while(BYTE_TOTAL > BYTE_CURRENT)

    prc_ =  floor(double(BYTE_CURRENT)*100/double(BYTE_TOTAL));
    if(prc_ > prc_dis)
        fprintf(1,'\b\b\b\b%3.0f%%',prc_);
        prc_dis = prc_;
    end
%--------------------------Parse Header Data-------------------------------
    %char(cmg_data(BYTE_CURRENT))
    Mode(icount) = uint8(ccg_data(BYTE_CURRENT+1));
    NbColorMap(icount) = uint8(ccg_data(BYTE_CURRENT+2));
    Class(icount) = bytes2type(ccg_data( (BYTE_CURRENT+3):(BYTE_CURRENT+6) ), 'uint32');
    Screenx(icount) = bytes2type(ccg_data( (BYTE_CURRENT+7):(BYTE_CURRENT+10) ), 'uint32');
    Screeny(icount) = bytes2type(ccg_data( (BYTE_CURRENT+11):(BYTE_CURRENT+14) ), 'uint32');
    Stagex(icount) = bytes2type(ccg_data( (BYTE_CURRENT+15):(BYTE_CURRENT+22) ), 'uint64');
    Stagey(icount) = bytes2type(ccg_data( (BYTE_CURRENT+23):(BYTE_CURRENT+30) ), 'uint64');
    Stagez(icount) = bytes2type(ccg_data( (BYTE_CURRENT+31):(BYTE_CURRENT+38) ), 'uint64');
    Resolution(icount) = bytes2type(ccg_data( (BYTE_CURRENT+39):(BYTE_CURRENT+42) ), 'single');
    LowThreshold(icount) = bytes2type(ccg_data( (BYTE_CURRENT+43):(BYTE_CURRENT+44) ), 'uint16');
    MidThreshold(icount) = bytes2type(ccg_data( (BYTE_CURRENT+45):(BYTE_CURRENT+46) ), 'uint16');
    Group(icount) = uint8(ccg_data(BYTE_CURRENT+47));
    Width(icount) = bytes2type(ccg_data( (BYTE_CURRENT+48):(BYTE_CURRENT+51) ), 'uint32');
    Height(icount) = bytes2type(ccg_data( (BYTE_CURRENT+52):(BYTE_CURRENT+55) ), 'uint32');
    imagesize = Width(icount)*Height(icount);    
    Accession(icount) = bytes2type(ccg_data( (BYTE_CURRENT+56):(BYTE_CURRENT+59) ), 'uint32');
    Iod(icount) = bytes2type(ccg_data( (BYTE_CURRENT+60):(BYTE_CURRENT+63) ), 'single');
    Fluor(icount) = uint8(ccg_data(BYTE_CURRENT+64));
    Diagnosis(icount) = bytes2type(ccg_data( (BYTE_CURRENT+65):(BYTE_CURRENT+66) ), 'uint16');
    RedFaction(icount) = bytes2type(ccg_data( (BYTE_CURRENT+67):(BYTE_CURRENT+70) ), 'single');
    GreenFaction(icount) = bytes2type(ccg_data( (BYTE_CURRENT+71):(BYTE_CURRENT+74) ), 'single');
    BlueFaction(icount) = bytes2type(ccg_data( (BYTE_CURRENT+75):(BYTE_CURRENT+78) ), 'single');
    Index(icount) = bytes2type(ccg_data( (BYTE_CURRENT+79):(BYTE_CURRENT+82) ), 'uint32');
    Objective(icount) = bytes2type(ccg_data( (BYTE_CURRENT+83):(BYTE_CURRENT+86) ), 'uint32');
    Calibrated(icount) = uint8(ccg_data(BYTE_CURRENT+87));
    StackX_int(icount) = bytes2type(ccg_data( (BYTE_CURRENT+88):(BYTE_CURRENT+91) ), 'uint32');
    StackY_int(icount) = bytes2type(ccg_data( (BYTE_CURRENT+92):(BYTE_CURRENT+95) ), 'uint32');
    NbBitMap(icount) = uint8(ccg_data(BYTE_CURRENT+96));
    CassettePosition(icount) = uint8(ccg_data(BYTE_CURRENT+97));
    vorx(icount) = bytes2type(ccg_data( (BYTE_CURRENT+98):(BYTE_CURRENT+101) ), 'uint32');
    vory(icount) = bytes2type(ccg_data( (BYTE_CURRENT+102):(BYTE_CURRENT+105) ), 'uint32');
    BestFocusFrame(icount) = uint8(ccg_data(BYTE_CURRENT+106));
    BackgroundFloat(icount) = bytes2type(ccg_data( (BYTE_CURRENT+107):(BYTE_CURRENT+110) ), 'single');
    PrimaryColourChannel(icount) = uint8(ccg_data(BYTE_CURRENT+111));
    
    % for i=1:2
    %     Layer{icount}(i) = uint8(ccg_data(BYTE_CURRENT+111+i));
    % end
    Layer(icount) = bytes2type(ccg_data( (BYTE_CURRENT+112):(BYTE_CURRENT+113) ), 'uint16');
    
    for i=1:9
        Points{icount}(i) = uint8(ccg_data(BYTE_CURRENT+113+i));
    end
    
    NumFeature(icount) = uint8(ccg_data(BYTE_CURRENT+123));
    
    for i=1:3
        RGB_Order{icount}(i) = uint8(ccg_data(BYTE_CURRENT+123+i));
    end
    
    % metapadding = 0;
    % isMeta = true;
    % while(isMeta)
    %     if(ccg_data(BYTE_CURRENT+127+metapadding) ~= 36)
    %         metapadding = metapadding+1;
    %     else
    %         isMeta = false;
    %     end
    % end

    %BYTE_CURRENT = BYTE_CURRENT + 128 + metapadding;
    BYTE_CURRENT = BYTE_CURRENT + 128;
%--------------------------------------------------------------------------




%------------------------Parse Mask Data-----------------------------------
    %Get Number of images for this cell
    ImageCount{icount} = bytes2type(ccg_data(BYTE_CURRENT:BYTE_CURRENT+1), 'uint16');

    %Get Number of masks for this cell
    MaskCount{icount} = bytes2type(ccg_data(BYTE_CURRENT+2:BYTE_CURRENT+3), 'uint16');
    Masks{icount} = zeros(Height(icount), Width(icount), MaskCount{icount}, 'uint8');

    %Advance current byte
    BYTE_CURRENT = BYTE_CURRENT + 4;

    %Get Cell interaction list
    ilen = bytes2type(ccg_data( (BYTE_CURRENT):(BYTE_CURRENT) ), 'uint8');
    BYTE_CURRENT = BYTE_CURRENT + 1;
    ilist{icount} = zeros(1, ilen, 'single');
    for il=1:ilen
        ilist{icount}(il) = bytes2type(ccg_data( (BYTE_CURRENT):(BYTE_CURRENT+3) ), 'single');
        BYTE_CURRENT = BYTE_CURRENT + 4;
    end
    
    %Temp space to store mask labels
    readStr = char(ones(MaskCount{icount}, 16, 'uint8')*32);

    %For each mask
    for n=1:MaskCount{icount}
        
        %Read mask label
        labelbytes = bytes2type(ccg_data( (BYTE_CURRENT):(BYTE_CURRENT+15) ), 'uint8');
        readStr(n, :) = char(labelbytes)';

        %Read chaincount
        ChainCount = bytes2type(ccg_data( (BYTE_CURRENT+16):(BYTE_CURRENT+17) ), 'uint16');

        %Advance current byte
        BYTE_CURRENT = BYTE_CURRENT + 18;

        %For each chaincode
        for m=1:ChainCount

            %Read starting X
            StartX = bytes2type(ccg_data( (BYTE_CURRENT):(BYTE_CURRENT+3) ), 'uint32');

            %Read starting Y
            StartY = bytes2type(ccg_data( (BYTE_CURRENT+4):(BYTE_CURRENT+7) ), 'uint32');

            %Read chain length
            CLength = bytes2type(ccg_data( (BYTE_CURRENT+8):(BYTE_CURRENT+11) ), 'uint32'); 

            %If no chaincode
            if(CLength==0)

                %Set starting pixel
                Masks{icount}(StartY,StartX,n) = 1;

                %Update current byte
                BYTE_CURRENT = BYTE_CURRENT + 12;

            else

                %Read in chaincode bytes
                for i=1:CLength
                    CCode(i) = uint8(ccg_data(BYTE_CURRENT+11+i));
                end

                %Convert chaincode to mask
                Masks{icount}(:,:,n) = generateMaskfromCCode(Masks{icount}(:,:,n), CCode, StartX, StartY);

                %Update current byte
                BYTE_CURRENT = BYTE_CURRENT + 12 + CLength;

                %Delete CCode so it can be a different size next chain
                clear CCode StartX StartY CLength
            end

        end
    end
    MaskLabel{icount} = readStr;
    
%--------------------------------------------------------------------------

    icount=icount+1;
    
end
fprintf(1,'\b\b\b\b%3.0f%%\n',100);
disp('')
clear ccg_data

%--------------------Read Multi page TIFF----------------------------------
%Get TIFF information
info = imfinfo([path, filename, '.TIFF']);

%Get number of pages
numberOfPages = length(info);

prc_dis=0;
prc_ =  floor(double(0)*100/double(numberOfPages));
fprintf(1,'Reading TIFF...%3d%%\n',prc_);
%For each Page
for k = 1 : numberOfPages
    prc_ =  floor(double(k)*100/double(numberOfPages));
    if(prc_ > prc_dis)
        fprintf(1,'\b\b\b\b%3.0f%%',prc_);
        prc_dis = prc_;
    end
    %Read the kth image in this multipage tiff file.
    Frame{k} = imread([path, filename, '.TIFF'], k);
end
fprintf(1,'\b\b\b\b%3.0f%%\n',100);
%--------------------------------------------------------------------------

%------------------------Parse Image Data----------------------------------
%     frame = zeros(Height(icount), Width(icount), NbColorMap(icount), 'uint8');
%     for n=1:NbColorMap(icount)
%         frame(:, :, n) = reshape(  cmg_data(  uint64( BYTE_CURRENT + double(imagesize*(uint32(n)-1)) )  :  uint64(BYTE_CURRENT + double(imagesize*uint32(n)) - 1)    ), Width(icount), Height(icount))';
%     end
prc_dis=0;
prc_ =  0;
fprintf(1,'Parsing Image data...%3d%%\n',prc_);

for i=1:(icount-1)

    prc_ =  floor(double(i)*100/double(icount-1));
    if(prc_ > prc_dis)
        fprintf(1,'\b\b\b\b%3.0f%%',prc_);
        prc_dis = prc_;
    end

    for k = 1 : numberOfPages
        img(:,:,k) = Frame{k}( (Screeny(i)+1):(Screeny(i)+Height(i)), (Screenx(i)+1):(Screenx(i)+Width(i)), :);
    end
    Images{i} = img;
    clear img
end
%--------------------------------------------------------------------------

%Make the CMG header
MakeHeader();

%Close the CCG
fclose('all');
Images = Images';
Masks = Masks';
%ImageLabel = ImageLabel';
MaskLabel = MaskLabel';
ilist = ilist';


end_time = toc(ccg_time);
fprintf(1,'\b\b\b\b%3.0f%%\n',100);
%fprintf(1,'\b\b\b\b%5.3f%',end_time);
%disp('s')
disp(['Done...',num2str(end_time),'s'])
disp('-------------------------------------------------')

    %Convert bytes to a data type
    function [returntype] = bytes2type(bytes, type)
        bytes = uint8(bytes);
        returntype = typecast( fliplr(bytes), type);%%flip bytes to change endianness
    end
    
    %Makes a CMG header
    function MakeHeader()
        header.Mode=Mode';
        header.NbColorMap=NbColorMap';        
        header.Class=Class';
        header.Screenx=Screenx';
        header.Screeny=Screeny';
        header.Stagex=Stagex';
        header.Stagey=Stagey';
        header.Stagez=Stagez';
        header.Resolution=Resolution';
        header.LowThreshold=LowThreshold';
        header.MidThreshold=MidThreshold';
        header.Group=Group';
        header.Width=Width';
        header.Height=Height';
        header.Accession=Accession';
        header.Iod=Iod';
        header.Fluor=Fluor';
        header.Diagnosis=Diagnosis';
        header.RedFaction=RedFaction';
        header.GreenFaction=GreenFaction';
        header.BlueFaction=BlueFaction';
        header.Index=Index';
        header.Objective=Objective';
        header.Calibrated=Calibrated';
        header.StackX_int=StackX_int';
        header.StackY_int=StackY_int';
        header.NbBitMap=NbBitMap';
        header.CassettePosition=CassettePosition';
        header.vorx=vorx';
        header.vory=vory';
        header.BestFocusFrame=BestFocusFrame';
        header.BackgroundFloat=BackgroundFloat';
        header.PrimaryColourChannel=PrimaryColourChannel';
        header.Layer=Layer';
        header.Points=Points';
        header.NumFeatures=NumFeature';
        header.RGB_Order=RGB_Order';

%         header.MaskCount=MaskCount;
%         header.StartX=StartX;
%         header.StartY=StartY;
%         header.CLength=CLength;
%         header.CCode=CCode;
    end

    function [mask] = generateMaskfromCCode(mask, CCode, starty, startx)
        move_index = [1, 0,-1,-1,-1,0,1,1;  % X / horizontal move
                     -1,-1,-1, 0, 1,1,1,0]; % Y / vertical move


        Cx = double(startx);
        Cy = double(starty);
        mask(Cy, Cx) = 1;

        Clen = numel(CCode);
        for c=1:Clen
            Cx = Cx + move_index(1, CCode(c));
            Cy = Cy + move_index(2, CCode(c));
            mask(Cy, Cx) = 1;
        end

        mask = imfill(mask, "holes");
        

    end

end