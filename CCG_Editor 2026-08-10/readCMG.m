function [ Images, Masks, header ] = readCMG( path, filename )
if( ~strcmp(path(end), '\') )
    path = [path,'\'];
end
cmg_time = tic;

disp('-------------------------------------------------')
disp([path, filename, '.cmg'])

s = dir([path, filename, '.cmg']);
BYTE_TOTAL = s.bytes;
CMG_path=[path,filename,'.cmg'];
cmg_file=fopen(CMG_path,'r');


BYTE_CURRENT = 1;
BYTE_FILEPOSITION = 1;
icount = 1;

prc_dis=0;
prc_ =  floor(double(BYTE_FILEPOSITION)*100/double(BYTE_TOTAL));
fprintf(1,'Reading CMG...%3d%%\n',prc_);
while(BYTE_TOTAL > BYTE_FILEPOSITION)

    prc_ =  floor(double(BYTE_FILEPOSITION)*100/double(BYTE_TOTAL));
    if(prc_ > prc_dis)
        fprintf(1,'\b\b\b\b%3.0f%%',prc_);
        prc_dis = prc_;
    end

    header_data=fread(cmg_file, 128, 'uint8');

%--------------------------Parse Header Data-------------------------------
    %char(cmg_data(BYTE_CURRENT))
    Mode(icount) = uint8(header_data(BYTE_CURRENT+1));
    NbColorMap(icount) = uint8(header_data(BYTE_CURRENT+2));
    Class(icount) = bytes2type(header_data( (BYTE_CURRENT+3):(BYTE_CURRENT+6) ), 'uint32');
    Screenx(icount) = bytes2type(header_data( (BYTE_CURRENT+7):(BYTE_CURRENT+10) ), 'uint32');
    Screeny(icount) = bytes2type(header_data( (BYTE_CURRENT+11):(BYTE_CURRENT+14) ), 'uint32');
    Stagex(icount) = bytes2type(header_data( (BYTE_CURRENT+15):(BYTE_CURRENT+22) ), 'uint64');
    Stagey(icount) = bytes2type(header_data( (BYTE_CURRENT+23):(BYTE_CURRENT+30) ), 'uint64');
    Stagez(icount) = bytes2type(header_data( (BYTE_CURRENT+31):(BYTE_CURRENT+38) ), 'uint64');
    Resolution(icount) = bytes2type(header_data( (BYTE_CURRENT+39):(BYTE_CURRENT+42) ), 'single');
    LowThreshold(icount) = bytes2type(header_data( (BYTE_CURRENT+43):(BYTE_CURRENT+44) ), 'uint16');
    MidThreshold(icount) = bytes2type(header_data( (BYTE_CURRENT+45):(BYTE_CURRENT+46) ), 'uint16');
    Group(icount) = uint8(header_data(BYTE_CURRENT+47));
    Width(icount) = bytes2type(header_data( (BYTE_CURRENT+48):(BYTE_CURRENT+51) ), 'uint32');
    Height(icount) = bytes2type(header_data( (BYTE_CURRENT+52):(BYTE_CURRENT+55) ), 'uint32');
    imagesize = Width(icount)*Height(icount);    
    Accession(icount) = bytes2type(header_data( (BYTE_CURRENT+56):(BYTE_CURRENT+59) ), 'uint32');
    Iod(icount) = bytes2type(header_data( (BYTE_CURRENT+60):(BYTE_CURRENT+63) ), 'single');
    Fluor(icount) = uint8(header_data(BYTE_CURRENT+64));
    Diagnosis(icount) = bytes2type(header_data( (BYTE_CURRENT+65):(BYTE_CURRENT+66) ), 'uint16');
    RedFaction(icount) = bytes2type(header_data( (BYTE_CURRENT+67):(BYTE_CURRENT+70) ), 'single');
    GreenFaction(icount) = bytes2type(header_data( (BYTE_CURRENT+71):(BYTE_CURRENT+74) ), 'single');
    BlueFaction(icount) = bytes2type(header_data( (BYTE_CURRENT+75):(BYTE_CURRENT+78) ), 'single');
    Index(icount) = bytes2type(header_data( (BYTE_CURRENT+79):(BYTE_CURRENT+82) ), 'uint32');
    Objective(icount) = bytes2type(header_data( (BYTE_CURRENT+83):(BYTE_CURRENT+86) ), 'uint32');
    Calibrated(icount) = uint8(header_data(BYTE_CURRENT+87));
    StackX_int(icount) = bytes2type(header_data( (BYTE_CURRENT+88):(BYTE_CURRENT+91) ), 'uint32');
    StackY_int(icount) = bytes2type(header_data( (BYTE_CURRENT+92):(BYTE_CURRENT+95) ), 'uint32');
    NbBitMap(icount) = uint8(header_data(BYTE_CURRENT+96));
    CassettePosition(icount) = uint8(header_data(BYTE_CURRENT+97));
    vorx(icount) = bytes2type(header_data( (BYTE_CURRENT+98):(BYTE_CURRENT+101) ), 'uint32');
    vory(icount) = bytes2type(header_data( (BYTE_CURRENT+102):(BYTE_CURRENT+105) ), 'uint32');
    BestFocusFrame(icount) = uint8(header_data(BYTE_CURRENT+106));
    BackgroundFloat(icount) = bytes2type(header_data( (BYTE_CURRENT+107):(BYTE_CURRENT+110) ), 'single');
    PrimaryColourChannel(icount) = uint8(header_data(BYTE_CURRENT+111));

    % for i=1:2
    %     Layer{icount}(i) = uint8(header_data(BYTE_CURRENT+111+i));
    % end
    Layer(icount) = bytes2type(header_data( (BYTE_CURRENT+112):(BYTE_CURRENT+113) ), 'uint16');

    for i=1:9
        Points{icount}(i) = uint8(header_data(BYTE_CURRENT+113+i));
    end

    NumFeature(icount) = uint8(header_data(BYTE_CURRENT+123));

    for i=1:3
        RGB_Order{icount}(i) = uint8(header_data(BYTE_CURRENT+123+i));
    end

    % metapadding = 0;
    % isMeta = true;
    % while(isMeta)
    %     if(cmg_data(BYTE_CURRENT+127+metapadding) ~= 36)
    %         metapadding = metapadding+1;
    %     else
    %         isMeta = false;
    %     end
    % end
    % BYTE_CURRENT = BYTE_CURRENT + 128 + metapadding;
    BYTE_FILEPOSITION = BYTE_FILEPOSITION + 128;
%--------------------------------------------------------------------------

%------------------------Parse Image Data----------------------------------
    image_data=fread(cmg_file, double(imagesize) * double(NbColorMap(icount)), 'uint8');
    frame = zeros(Height(icount), Width(icount), NbColorMap(icount), 'uint8');
    for n=1:NbColorMap(icount)
        frame(:, :, n) = reshape(  image_data(  uint64( BYTE_CURRENT + double(imagesize*(uint32(n)-1)) )  :  uint64(BYTE_CURRENT + double(imagesize*uint32(n)) - 1)    ), Width(icount), Height(icount))';
    end
    Images{icount} = frame;
    BYTE_FILEPOSITION = BYTE_FILEPOSITION + double(imagesize)*double(NbColorMap(icount));
%--------------------------------------------------------------------------

%------------------------Parse Mask Data-----------------------------------
    mask_data=fread(cmg_file, double(imagesize) * double(NbBitMap(icount)), 'uint8');
    mask = zeros(Height(icount), Width(icount), NbBitMap(icount), 'uint8');
    for n=1:NbBitMap(icount)
        mask(:, :, n) = reshape(mask_data( uint64( BYTE_CURRENT + double(imagesize*(uint32(n)-1)) )  :  uint64(BYTE_CURRENT + double(imagesize*uint32(n)) - 1)   ), Width(icount), Height(icount))';
    end
    Masks{icount} = mask;
    BYTE_FILEPOSITION = BYTE_FILEPOSITION + double(imagesize)*double(NbBitMap(icount));
%--------------------------------------------------------------------------
    icount=icount+1;
    clear header_data image_data mask_data
end

MakeHeader();

fclose('all');
Images = Images';
Masks = Masks';

end_time = toc(cmg_time);
fprintf(1,'\b\b\b\b%3.0f%%\n',100);
%fprintf(1,'\b\b\b\b%5.3f%',end_time);
%disp('s')
disp(['Done...',num2str(end_time),'s'])
disp('-------------------------------------------------')


    function [returntype] = bytes2type(bytes, type)
        bytes = uint8(bytes);
        returntype = typecast( fliplr(bytes), type);%%flip bytes to change endianness
    end

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
    end

end