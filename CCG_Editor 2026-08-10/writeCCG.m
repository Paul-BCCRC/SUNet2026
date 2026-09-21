function [  ] = writeCCG( header_CMG, ImageLabel, MaskLabel, Images, Masks, path, filename )
if( ~strcmp(path(end), '\') )
    path = [path,'\'];
end
ccg_time = tic;

disp('-------------------------------------------------')
disp([path, filename, '.ccg'])
%disp('Writing CCG...')

[h_minX] = min(header_CMG.Screenx);
[h_minY] = min(header_CMG.Screeny);

header_CMG.Screenx = header_CMG.Screenx - (h_minX-1);
header_CMG.Screeny = header_CMG.Screeny - (h_minY-1);

header_CMG.vorx = header_CMG.vorx - (h_minX-1);
header_CMG.vory = header_CMG.vory - (h_minY-1);

%------------------------Save Image data to TIFF---------------------------
%Construct full frame image from CMG objects
disp('Constructing TIFF...')
[DATA_frame, numChannels, header_CMG]= createFullData(Images, header_CMG);

%Save frame as multi-page TIFF
disp('Writing TIFF...')
for k=1:numChannels
    if(k==1)
        imwrite( DATA_frame(:,:,k), [path, filename,'.TIFF'], 'Compression', 'lzw');
    else
        imwrite( DATA_frame(:,:,k), [path, filename,'.TIFF'], 'Compression', 'lzw', 'WriteMode', 'append');
    end
end
%--------------------------------------------------------------------------

%Set CCG location
CCG_path=[path,'\',filename,'.ccg'];

%Set CCG file for writing
ccg_file=fopen(CCG_path,'w');

%Get total number of Cells
N=numel(Masks);

CURRENT_BYTE = uint64(0);

%
disp('Finding Interacting Cells...')
[I_Obj] = findOverlap(Masks, header_CMG);

%-------------------------------CCG header---------------------------------
%Write number of cell objects in CCG
fwrite(ccg_file, uint64(N), 'ubit64');
CURRENT_BYTE = CURRENT_BYTE + uint64(8);

%Write number of image planes
fwrite(ccg_file, uint16(numChannels), 'ubit16');
CURRENT_BYTE = CURRENT_BYTE + uint64(2);

%Write each image label
for n=1:numChannels
    try
        labelbytes = ones(1,16, 'uint8')*32;
        numchar = numel(ImageLabel{n});
        if(numchar <= 16)
            labelbytes(1, 1:numchar) = ImageLabel{n};
        else
            labelbytes(1, 1:16) = ImageLabel{n}(1:16);
        end

        %labelbytes = uint8( ImageLabel{n} );
    catch
        %String of spaces 32=space in ascii
        labelbytes = ones(1,16, 'uint8')*32;
    end
    for s=1:16
        fwrite(ccg_file, labelbytes(s), 'ubit8');
    end
    CURRENT_BYTE = CURRENT_BYTE + uint64(16);
    clear labelbytes
end

%Write vorx for each object
for n=1:N
    fwrite(ccg_file, header_CMG.vorx(n), 'ubit32');
    CURRENT_BYTE = CURRENT_BYTE + uint64(4);
end

%Write vory for each object
for n=1:N
    fwrite(ccg_file, header_CMG.vory(n), 'ubit32');
    CURRENT_BYTE = CURRENT_BYTE + uint64(4);
end

MEMSTART = CURRENT_BYTE;
%Allocate space for memory pointers
memlist = zeros(N, 1, 'uint64');
for n=1:N
    fwrite(ccg_file, 0,'ubit64');
    CURRENT_BYTE = CURRENT_BYTE + uint64(8);
end

%--------------------------------------------------------------------------

%For each Cell
prc_dis=0;
prc_ =  floor(double(1)*100/double(N));
fprintf(1,'Writing Chain Code...%3d%%\n',prc_);

for n=1:N
    %disp(['Cell ',num2str(n),'/',num2str(N), ' : ',num2str( uint8(double(n)*100/double(N)) ),'%'])

    prc_ =  floor(double(n)*100/double(N));
    if(prc_ > prc_dis)
        fprintf(1,'\b\b\b\b%3.0f%%',prc_);
        prc_dis = prc_;
    end

    memlist(n) = CURRENT_BYTE;
    %disp([num2str(n), ' of ', num2str(N)])

    %Get mask shape
    [Y, X, Z] = size(Masks{n});

    %Copy/Paste CMG header to the CCG header
    writecmgheader(header_CMG, n, ccg_file, numChannels, X, Y, Z)
    CURRENT_BYTE = CURRENT_BYTE + uint64(128);

    %Write image count in CCG header
    fwrite(ccg_file, uint16(numChannels), 'ubit16');
    CURRENT_BYTE = CURRENT_BYTE + uint64(2);

    %Write mask count in CCG header
    fwrite(ccg_file, uint16(Z), 'ubit16');
    CURRENT_BYTE = CURRENT_BYTE + uint64(2);

    %Write overlap mask count and values
    I_Obj_num = numel(I_Obj{n});
    fwrite(ccg_file, uint8(I_Obj_num), 'ubit8');
    for i=1:I_Obj_num
        fwrite(ccg_file, I_Obj{n}(i), 'single');
    end
    CURRENT_BYTE = CURRENT_BYTE + uint64(4*I_Obj_num + 1);

    %For each mask
    for z=1:Z
        %Convert to logical
        mask = imbinarize(Masks{n}(:,:,z));

        %Write mask label
        try
            labelbytes = ones(1,16, 'uint8')*32;
            numchar = numel(MaskLabel{z});
            if(numchar <= 16)
                labelbytes(1, 1:numchar) = MaskLabel{z};
            else
                labelbytes(1, 1:16) = MaskLabel{z}(1:16);
            end

            %labelbytes = uint8( MaskLabel{n}(z, :) );
        catch
            %String of spaces 32=space in ascii
            labelbytes = ones(1,16, 'uint8')*32;
        end
        for s=1:16
            fwrite(ccg_file, labelbytes(s), 'ubit8');
        end
        clear labelbytes
        %fwrite(ccg_file, uint8(0), 'ubit8');
        %fwrite(ccg_file, uint8(0), 'ubit8');
        CURRENT_BYTE = CURRENT_BYTE + uint64(16);

        %Count pixels
        masksum = sum(mask(:));

        if(masksum==0)
            %Mask is blank

            %Write chain count = 0
            fwrite(ccg_file, uint16(0), 'ubit16');
            CURRENT_BYTE = CURRENT_BYTE + uint64(2);

        else
            %Number each seperate object in mask and get a count
            [labelmask, M] = bwlabel(mask, 4);

            %Write chain count
            fwrite(ccg_file, uint16(M), 'ubit16');
            CURRENT_BYTE = CURRENT_BYTE + uint64(2);

            %For each seperate mask object
            for m=1:M
                %Get sub mask
                submask = zeros(Y, X, 'logical');
                submask(find(labelmask==m))=1;

                %Count pixels of submask
                submasksum = sum(submask(:));

                %Check for unchainable single pixel mask
                if(submasksum==1)

                    %Find the pixel location
                    [I,J] = find(submask, 1);

                    %Write StartX
                    fwrite(ccg_file, uint32(J), 'ubit32');
                    CURRENT_BYTE = CURRENT_BYTE + uint64(4);
                    %Write StartY
                    fwrite(ccg_file, uint32(I), 'ubit32');
                    CURRENT_BYTE = CURRENT_BYTE + uint64(4);
                    %Write Chain Length = 0
                    fwrite(ccg_file, uint32(0), 'ubit32');
                    CURRENT_BYTE = CURRENT_BYTE + uint64(4);

                else

                    %Get submask chain code
                    [~, X0, Code, ~, ~] = Freeman_chain_code(submask, 0);

                    %remove chaincode from cell array
                    Code = Code{1};

                    %Get chain code length
                    C = uint32(numel(Code));

                    %Cast Start points to uint32
                    X0 = uint32(X0);

                    %Write StartX
                    fwrite(ccg_file, X0(1), 'ubit32');
                    CURRENT_BYTE = CURRENT_BYTE + uint64(4);
                    %Write StartY
                    fwrite(ccg_file, X0(2), 'ubit32');
                    CURRENT_BYTE = CURRENT_BYTE + uint64(4);
                    %Write Chain Length
                    fwrite(ccg_file, C, 'ubit32');
                    CURRENT_BYTE = CURRENT_BYTE + uint64(4);
                    %Write Chain Data
                    for c=1:C
                        fwrite(ccg_file, uint8(Code(c)), 'ubit8');
                        CURRENT_BYTE = CURRENT_BYTE + uint64(1);
                    end

                end

            end

        end

    end

end

%Set write location to begining of Memory locations.
fseek(ccg_file, MEMSTART, -1);

%Overwrite dummy(zeros) memory locations
for n=1:N
    fwrite(ccg_file, memlist(n),'ubit64');
end


fclose('all');

end_time = toc(ccg_time);
fprintf(1,'\b\b\b\b%3.0f%%\n',100);
%fprintf(1,'\b\b\b\b%5.3f%',end_time);
%disp('s')
disp(['Done...',num2str(end_time),'s'])
disp('-------------------------------------------------')

    function [] = writecmgheader(header, i, cmg_file, N, X, Y, Z)
        fwrite(cmg_file,'c','ubit8');
        fwrite(cmg_file,header.Mode(i),'ubit8');
        fwrite(cmg_file,uint8(N),'ubit8');
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
        fwrite(cmg_file,uint32(X),'ubit32');%width of the image
        fwrite(cmg_file,uint32(Y),'ubit32');%height of the image
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
        fwrite(cmg_file,uint8(Z),'ubit8');
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
    end

    function [DATA_frame, maxZ, header] = createFullData(Images, header)
        NN = numel(header.Mode);


        [minX] = min(header.Screenx);
        [minY] = min(header.Screeny);

        [maxX, I_X] = max(header.Screenx);
        [maxY, I_Y] = max(header.Screeny);
        Nimg = numel(Images);
        if(Nimg > 100)
            Nbackimg = 100;
        else
            Nbackimg = Nimg;
        end
        maxZ = 1;
        borderpixels = 128;
        for nn=1:Nimg
            [~, ~, Zimg] = size(Images{nn});
            if(Zimg > maxZ)
                maxZ = Zimg;
            end
            if(Nbackimg < nn)
                borderpixels = cat(1, borderpixels, Images{nn}(1,:,1)' );
                borderpixels = cat(1, borderpixels, Images{nn}(:,1,1) );
                borderpixels = cat(1, borderpixels, Images{nn}(end,:,1)' );
                borderpixels = cat(1, borderpixels, Images{nn}(:,end,1) );
            end
        end

        AvgBP = mean(borderpixels);

        XX = maxX - minX + header.Width(I_X);
        YY = maxY - minY + header.Height(I_Y);

        if(AvgBP>=128)
            DATA_frame = ones(YY, XX, maxZ, 'uint8')*255;
        else
            DATA_frame = zeros(YY, XX, maxZ, 'uint8');
        end
        for nn=1:NN
            posx = header.Screenx(nn) - minX;
            posy = header.Screeny(nn) - minY;
            [YY, XX, ZZ] = size(Images{nn});

            for zz=1:ZZ
                DATA_frame(posy+1:posy+YY, posx+1:posx+XX, zz) = Images{nn}(:,:,zz);
            end

        end

        header.Screenx = header.Screenx - minX;
        header.Screeny = header.Screeny - minY;

    end

    function [total_list] = findOverlap(Masks, header)
        NN = numel(Masks);
        for nn=1:NN
            SX = double( header.Screenx(nn) );
            SY = double( header.Screeny(nn) );
            XX = double( header.Width(nn) );
            YY = double( header.Height(nn) );

            top = floor(SY);
            bottom = floor(SY) + ceil(YY);
            left = floor(SX);
            right = floor(SX) + ceil(XX);

            xp = header.vorx;
            yp = header.vory;

            xrange = find(xp >= left & xp <= right);
            yrange = find(yp >= top & yp <= bottom);

            L = ismember(xrange, yrange);
            L = L .* xrange;
            L(L==0) = [];
            List = L;
            LL = numel(L);

            Xmax=0;
            Ymax=0;
            for l=1:LL
                if(Xmax < header.Width(L(l)))
                    Xmax = header.Width(L(l));
                end
                if(Ymax < header.Height(L(l)))
                    Ymax = header.Height(L(l));
                end
            end

            for l=1:LL
                ll = LL - (l-1);
                nSX = double( header.Screenx(L(ll)) );
                nSY = double( header.Screeny(L(ll)) );
                nXX = double( header.Width(L(ll)) );
                nYY = double( header.Height(L(ll)) );

                OX = Xmax + (nSX - SX);
                OY = Ymax + (nSY - SY);

                Omask = zeros(YY+2*Ymax, XX+2*Xmax, 'uint8');
                Omask( (1+Ymax):(YY+Ymax) , (1+Xmax):(XX+Xmax) ) = Masks{nn}(:,:,1);
                Omask( (1+OY):(nYY+OY), (1+OX):(nXX+OX) ) = Omask( (1+OY):(nYY+OY), (1+OX):(nXX+OX) ) + Masks{L(ll)}(:,:,1);

                if(max(Omask,[],"all") <= 1 || List(ll)==nn)
                    List(ll) = [];
                end

            end

            total_list{nn} = single(List);
        end

    end

end