%function [handles, currentCell, sMode, isAbort] = selectCell(app, event, hh, ff, handles)
function [handles, currentCell, sMode, isAbort] = selectCell(app, event)
[hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event);
hh = handles.MapFig;
ff = handles.MapFigHandle;
%oframe = get(hh, 'CData');
frame = get(hh, 'CData');
[yy, xx, ~] = size(frame);
%masks=zeros(yy,xx,'uint8');

currentCell = -1;
sMode = -1;
isAbort=0;
rlist=0;
%r={};

isDone=0;
allowkeys=1;

xp=[];
yp=[];

%imgHandle = findobj(ancestor(hh, 'figure'), 'Type', 'image');
% set(handles.MapFigViewer, 'ButtonDownFcn', @start_pencil);

set(ff,'WindowKeyPressFcn',@key_pressed)
set(ff,'WindowButtonDownFcn',@start_pencil)
setptr(ff,'datacursor');
%count=0
while(~isDone)
    pause(0.1)
    %[hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event);
    if(findobj('type','figure','name','CMG Map')==1)
        
    else
        isDone = true;
        isAbort = true;
    end

%    disp(count)
%    count=count+1;
end


    function key_pressed(src,eventdata)
        if(allowkeys)
        key = get(ff,'CurrentKey');
        %disp(key)
            switch key
                case 'insert'
                    isDone=1;
                    sMode=0;                
                case 'return'
                    isDone=1;
                    sMode=1;
                    %handles.full_mask = frame;
%                    delete(t)
%                    t = text(2,2,'Saving...','FontSize',17,'FontWeight','bold','Color',[1 0 0]);
%                    setptr(ff,'watch');
                    %close(ff);
%                 case 'escape'
%                         isDone=1;
%                         isAbort=1;
%                         close(ff);
                case 'home'
                    isDone=1;
                    sMode=2;
                case 'equal'
                    isDone=1;
                    sMode=3;
                case 'hyphen'
                    isDone=1;
                    sMode=4;
                case 'backspace'
                    
                case 'space'

                case 'delete'
                    isDone=1;
                    sMode=5;
                otherwise

            end
        end
    end


    function start_pencil(src,eventdata)
        allowkeys=0;
        

        coords=get(gca,'currentpoint'); %since this is the axes callback, src=gca
        x=coords(1,1,1);
        y=coords(1,2,1);

        [ I ] = getNearestCell(x, y);
        if(I~=currentCell)
            updateMap( I );
            currentCell=I;
            if(ishandle(handles.FeatureViewFigHandle))
                handles.FeatureViewFigHandle.Children(1).Data(:,2) = num2cell(handles.featureData(:, I));
            end
            if(ishandle(handles.FeatureFigHandle))
                handles.FeatureFigHandle.Children(6).Data(:,2) = num2cell(handles.featureData(:, I));
            end
        end

        set(ff,'WindowButtonMotionFcn',{@continue_pencil})
        set(ff,'WindowButtonUpFcn',{@done_pencil})
    end

    function continue_pencil(src,eventdata)
        %Note: src is now the figure handle, not the axes, so we need to use gca.
        coords=get(gca,'currentpoint'); %this updates every time i move the mouse
        x=coords(1,1,1);
        y=coords(1,2,1);
        %disp(['X:',num2str(x),'Y:',num2str(y)])

        [ I ] = getNearestCell(x, y);
        if(I~=currentCell)
            updateMap( I );
            currentCell=I;
            if(ishandle(handles.FeatureViewFigHandle))
                handles.FeatureViewFigHandle.Children(1).Data(:,2) = num2cell(handles.featureData(:, I));
            end
            if(ishandle(handles.FeatureFigHandle))
                handles.FeatureFigHandle.Children(6).Data(:,2) = num2cell(handles.featureData(:, I));
            end
        end

    end

    function [ ] = done_pencil(src,evendata)
        %all this funciton does is turn the motion function off
        set(ff,'WindowButtonMotionFcn','')
        set(ff,'WindowButtonupFcn','')   
        
        coords=get(gca,'currentpoint');
        x=coords(1,1,1);
        y=coords(1,2,1);

        [ I ] = getNearestCell(x, y);
        if(I~=currentCell)
            updateMap( I );
            currentCell=I;
            if(ishandle(handles.FeatureViewFigHandle))
                handles.FeatureViewFigHandle.Children(1).Data(:,2) = num2cell(handles.featureData(:, I));
            end
            if(ishandle(handles.FeatureFigHandle))
                handles.FeatureFigHandle.Children(6).Data(:,2) = num2cell(handles.featureData(:, I));
            end
        end

        allowkeys=1;
        
    end

    function [ I ] = getNearestCell(xx, yy)
        dx = double(handles.vorx) - xx;
        dy = double(handles.vory) - yy;
        d = sqrt(abs(dx).^2 + abs(dy).^2);
        [~, I] = min(d);
    end

    function [] = updateMap( I )
        %[hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event);
        isMask = get(app.showmask_checkbox, 'Value');

        if(isMask)
            handles = guidata(hObject);
            frame = handles.full_mask;
        else
            handles = guidata(hObject);
            frame = handles.full_raw;
        end
        %SE = strel('square',3);
            posx = handles.Screenx(I);
            posy = handles.Screeny(I);
            if(numel(handles.ilist)>0)
                List = handles.ilist{I};
                L = numel(List);
            end

        for n=1:handles.header.NbBitMap(I)
            if(handles.maskNum(n))

                if(numel(handles.ilist)>0)
                    for l=1:L
                        II = List(l);
                        Iposx = handles.Screenx(II);
                        Iposy = handles.Screeny(II);
                        Imask = handles.dMasks{II}(:,:,n)*255;
                        [Y, X, Z] = size(Imask);

                        frame(Iposy+1:Iposy+Y, Iposx+1:Iposx+X,:) = frame(Iposy+1:Iposy+Y, Iposx+1:Iposx+X,:) - Imask;
            
                        %frame(Iposy+1:Iposy+Y, Iposx+1:Iposx+X,1) = frame(Iposy+1:Iposy+Y, Iposx+1:Iposx+X,1) + Imask;
                        %frame(Iposy+1:Iposy+Y, Iposx+1:Iposx+X,2) = frame(Iposy+1:Iposy+Y, Iposx+1:Iposx+X,2) + Imask;
                        %frame(Iposy+1:Iposy+Y, Iposx+1:Iposx+X,3) = frame(Iposy+1:Iposy+Y, Iposx+1:Iposx+X,3) - (Imask-105);
                        
                        frame(Iposy+1:Iposy+Y, Iposx+1:Iposx+X,:) = frame(Iposy+1:Iposy+Y, Iposx+1:Iposx+X,:) + (Imask-128);%205

                    end                    
                end

                mask = handles.dMasks{I}(:,:,n)*255;
                [Y, X, Z] = size(mask);
                %masks(posy+1:posy+Y, posx+1:posx+X) = masks(posy+1:posy+Y, posx+1:posx+X) + mask*255;

                frame(posy+1:posy+Y, posx+1:posx+X,:) = frame(posy+1:posy+Y, posx+1:posx+X,:) - mask;
    
                frame(posy+1:posy+Y, posx+1:posx+X,1) = frame(posy+1:posy+Y, posx+1:posx+X,1) + (mask-105);
                frame(posy+1:posy+Y, posx+1:posx+X,2) = frame(posy+1:posy+Y, posx+1:posx+X,2) + mask;
                %frame(posy+1:posy+Y, posx+1:posx+X,3) = frame(posy+1:posy+Y, posx+1:posx+X,3);
            end
        end
        
        set(hh, 'CData', frame) 



        % if(handles.header.NbBitMap(I) >= numMask)
        %     SE = strel('square',3);
        %     posx = handles.header.Screenx(I)-handles.X_Min;
        %     posy = handles.header.Screeny(I)-handles.Y_Min;
        %     %masks=zeros(yy,xx,'uint8');
        % 
        %     %mask = imdilate(handles.Masks{I}, SE) - handles.Masks{I};
        %     %mask = mask(:,:,numMask)*255;
        %     mask = handles.dMasks{I}(:,:,numMask)*255;
        % 
        % 
        %     [Y, X, Z] = size(mask);
        %     %masks(posy+1:posy+Y, posx+1:posx+X) = masks(posy+1:posy+Y, posx+1:posx+X) + mask*255;
        % 
        %     if(isMask)
        %         handles = guidata(hObject);
        %         frame = handles.full_mask;
        %     else
        %         handles = guidata(hObject);
        %         frame = handles.full_raw;
        %     end
        % 
        %     frame(posy+1:posy+Y, posx+1:posx+X,1) = frame(posy+1:posy+Y, posx+1:posx+X,1) - mask;
        %     frame(posy+1:posy+Y, posx+1:posx+X,2) = frame(posy+1:posy+Y, posx+1:posx+X,2) - mask;
        %     frame(posy+1:posy+Y, posx+1:posx+X,3) = frame(posy+1:posy+Y, posx+1:posx+X,3) - mask;
        % 
        %     frame(posy+1:posy+Y, posx+1:posx+X,1) = frame(posy+1:posy+Y, posx+1:posx+X,1) + (mask-105);
        %     frame(posy+1:posy+Y, posx+1:posx+X,2) = frame(posy+1:posy+Y, posx+1:posx+X,2) + mask;
        %     frame(posy+1:posy+Y, posx+1:posx+X,3) = frame(posy+1:posy+Y, posx+1:posx+X,3) - (mask-105);
        % 
        % 
        %     set(hh, 'CData', frame) 
        % 
        % end
    end

end