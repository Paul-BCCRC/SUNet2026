function [bitmap, nuclei_mask, overlap_mask, isAbort, useMask] = getMultiFreeHand_snap2pixels(yy, xx, hh, ff)
iyy = uint64(yy);
ixx = uint64(xx);

oframe = get(hh, 'CData');
frame = get(hh, 'CData');
bitmap = -1;
nuclei_mask = -1;
overlap_mask = -1;
rlist=0;
r={};
masks=zeros(iyy,ixx);
isAbort = 0;
useMask = 0;
maskcount = 0;
isDone=0;
allowkeys=1;

xp=[];
yp=[];

t = text(40,2,'Select a Nuclei...(Escape=Abort/Enter=Done)','FontSize',17,'FontWeight','bold','Color',[1 0 0]);
tt = text(40,4,['MaskCount = ',num2str(maskcount)],'FontSize',17,'FontWeight','bold','Color',[1 0 0]);
set(ff,'WindowKeyPressFcn',@key_pressed)
set(ff,'WindowButtonDownFcn',@start_pencil)
setptr(ff,'datacursor');

while(~isDone && ishandle(ff))
    pause(0.1)
    %disp('.')
end

for m=1:maskcount
    lastx=get(r{maskcount-(m-1)},'xdata');
    lasty=get(r{maskcount-(m-1)},'ydata');
    xr=lastx( 1:(end-rlist(maskcount-(m-1))) );
    yr=lasty( 1:(end-rlist(maskcount-(m-1))) );
    
    set(r{maskcount-(m-1)},'xdata',xr,'ydata',yr);
end


    function key_pressed(src,eventdata)
        if(allowkeys)
        key = get(ff,'CurrentKey');
            switch key
                case 'return'
                    isDone=1;
                    useMask=1;
                    delete(t)
                    t = text(2,2,'Saving...','FontSize',17,'FontWeight','bold','Color',[1 0 0]);
                    setptr(ff,'watch');
                    %close(ff);
                case 'escape'
                    if(maskcount==0)
                        isDone=1;
                        isAbort=1;
                        %close(ff);
                    end
                case 'backspace'
                    if(maskcount>0)
                        masks(:,:,end) = [];
                        [ bitmap, nuclei_mask, overlap_mask ] = maskArray2BitMap( masks );
%                         frame(:,:,3) = nuclei_mask/3;
%                         frame(:,:,1) = overlap_mask/3;
                        frame(:,:,1) = oframe(:,:,1) + overlap_mask/8;
%                        frame(:,:,1) = oframe(:,:,1) + overlap_mask/4;
                         frame(:,:,2) = oframe(:,:,2) + nuclei_mask/8;
                         frame(:,:,3) = oframe(:,:,3) + nuclei_mask/8;
                        
                        %rlist{maskcount}=[];
                        lastx=get(r{maskcount},'xdata');
                        lasty=get(r{maskcount},'ydata');
                        xr=lastx( 1:(end-rlist(maskcount)) );
                        yr=lasty( 1:(end-rlist(maskcount)) );
                    
                        set(r{maskcount},'xdata',xr,'ydata',yr);
                        set(hh, 'CData', frame)
                        rlist(maskcount)=[];                       
                        maskcount=maskcount-1;
                        if(maskcount==0)
                            delete(t)
                            t = text(2,2,'Select a Nuclei...(Escape=Abort/Enter=Done)','FontSize',17,'FontWeight','bold','Color',[1 0 0]);
                        end
                        delete(tt)
                        tt = text(2,4,['MaskCount = ',num2str(maskcount)],'FontSize',17,'FontWeight','bold','Color',[1 0 0]);
                    end
                otherwise

            end
        end
    end

    function start_pencil(src,eventdata)
        allowkeys=0;
        delete(t)
        t = text(2,2,'Selecting Nuclei...','FontSize',17,'FontWeight','bold','Color',[1 0 0]);
        
        coords=get(gca,'currentpoint'); %since this is the axes callback, src=gca
        x=floor(coords(1,1,1))+0.5;
        y=floor(coords(1,2,1))+0.5;
        %disp(['X:',num2str(x),' Y:',num2str(y)])
        %r=line(x, y, 'color', [0 .5 1], 'LineWidth', 2, 'hittest', 'off'); %turning     hittset off allows you to draw new lines that start on top of an existing line.
        
        r{maskcount+1}=line(x, y, 'color', [1 0 1], 'LineWidth',1);

        set(ff,'WindowButtonMotionFcn',{@continue_pencil})
        set(ff,'WindowButtonUpFcn',{@done_pencil})
    end

    function continue_pencil(src,eventdata)
        %Note: src is now the figure handle, not the axes, so we need to use gca.
        coords=get(gca,'currentpoint'); %this updates every time i move the mouse
        x=floor(coords(1,1,1))+0.5;
        y=floor(coords(1,2,1))+0.5;
        %disp(['X:',num2str(x),' Y:',num2str(y)])
        %get the line's existing coordinates and append the new ones.
        lastx=get(r{maskcount+1},'xdata');
        lasty=get(r{maskcount+1},'ydata');
        xp=[lastx x];
        yp=[lasty y];
        set(r{maskcount+1},'xdata',xp,'ydata',yp);
    end

    function [ ] = done_pencil(src,evendata)
        %all this funciton does is turn the motion function off
        set(ff,'WindowButtonMotionFcn','')
        set(ff,'WindowButtonupFcn','')   
        
        maskcount=maskcount+1;
        
        %%%Close selection
        lastx=get(r{maskcount},'xdata');
        lasty=get(r{maskcount},'ydata');
        xp=[lastx lastx(1)];
        yp=[lasty lasty(1)];
        set(r{maskcount},'xdata',xp,'ydata',yp);

        %selected_mask = polylist2Mask(xp, yp, yy, xx, 1);
        %selected_mask = poly2mask(xp, yp, iyy, ixx);
        selected_mask = images.internal.builtins.poly2mask(xp,yp, uint32(yy), uint32(xx));
        masks(:,:,maskcount) = selected_mask;
        [ bitmap, nuclei_mask, overlap_mask ] = maskArray2BitMap( masks );
%         frame(:,:,3) = nuclei_mask/3;
%         frame(:,:,1) = overlap_mask/3;
        frame(:,:,1) = oframe(:,:,1) + overlap_mask/8;
%        frame(:,:,1) = oframe(:,:,1) + overlap_mask/4;
        frame(:,:,2) = oframe(:,:,2) + nuclei_mask/8;
        frame(:,:,3) = oframe(:,:,3) + nuclei_mask/8;
                        
        set(hh, 'CData', frame)  
        
        rlist(maskcount) = numel(xp);
        delete(t);
        delete(tt);
        t = text(2,2,'Select New Nuclei...(Backspace=Undo/Enter=Done)','FontSize',17,'FontWeight','bold','Color',[1 0 0]);
        tt = text(2,4,['MaskCount = ',num2str(maskcount)],'FontSize',17,'FontWeight','bold','Color',[1 0 0]);
        allowkeys=1;
        
        
        
    end

end