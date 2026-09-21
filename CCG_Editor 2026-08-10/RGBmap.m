function varargout = RGBmap(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RGBmap_OpeningFcn, ...
                   'gui_OutputFcn',  @RGBmap_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
end


function RGBmap_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
res = get(0,'screensize');
resX = res(3);
resY = res(4);
clear res

handles.isEdit = false;
handles.colormap = varargin{1};
%colormap = colormap{1};
[maxImages, ~] = size(handles.colormap);

%channellabels = {'R', 'G', 'B'};
for k=1:maxImages
    imagelabels{k} = ['Channel: ', num2str(k)];
end

%figHandles=findobj('Type','figure');

w = 180 + 30*3;
h = 20 + 20*maxImages;

set(handles.rgbselect_figure, 'Position', [((resX/2)-(w/2)) ((resY/2)-(h/2)) w h])

% for n=1:3
%     handles.channeltxt(n) = uicontrol('Style','text','Units','pixels','HorizontalAlignment','left','Position',[100+(n-1)*30 h-15 20 15],'String',channellabels(n));
% end

for k=1:maxImages
%   handles.spectracheckbox(k) = uicontrol('Style','checkbox','Units','pixels','Position',[10 (h-30)+(k-1)*-15 20 15],'String',' ','Callback', {@spectracheckboxCallback, guidata(hObject), k});
    handles.imagelabel(k) = uicontrol('Style','text','Units','pixels','HorizontalAlignment','left','Position',[25 (h-30)+(k-1)*-15 75 15],'String',imagelabels(k));
    % 
    % for n=1:3
    %     if(handles.colormap(k,n))
    %         handles.matrixcheckbox(k, n) = uicontrol('Style','checkbox','Units','pixels','HorizontalAlignment','center','Position',[100+(n-1)*30 (h-30)+(k-1)*-15 20 15],'String','','Value',1);
    %     else
    %         handles.matrixcheckbox(k, n) = uicontrol('Style','checkbox','Units','pixels','HorizontalAlignment','center','Position',[100+(n-1)*30 (h-30)+(k-1)*-15 20 15],'String','','Value',0);
    %     end
    % end
    handles.colorpickers(k) = uicolorpicker(handles.rgbselect_figure,'Position',[100 (h-30)+(k-1)*-15 20 15],'Value',[0 0 0]);
    
end

handles.OKbutton = uicontrol('Style','pushbutton','Units','pixels','HorizontalAlignment','center','Position',[w-50 5 50 20],'String','OK','Callback', {@(hObject,eventdata)RGBmap('OKbutton_Callback',hObject,eventdata,guidata(hObject))});



guidata(hObject, handles);
uiwait(hObject);
end



function OKbutton_Callback(hObject, eventdata, handles, varargin)
%colormap = colormap{1};
[maxImages, ~] = size(handles.colormap);
    
    for k=1:maxImages
        for n=1:3
            if(get(handles.matrixcheckbox(k, n),'Value')==1)
                handles.selected(k, n) = 1;
            else
                handles.selected(k, n) = 0;
            end
        end
    end

    if(sum(sum(handles.selected == handles.colormap)) == numel(handles.colormap) )
        handles.isEdit = false;
    else
        handles.isEdit = true;
    end
%     pos = get(handles.rgbselect_figure, 'Position');
%     
%     disp(pos)
    
    guidata(hObject, handles);
    close
    
end

function varargout = RGBmap_OutputFcn(hObject, eventdata, handles) 
%get(handles.rgbselect_figure, 'Position')
%varargout{1} = handles.output;

varargout{1} = handles.selected;
varargout{2} = handles.isEdit;
delete(hObject);
end

function rgbselect_figure_CloseRequestFcn(hObject, eventdata, handles)
    if isequal(get(hObject, 'waitstatus'), 'waiting')
        uiresume(hObject);
    else
        delete(hObject);
    end
end
