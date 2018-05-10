function varargout = MyProgram(varargin)
% UNTITLED MATLAB code for untitled.fig
%      UNTITLED, by itself, creates a new UNTITLED or raises the existing
%      singleton*.
%
%      H = UNTITLED returns the handle to a new UNTITLED or the handle to
%      the existing singleton*.
%
%      UNTITLED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED.M with the given input arguments.
%
%      UNTITLED('Property','Value',...) creates a new UNTITLED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled

% Last Modified by GUIDE v2.5 10-May-2018 23:12:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before untitled is made visible.
function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled (see VARARGIN)

% Choose default command line output for untitled
handles.output = hObject;

% All the user data goes here
handles.colImg = [];
handles.inImg = [];
handles.enImg = [];
handles.labelImg = [];
handles.cropRec = [];
handles.rstImg = [];

set(handles.contrastEnhance,'enable','off');
set( handles.doCrop, 'enable', 'off' );
set( handles.localGrowing, 'enable', 'off' );
set( handles.labelMenu, 'enable', 'off' );

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadImg.
function loadImg_Callback(hObject, eventdata, handles)
% hObject    handle to loadImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    [colImg,inImg,~, ~] = StdIP.readImg2D(-1);
    inImg = inImg/max(inImg(:));
catch 
    set( handles.currentState, 'String', ...
        'Image is not loaded rightly, check and try again' );
    return;
end

% inImg = imadjust(inImg);
handles.colImg = colImg;
handles.inImg = inImg;
handles.enImg = inImg;
handles.labelImg = [];
handles.rstImg = [];

[m,n] = size( inImg );
handles.cropRec = [1, 1, n, m];

set(handles.disOption,'Value',1);
disOption_Callback(hObject, eventdata, handles);
set( handles.contrastEnhance, 'enable', 'on' );
set( handles.doCrop, 'enable', 'off' );
set( handles.localGrowing, 'enable', 'off' );
set( handles.labelMenu, 'enable', 'off' );
set( handles.currentState, 'String','Ready for processing' );
guidata(hObject, handles);



% --- Executes on selection change in disOption.
function disOption_Callback(hObject, eventdata, handles)
% hObject    handle to disOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns disOption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from disOption
v = get(handles.disOption,'Value');
colImg = handles.colImg;
enImg = handles.enImg;
rstImg = handles.rstImg;

switch  v
    case    1    % show original image
        if ~isempty( colImg )
            imshow( colImg, 'parent', handles.figAxis );
        else
            set( handles.currentState, 'String', 'No image selected');
        end
    case    2    % show contrast enhancement image
        if ~isempty( enImg )
            imshow( enImg, 'parent', handles.figAxis );
        else
            set( handles.currentState, 'String', ...
                'Not contrast-enhanced yet' );
        end
    case    3    % show result image
        if ~isempty( rstImg )
            imshow( rstImg, 'parent', handles.figAxis );
        else
            set( handles.currentState, 'String', 'No result image found');
        end
    otherwise
        set(handles.currentState,'String','Not yet implemented!');
end


% --- Executes during object creation, after setting all properties.
function disOption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to disOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in contrastEnhance.
function contrastEnhance_Callback(hObject, eventdata, handles)
% hObject    handle to contrastEnhance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if( get( handles.vesselBright, 'value' ) )
    handles.enImg = WD.contrastEnhance( handles.inImg );
else
    handles.enImg = WD.contrastEnhance( 1 - handles.inImg );
end

set(handles.disOption,'Value',2);
disOption_Callback(hObject, eventdata, handles);
set(handles.currentState,'String','Contrast Enhancement is done');
set( handles.doCrop, 'enable', 'on' );
set( handles.localGrowing, 'enable', 'on' );
guidata(hObject, handles);


% --- Executes on button press in vesselBright.
function vesselBright_Callback(hObject, eventdata, handles)
% hObject    handle to vesselBright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vesselBright


% --- Executes on button press in vesselDark.
function vesselDark_Callback(hObject, eventdata, handles)
% hObject    handle to vesselDark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vesselDark


% --- Executes on button press in localGrowing.
function localGrowing_Callback(hObject, eventdata, handles)
% hObject    handle to localGrowing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.currentState,'String',...
    'A little time is needed to do localGrowing, be paitient');
guidata(hObject, handles);
drawnow;

str_iter = get( handles.Iteration, 'String' );
iter = str2double( str_iter );
if  isnan(iter) || iter < 0
    set( handles.currentState, 'String', ...
        '"iteration" should be a positive integer' );
    return;
end
[m,n] = size( handles.inImg );
rec = handles.cropRec;
mask = zeros( m, n );
mask( rec(2):rec(2)+rec(4)-1, rec(1):rec(1)+rec(3)-1 ) = 1;
inImg = handles.inImg .* mask;
enImg = handles.enImg .* mask;
labelImg = WD.localGrowing( inImg, enImg, iter );
% assert( max(max(labelImg)) == length(unique(labelImg)) - 1 );

colorI = label2rgb( labelImg, @jet, [0, 0, 0] );
colorI = im2double( colorI );
mask = labelImg == 0;
bkgImg = handles.inImg .* mask * 0.5;
grayImg = zeros( m,n,3 );
for i = 1 : 3
    grayImg(:,:,i) = bkgImg;
end
rstImg = imadd( grayImg, colorI );
handles.labelImg = labelImg;
handles.rstImg = im2uint8( rstImg );

set(handles.disOption,'Value',3);
disOption_Callback(hObject, eventdata, handles);
set(handles.currentState,'String','Vessel Growing is done');
set(handles.labelMenu, 'String', 0:max(max(labelImg)) );
set(handles.labelMenu, 'value', 1 );
set(handles.labelMenu, 'enable', 'on' );
guidata(hObject, handles);


% --- Executes on selection change in labelMenu.
function labelMenu_Callback(hObject, eventdata, handles)
% hObject    handle to labelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns labelMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from labelMenu
label = get( handles.labelMenu, 'value') - 1;

colorI = label2rgb( handles.labelImg, @jet, [0, 0, 0] );
colorI = im2double( colorI );
labelImg = handles.labelImg;
mask = labelImg == 0;
bkgImg = handles.inImg .* mask * 0.5;
[m,n] = size(labelImg);
grayImg = zeros( m,n,3 );
for i = 1 : 3
    grayImg(:,:,i) = bkgImg;
end
rstImg = imadd( grayImg, colorI );
handles.rstImg = im2uint8( rstImg );
if( label ~= 0 )
    ind = find( labelImg == label );
    [X,Y] = ind2sub( size(labelImg), ind );
    
    for i = 1 : length(ind)
        handles.rstImg( X(i), Y(i), : ) = [255, 0, 0];
        handles.rstImg( X(i)-1, Y(i), : ) = [255, 0, 0];
        handles.rstImg( X(i)+1, Y(i), : ) = [255, 0, 0];
        handles.rstImg( X(i), Y(i)-1, : ) = [255, 0, 0];
        handles.rstImg( X(i), Y(i)+1, : ) = [255, 0, 0];
    end
end
set(handles.disOption,'Value',3);
disOption_Callback(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function labelMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to labelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in doCrop.
function doCrop_Callback(hObject, eventdata, handles)
% hObject    handle to doCrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disOption_Callback(hObject, eventdata, handles );
set( handles.localGrowing, 'enable', 'off' );
set( handles.currentState, 'String', ...
    'If everything is done, please double-left-click' );
guidata( hObject, handles );

disImg = getimage( handles.figAxis );
if isempty(disImg)
    return;
end

[cropImg,rec] = imcrop( disImg );
if ~isempty(rec)    
    handles.cropRec = round( rec );
    imshow( cropImg, 'parent', handles.figAxis );
    set( handles.currentState, 'String', 'Area select is done' );
else
    set( handles.currentState, 'String', 'Area select has been canceled' );
end

set( handles.localGrowing, 'enable', 'on' );
guidata( hObject, handles );



function Iteration_Callback(hObject, eventdata, handles)
% hObject    handle to Iteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Iteration as text
%        str2double(get(hObject,'String')) returns contents of Iteration as a double


% --- Executes during object creation, after setting all properties.
function Iteration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Iteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
