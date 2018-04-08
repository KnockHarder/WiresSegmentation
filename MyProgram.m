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

% Last Modified by GUIDE v2.5 07-Apr-2018 12:56:09

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
handles.rstImg = [];

set(handles.contrastEnhance,'enable','off');

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
catch ME
    set( handles.currentState, 'String', ...
        strcat('Failed to load,erro code: ',ME.identifier ) );
    return;
end

% inImg = imadjust(inImg);
handles.colImg = colImg;
handles.inImg = inImg;
handles.enImg = inImg;
handles.rstImg = [];

set(handles.disOption,'Value',1);
disOption_Callback(hObject, eventdata, handles);
set( handles.contrastEnhance, 'enable', 'on' );
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
resize = get( handles.resizeSlider, 'value' );
colImg = handles.colImg;
enImg = handles.enImg;
rstImg = handles.rstImg;

switch  v
    case    1    % show original image
        if ~isempty( colImg )
            colImg = imresize( colImg , resize );
            imshow( colImg, 'parent', handles.figAxis );
        else
            set( handles.currentState, 'String', 'No image selected');
        end
    case    2    % show contrast enhancement image
        if ~isempty( enImg )
            enImg = imresize( enImg, resize );
            imshow( enImg, 'parent', handles.figAxis );
        else
            set( handles.currentState, 'String', ...
                'Not contrast-enhanced yet' );
        end
    case    3    % show result image
        if ~isempty( rstImg )
            rstImg = imresize( rstImg , resize );
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


% --- Executes on slider movement.
function resizeSlider_Callback(hObject, eventdata, handles)
% hObject    handle to resizeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
value = get( handles.resizeSlider,'value' );
value = round( value * 10 ) / 10;
set( handles.resizeSlider, 'value', value );
set( handles.resizeValue, 'String',  value );
disOption_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function resizeSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resizeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
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
