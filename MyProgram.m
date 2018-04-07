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
handles.CEImg = [];
handles.BEImg = [];
handles.rstImg = [];

set(handles.contrastEnhance,'enable','off');
set(handles.binaryEnhance,'enable','off');

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
handles.CEImg = inImg;
handles.BEImg = [];
handles.rstImg = [];

set(handles.disOption,'Value',1);
disOption_Callback(hObject, eventdata, handles);
set( handles.contrastEnhance, 'enable', 'on' );
set( handles.binaryEnhance, 'enable', 'on' );
set( handles.currentState, 'String','Ready for processing' );
guidata(hObject, handles);


% --- Executes on button press in binaryEnhance.
function binaryEnhance_Callback(hObject, eventdata, handles)
% hObject    handle to binaryEnhance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

minValue = str2double( get( handles.valueMin, 'String' ) );
maxValue = str2double(get( handles.valueMax, 'String' ) );
phases = str2double(get( handles.valueP, 'String' ) );

if( isnan(minValue) || isnan(maxValue) || isnan(phases) )
    set( handles.currentState, 'String',...
        'Values (min,max,phase) are not numbers.' );
    return;
end

minValue = round( minValue );
maxValue = round( maxValue );
phases = round( phases );
if( minValue < 0 || maxValue < 0 || minValue == maxValue ...
        || phases < 0 ||  phases > 100 )
    set( handles.currentState, 'String', ...
        'Values (min,max,phase) are out of place.' );
    return;
end

if( minValue > maxValue )
    x = minValue; minValue = maxValue; maxValue =x;
end

set( handles.valueMin, 'String', minValue );
set( handles.valueMax, 'String', maxValue );
set( handles.valueP, 'String', phases );

handles.BEImg = WD.areaCut( handles.inImg, minValue, maxValue, phases );

set(handles.disOption,'Value',3);
disOption_Callback(hObject, eventdata, handles);
set(handles.currentState,'String','Binary Enhancement is done');
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
CEImg = handles.CEImg;
BEImg = handles.BEImg;
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
        if ~isempty( CEImg )
            CEImg = imresize( CEImg, resize );
            imshow( CEImg, 'parent', handles.figAxis );
        else
            set( handles.currentState, 'String', ...
                'Not contrast-enhanced yet' );
        end
    case    3    % show binary enhancement image
        if ~isempty( BEImg )
            BEImg = imresize( BEImg, resize );
            imshow( BEImg, 'parent', handles.figAxis );
        else
            set( handles.currentState, 'String', ...
                'Not binary-enhanced yet' );
        end
    case    4    % show result image
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



function valueMin_Callback(hObject, eventdata, handles)
% hObject    handle to valueMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of valueMin as text
%        str2double(get(hObject,'String')) returns contents of valueMin as a double



% --- Executes during object creation, after setting all properties.
function valueMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to valueMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function valueMax_Callback(hObject, eventdata, handles)
% hObject    handle to valueMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of valueMax as text
%        str2double(get(hObject,'String')) returns contents of valueMax as a double


% --- Executes during object creation, after setting all properties.
function valueMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to valueMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function valueP_Callback(hObject, eventdata, handles)
% hObject    handle to valueP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of valueP as text
%        str2double(get(hObject,'String')) returns contents of valueP as a double


% --- Executes during object creation, after setting all properties.
function valueP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to valueP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function alphaSlider_Callback(hObject, eventdata, handles)
% hObject    handle to alphaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
value = get( handles.alphaSlider,'value' );
value = round( value );
set( handles.alphaSlider, 'value', value );
set( handles.alphaValue, 'String',  value );
disOption_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function alphaSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function valueAlpha_Callback(hObject, eventdata, handles)
% hObject    handle to valueAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of valueAlpha as text
%        str2double(get(hObject,'String')) returns contents of valueAlpha as a double


% --- Executes during object creation, after setting all properties.
function valueAlpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to valueAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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
    handles.CEImg = WD.contrastEnhance( handles.inImg );
else
    handles.CEImg = WD.contrastEnhance( 1 - handles.inImg );
end

set(handles.disOption,'Value',2);
disOption_Callback(hObject, eventdata, handles);
set(handles.currentState,'String','Contrast Enhancement is done');
set(handles.binaryEnhance, 'enable', 'on' );
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
