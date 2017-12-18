function varargout = SSPM(varargin)
% SSPM M-file for SSPM.fig
%      SSPM, by itself, creates a new SSPM or raises the existing
%      singleton*.
%
%      H = SSPM returns the handle to a new SSPM or the handle to
%      the existing singleton*.
%
%      SSPM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SSPM.M with the given input arguments.
%
%      SSPM('Property','Value',...) creates a new SSPM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SSPM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SSPM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SSPM

% Last Modified by GUIDE v2.5 12-Sep-2014 10:52:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SSPM_OpeningFcn, ...
                   'gui_OutputFcn',  @SSPM_OutputFcn, ...
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


% --- Executes just before SSPM is made visible.
function SSPM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SSPM (see VARARGIN)




set(gcf,'Color',[0.84,0.84,0.84]);
set(gcf, 'Units','characters');
set(gcf, 'Name','SSPM','Position',[10, 60, 33, 28]);
set(handles.text1,'BackgroundColor',[0.84,0.84,0.84]);
set(handles.text2,'BackgroundColor',[0.84,0.84,0.84]);
set(handles.MAGEE,'BackgroundColor',[0.84,0.84,0.84]);
set(handles.FADTTS,'BackgroundColor',[0.84,0.84,0.84]);
set(handles.FMPM,'BackgroundColor',[0.84,0.84,0.84]);

% Choose default command line output for SSPM
handles.output = hObject;

im=imread('uncBIASlogo.jpg');
handles.h1 = image(im,'Parent',handles.axes1);
set(handles.axes1,'Visible','off')  


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SSPM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SSPM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in MAGEE.
function MAGEE_Callback(hObject, eventdata, handles)
% hObject    handle to MAGEE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath(pwd,'MAGEE');
mageegui;

% --- Executes on button press in FADTTS.
function FADTTS_Callback(hObject, eventdata, handles)
% hObject    handle to FADTTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath(pwd,'FADTTS');
FADTTS_GUI;


% --------------------------------------------------------------------
function MAGEE_menu_Callback(hObject, eventdata, handles)
% hObject    handle to MAGEE_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hint: get(hObject,'Value') returns toggle state of FMPM


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FMPM_menu_Callback(hObject, eventdata, handles)
% hObject    handle to FMPM_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in FMPM.
function FMPM_Callback(hObject, eventdata, handles)
% hObject    handle to FMPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

addpath(pwd,'FMPM');
FMPMgui;
