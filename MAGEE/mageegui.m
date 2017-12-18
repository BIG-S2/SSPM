function varargout = mageegui(varargin)
% mageegui M-file for mageegui.fig
%      mageegui, by itself, creates a new mageegui or raises the existing
%      singleton*.
%
%      H = mageegui returns the handle to a new mageegui or the handle to
%      the existing singleton*.
%
%      mageegui('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in mageegui.M with the given input arguments.
%
%      mageegui('Property','Value',...) creates a new mageegui or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mageegui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mageegui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mageegui

% Last Modified by GUIDE v2.5 05-Sep-2014 16:37:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mageegui_OpeningFcn, ...
                   'gui_OutputFcn',  @mageegui_OutputFcn, ...
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


% --- Executes just before mageegui is made visible.
function mageegui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mageegui (see VARARGIN)
clr=[0.84,0.84,0.84];
set(gcf,'Color',clr);
set(gcf, 'Units','characters');
set(gcf, 'Name','MAGEE','Position',[60, 60, 108, 31]);
set(handles.RUN,'BackgroundColor',clr);
set(handles.ImShow,'BackgroundColor',clr);
set(handles.clearbut,'BackgroundColor',clr);
set(handles.InputImageFilesDir,'BackgroundColor',clr);
set(handles.OutputFilesDir,'BackgroundColor',clr);
set(handles.text1,'BackgroundColor',clr);
set(handles.text2,'BackgroundColor',clr);
set(handles.text3,'BackgroundColor',clr);
set(handles.text4,'BackgroundColor',clr);
set(handles.text5,'BackgroundColor',clr);
set(handles.text6,'BackgroundColor',clr);
set(handles.text7,'BackgroundColor',clr);
set(handles.text8,'BackgroundColor',clr);
set(handles.uipanel1,'BackgroundColor',clr);
set(handles.uipanel2,'BackgroundColor',clr);
set(handles.uipanel3,'BackgroundColor',clr);
set(handles.TimePointFn,'BackgroundColor',clr);
set(handles.DesignMatrixFn,'BackgroundColor',clr);
set(handles.MaskFn,'BackgroundColor',clr);
set(handles.ConsMatrixFn,'BackgroundColor',clr);
set(handles.ImageFn,'BackgroundColor',clr);

% Choose default command line output for mageegui
handles.output = hObject;

set(hObject, 'Name','MAGEE for Windows (64-bit) V1.1');

im=imread('uncLogo.jpg');
handles.h1 = image(im,'Parent',handles.axes1);
set(handles.axes1,'Visible','off')  


im=imread('uncbiaslogo.jpg');
handles.h2 = image(im,'Parent',handles.axes2);
set(handles.axes2,'Visible','off')  

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mageegui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mageegui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in RUN.
function RUN_Callback(hObject, eventdata, handles)
% hObject    handle to RUN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Get the five numbers from handles %%%%%%%%%%%%%%%%%%%%

    Dim_Response_str = get(handles.Dim_Response,'String');
    DimSPD=str2double(Dim_Response_str);
    if DimSPD<=0 || mod(DimSPD,1)~=0
       warndlg('Please enter a positive integer in the "Dimension of Response" box.')
        return;
    end
    num_Subjects_str = get(handles.Num_Subject,'String');
    num_Subjects=str2double(num_Subjects_str);
    if num_Subjects<=0 || mod(num_Subjects,1)~=0
        warndlg('Please enter a positive integer in the "Number of Sujects" box.')
         return;
    end
        num_cov_str = get(handles.Num_covariates,'String');
        num_cov=str2double(num_cov_str);
    if num_cov<=0 || mod(num_cov,1)~=0
        warndlg('Please enter a positive integer in the "Number of Covariates" box.')
         return;
    end
        MaxTimePoints_str = get(handles.Num_TimePoints,'String');
        MaxTimePoints=str2double(MaxTimePoints_str);
    if MaxTimePoints<=0 || mod(MaxTimePoints,1)~=0
        warndlg('Please enter a positive integer in the "Maximum Number of Time Points" box.')
        return;
    end
       num_RowofConstrM_str = get(handles.num_RowofConstrM,'String');
       No_Row=str2double(num_RowofConstrM_str);
    if No_Row<=0 || mod(No_Row,1)~=0              
        warndlg('Please enter a positive integer in the "Number of Rows of Constraint Matrix" box.')
        return;
    end
     num_iteration_str = get(handles.NumberofIteration,'String');
     num_iteration=str2double(num_iteration_str);
    if num_iteration<=0 || mod(num_iteration,1)~=0
        warndlg('Please enter a positive integer in the "Number of iterations" box.')
         return;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Draw out filenames from handles
        
       TimePointsFname=get(handles.TimePointFn,'UserData');
       if isempty(TimePointsFname)
           warndlg('Please enter a filename by pressing "Filename Containing Time Points" button.')
       return;
       end
       
       DesignMatrixFname=get(handles.DesignMatrixFn,'UserData');
       if isempty(DesignMatrixFname)
           warndlg('Please enter a filename by pressing "Filename Containing Design Matrix" button.')
       return;
       end
       
       MaskFname=get(handles.MaskFn,'UserData');
       if isempty(MaskFname)
           warndlg('Please enter a filename by pressing "Mask Image Name" button.')
       return;
       end
       
       ConsMatrixFname=get(handles.ConsMatrixFn,'UserData');
       if isempty(ConsMatrixFname)
           warndlg('Please enter a filename by pressing "Filename Containing Constraint Matrix" button.')
       return;
       end
       
       
       InputImageFn=get(handles.ImageFn,'UserData');
       if isempty(InputImageFn)
           warndlg('Please enter a filename by pressing "Filename Containing Images" button')
       return;
       end
                     
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Draw input image directory from handle      
  
       InputImageFilesDirectory= get(handles.InputImageFilesDir,'UserData');
       if isempty(InputImageFilesDirectory)
           warndlg('Please select the input-image folder by pressing "Input Image Folder" button')
       return;
       end         
        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Draw output image directory from handle
  
      OutputFilesDirectory= get(handles.OutputFilesDir,'UserData'); %%% filename without "\" at the end
      if isempty(OutputFilesDirectory)
           warndlg('Please select the input-image folder by pressing "Output Image Folder" button')
      return;
      end
     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf('Dimension of response:                    %s\n',Dim_Response_str);
  fprintf('Number of sujects:                        %s\n',num_Subjects_str);
  fprintf('Number of covariates:                     %s\n',num_cov_str);
  fprintf('Maximum number of time points:            %s\n',MaxTimePoints_str); 
  fprintf('Number of rows of constraint matrix:      %s\n',num_RowofConstrM_str);
  fprintf('Number of iterations:                     %d\n\n',num_iteration);
  fprintf('Filename containing time points:\n%s\n',TimePointsFname);
  fprintf('Filename containing design matrix:\n%s\n',DesignMatrixFname);
  fprintf('Mask image name: \n%s\n',MaskFname);
  fprintf('Filename containing constraint matrix:\n%s\n',ConsMatrixFname);
  fprintf('Filename containing input images:\n%s\n\n',InputImageFn);
  fprintf('Input image file directory:\n%s\n',InputImageFilesDirectory);
  fprintf('Output image file directory:\n%s\n\n',OutputFilesDirectory);
    
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Mask_im=load_nii(MaskFname);
  Mask_matrix=Mask_im.img;
  [dimX,dimY,dimZ]=size(Mask_matrix);
  voxels=dimX*dimY*dimZ;
  [path,name,ext] = fileparts(MaskFname);  
  MaskImName=strcat(path,'\',name,'.img');
   
  fprintf('The mask image size : [%d, %d, %d]\n',dimX,dimY,dimZ);
  fprintf('Number of voxels: %d\n', voxels); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 InputImageFilesFullName = strcat(InputImageFilesDirectory,'\InputImageFilesname.txt');
 fprintf('InputImageFilesFullName:%s\n',InputImageFilesFullName);
 fid=fopen(InputImageFilesFullName,'wt');
 fid1=fopen(InputImageFn,'rt');
 tline = fgetl(fid1);
 ii=1;
 index=0;
while ischar(tline)
   str=strcat(InputImageFilesDirectory,'\',tline,'.img');
   fprintf(fid,'%s\n',str);
   im=load_nii(str);
   matr=im.img;
   [dimxr,dimyr,dimzr]=size(matr);
   fprintf('The %dth real image data size : [%d, %d, %d]\n',ii,dimxr,dimyr,dimzr);
   
   if dimX~=dimxr ||dimY~=dimyr || dimZ~=dimzr
      fprintf('The %dth real image data size does not match those of the mask image!\n',ii);
      index=1;
      break;
   end
   tline = fgetl(fid1);   
   ii=ii+1;   
end
fclose(fid1);
fclose(fid);
clear matr
if index==0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%% mex MAGEE.cpp  predef.cpp newdelete.cpp 
%%% compile mex file MAGEE main program
fprintf('\nData processing, please wait...\n');
MAGEE(DimSPD,num_Subjects,num_cov,MaxTimePoints,TimePointsFname,DesignMatrixFname,...
      MaskImName,InputImageFilesFullName,No_Row,ConsMatrixFname,OutputFilesDirectory,dimX,dimY,dimZ,num_iteration);

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% transform output data into imge file
 fprintf('\nTransforming output data into imges,please wait...\n\n');
 OutputImageMatrixFilesFullName=strcat(OutputFilesDirectory,'\OutputFilenames.txt'); 
 fid=fopen(OutputImageMatrixFilesFullName,'rt');
 Matr_temp=zeros(dimX,dimY,dimZ,'single');
 tline = fgetl(fid);
 mm=1;
 while ischar(tline)
   fidinside=fopen(tline,'r');
   for kk=1:dimZ
       for jj=1:dimY
           Matr_temp(:,jj,kk)=fread(fidinside,dimX,'float32');
       end
   end
   fclose(fidinside);
   fprintf('%dth output images is being created...\n',mm); 
   [Path,Name]=fileparts(tline);
   imName=fullfile(Path,Name);  
   imName=strcat(imName,'.nii');        
   IM=make_nii(Matr_temp,[1,1,1],[1.0,1.0,1.0],16,'');
   save_nii(IM,imName);
   tline = fgetl(fid);
   mm=mm+1;
   Matr_temp=zeros(dimX,dimY,dimZ,'single');
 end
 fclose(fid); 
 
  str=strcat(OutputFilesDirectory,'\*.dat');
  delete(str)
  str=strcat(OutputFilesDirectory,'\OutputFilenames.txt');
  delete(str)
 
toc
ETstr=num2str(toc);
set(handles.ElapsedTime,'String',ETstr);
fprintf('+++++++Your job is finished!!+++++++\n');
 end
guidata(hObject,handles);


function Dim_Response_Callback(hObject, eventdata, handles)
% hObject    handle to Dim_Response (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

% Hints: get(hObject,'String') returns contents of Dim_Response as text
%        str2double(get(hObject,'String')) returns contents of Dim_Response as a double


% --- Executes during object creation, after setting all properties.
function Dim_Response_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dim_Response (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Num_Subject_Callback(hObject, eventdata, handles)
% hObject    handle to Num_Subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);


% Hints: get(hObject,'String') returns contents of Num_Subject as text
%        str2double(get(hObject,'String')) returns contents of Num_Subject as a double


% --- Executes during object creation, after setting all properties.
function Num_Subject_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Num_Subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Num_covariates_Callback(hObject, eventdata, handles)
% hObject    handle to Num_covariates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

% Hints: get(hObject,'String') returns contents of Num_covariates as text
%        str2double(get(hObject,'String')) returns contents of Num_covariates as a double


function Num_TimePoints_Callback(hObject, eventdata, handles)
% hObject    handle to Num_TimePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);

% Hints: get(hObject,'String') returns contents of Num_TimePoints as text
%        str2double(get(hObject,'String')) returns contents of Num_TimePoints as a double


% --- Executes during object creation, after setting all properties.
function Num_TimePoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Num_TimePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_RowofConstrM_Callback(hObject, eventdata, handles)
% hObject    handle to num_RowofConstrM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject,handles);


% Hints: get(hObject,'String') returns contents of num_RowofConstrM as text
%        str2double(get(hObject,'String')) returns contents of num_RowofConstrM as a double


% --- Executes during object creation, after setting all properties.
function num_RowofConstrM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_RowofConstrM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on MaskFnInput and none of its controls.
function MaskFnInput_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to MaskFnInput (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


function ElapsedTime_Callback(hObject, eventdata, handles)
% hObject    handle to ElapsedTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 
% Hints: get(hObject,'String') returns contents of ElapsedTime as text
%        str2double(get(hObject,'String')) returns contents of ElapsedTime as a double


% --- Executes during object creation, after setting all properties.
function ElapsedTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ElapsedTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on ConsMatrixFnInput and none of its controls.
function ConsMatrixFnInput_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to ConsMatrixFnInput (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Num_covariates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Num_covariates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clearbut.
function clearbut_Callback(hObject, eventdata, handles)
% hObject    handle to clearbut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 set(handles.Dim_Response,'String','');   
 set(handles.Num_Subject,'String','');
 set(handles.Num_covariates,'String','');
 set(handles.Num_TimePoints,'String','');
 set(handles.num_RowofConstrM,'String','');
 set(handles.NumberofIteration,'String','');
 set(handles.TimePointFn,'UserData','');
 set(handles.DesignMatrixFn,'UserData','');
 set(handles.MaskFn,'UserData','');
 set(handles.ConsMatrixFn,'UserData','');
 set(handles.ImageFn,'UserData','');
 set(handles.ImShow,'UserData','');
 set(handles.ElapsedTime,'String','');
 set(handles.InputImageFilesDir,'UserData',''); 
 set(handles.OutputFilesDir,'UserData',''); 
 guidata(hObject,handles);
                          
 
% --- Executes on button press in TimePointFn.
function TimePointFn_Callback(hObject, eventdata, handles)
% hObject    handle to TimePointFn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
     DesignMatrixFname=get(handles.DesignMatrixFn,'UserData');
     MaskFname=get(handles.MaskFn,'UserData');
     ConsMatrixFname=get(handles.ConsMatrixFn,'UserData');
     InputImageFn=get(handles.ImageFn,'UserData');
     if ~isempty(DesignMatrixFname)
        [path,name,ext] = fileparts(DesignMatrixFname); 
     elseif ~isempty(MaskFname)
        [path,name,ext] = fileparts(MaskFname); 
     elseif ~isempty(ConsMatrixFname)
        [path,name,ext] = fileparts(ConsMatrixFname); 
     elseif ~isempty(InputImageFn)
        [path,name,ext] = fileparts(InputImageFn); 
     else
         path=pwd;
     end       
   [FileName,PathName] = uigetfile(fullfile(path,'*.txt'),'Select a file');
    if FileName~=0
    str=strcat(PathName,FileName);
    set(handles.TimePointFn,'UserData',str);
    end
end
guidata(hObject,handles);


% Hint: get(hObject,'Value') returns toggle state of TimePointFn


% --- Executes on button press in DesignMatrixFn.
function DesignMatrixFn_Callback(hObject, eventdata, handles)
% hObject    handle to DesignMatrixFn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max') 
     TimePointsFname=get(handles.TimePointFn,'UserData');
     MaskFname=get(handles.MaskFn,'UserData');
     ConsMatrixFname=get(handles.ConsMatrixFn,'UserData');
     InputImageFn=get(handles.ImageFn,'UserData');
     if ~isempty(TimePointsFname)
        [path,name,ext] = fileparts(TimePointsFname); 
     elseif ~isempty(MaskFname)
        [path,name,ext] = fileparts(MaskFname); 
     elseif ~isempty(ConsMatrixFname)
        [path,name,ext] = fileparts(ConsMatrixFname); 
     elseif ~isempty(InputImageFn)
        [path,name,ext] = fileparts(InputImageFn); 
     else
         path=pwd;
     end       
   [FileName,PathName] = uigetfile(fullfile(path,'*.txt'),'Select a file');
   if FileName~=0
    str=strcat(PathName,FileName);
    set(handles.DesignMatrixFn,'UserData',str);
   end
end
    guidata(hObject,handles);

% Hint: get(hObject,'Value') returns toggle state of DesignMatrixFn


% --- Executes on button press in MaskFn.
function MaskFn_Callback(hObject, eventdata, handles)
% hObject    handle to MaskFn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    workdir=pwd;
     TimePointsFname=get(handles.TimePointFn,'UserData');
     DesignMatrixFname=get(handles.DesignMatrixFn,'UserData');
     ConsMatrixFname=get(handles.ConsMatrixFn,'UserData');
     InputImageFn=get(handles.ImageFn,'UserData');
     if ~isempty(TimePointsFname)
        [path,name,ext] = fileparts(TimePointsFname); 
     elseif ~isempty(DesignMatrixFname)
        [path,name,ext] = fileparts(DesignMatrixFname); 
     elseif ~isempty(ConsMatrixFname)
        [path,name,ext] = fileparts(ConsMatrixFname); 
     elseif ~isempty(InputImageFn)
        [path,name,ext] = fileparts(InputImageFn); 
     else
         path=workdir;
     end
     
     if ~isempty(path)
         cd(path);  
         [FileName,PathName] = uigetfile({'*.hdr';'*.nii'},'Select a file');
         cd(workdir)
     else
         [FileName,PathName] = uigetfile({'*.hdr';'*.nii'},'Select a file');
     end
     
    if FileName~=0
    str=strcat(PathName,FileName);
    set(handles.MaskFn,'UserData',str);
    end
end
guidata(hObject,handles);

% Hint: get(hObject,'Value') returns toggle state of MaskFn


% --- Executes on button press in ConsMatrixFn.
function ConsMatrixFn_Callback(hObject, eventdata, handles)
% hObject    handle to ConsMatrixFn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
     TimePointsFname=get(handles.TimePointFn,'UserData');
     DesignMatrixFname=get(handles.DesignMatrixFn,'UserData');
     MaskFname=get(handles.MaskFn,'UserData');
     InputImageFn=get(handles.ImageFn,'UserData');
     if ~isempty(TimePointsFname)
        [path,name,ext] = fileparts(TimePointsFname); 
     elseif ~isempty(DesignMatrixFname)
        [path,name,ext] = fileparts(DesignMatrixFname); 
     elseif ~isempty(MaskFname)
        [path,name,ext] = fileparts(MaskFname); 
     elseif ~isempty(InputImageFn)
        [path,name,ext] = fileparts(InputImageFn); 
     else
         path=pwd;
     end       
   [FileName,PathName] = uigetfile(fullfile(path,'*.txt'),'Select a file');
    if FileName~=0
    str=strcat(PathName,FileName);
    set(handles.ConsMatrixFn,'UserData',str);
    end
end
guidata(hObject,handles);



% Hint: get(hObject,'Value') returns toggle state of ConsMatrixFn


% --- Executes on button press in ImageFn.
function ImageFn_Callback(hObject, eventdata, handles)
% hObject    handle to ImageFn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    TimePointsFname=get(handles.TimePointFn,'UserData');
     DesignMatrixFname=get(handles.DesignMatrixFn,'UserData');
     MaskFname=get(handles.MaskFn,'UserData');
     ConsMatrixFname=get(handles.ConsMatrixFn,'UserData');
     if ~isempty(TimePointsFname)
        [path,name,ext] = fileparts(TimePointsFname); 
     elseif ~isempty(DesignMatrixFname)
        [path,name,ext] = fileparts(DesignMatrixFname); 
     elseif ~isempty(MaskFname)
        [path,name,ext] = fileparts(MaskFname); 
     elseif ~isempty(ConsMatrixFname)
        [path,name,ext] = fileparts(ConsMatrixFname); 
     else
         path=pwd;
     end       
   [FileName,PathName] = uigetfile(fullfile(path,'*.txt'),'Select a file');
    if FileName~=0
    str=strcat(PathName,FileName);
    set(handles.ImageFn,'UserData',str);
    end
end
guidata(hObject,handles);


% Hint: get(hObject,'Value') returns toggle state of ImageFn


% --- Executes on button press in ImShow.
function ImShow_Callback(hObject, eventdata, handles)
% hObject    handle to ImShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uigetfile('*.nii','Select a text file');
if FileName~=0
    str=strcat(PathName,FileName,'.nii');
    set(handles.ImShow,'UserData',str);
    IM=load_nii(str);
    view_nii(IM);
else
    [FileName,PathName] = uigetfile('*.hdr','Select a text file');
    str=strcat(PathName,FileName);
    set(handles.ImShow,'UserData',str);
    IM=load_nii(str);
    view_nii(IM);
end
guidata(hObject,handles);


% --- Executes on key press with focus on DesignMatrixFn and none of its controls.
function DesignMatrixFn_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to DesignMatrixFn (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in InputImageFilesDir.
function InputImageFilesDir_Callback(hObject, eventdata, handles)
% hObject    handle to InputImageFilesDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
FilesDir= uigetdir('','Select a folder containing input images');
    if FilesDir~=0
       set(handles.InputImageFilesDir,'UserData',FilesDir);
    end
end
guidata(hObject,handles);

% Hint: get(hObject,'Value') returns toggle state of InputImageFilesDir


% --- Executes on button press in OutputFilesDir.
function OutputFilesDir_Callback(hObject, eventdata, handles)
% hObject    handle to OutputFilesDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
FilesDir= uigetdir('','Select a folder to output images');
    if FilesDir~=0
       set(handles.OutputFilesDir,'UserData',FilesDir);
    end
end
guidata(hObject,handles);

% Hint: get(hObject,'Value') returns toggle state of OutputFilesDir



function NumberofIteration_Callback(hObject, eventdata, handles)
% hObject    handle to NumberofIteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumberofIteration as text
%        str2double(get(hObject,'String')) returns contents of NumberofIteration as a double


% --- Executes during object creation, after setting all properties.
function NumberofIteration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumberofIteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function text1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2
