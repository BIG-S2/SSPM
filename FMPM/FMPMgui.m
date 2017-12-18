function varargout = FMPMgui(varargin)
% FMPMGUI MATLAB code for FMPMgui.fig
%      FMPMGUI, by itself, creates a new FMPMGUI or raises the existing
%      singleton*.
%
%      H = FMPMGUI returns the handle to a new FMPMGUI or the handle to
%      the existing singleton*.
%
%      FMPMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FMPMGUI.M with the given input arguments.
%
%      FMPMGUI('Property','Value',...) creates a new FMPMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FMPMgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FMPMgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FMPMgui

% Last Modified by GUIDE v2.5 08-Sep-2014 10:37:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FMPMgui_OpeningFcn, ...
                   'gui_OutputFcn',  @FMPMgui_OutputFcn, ...
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


% --- Executes just before FMPMgui is made visible.
function FMPMgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FMPMgui (see VARARGIN)

clr=[0.84,0.84,0.84];
set(gcf,'Color',clr);
set(gcf, 'Units','characters');
set(gcf, 'Name','FMPM','Position',[60, 60, 150, 17]);
set(handles.title,'BackgroundColor',clr);
set(handles.IODir,'BackgroundColor',clr);
set(handles.DataFn,'BackgroundColor',clr);
set(handles.CdesignFn,'BackgroundColor',clr);
set(handles.expvartitle,'BackgroundColor',clr);
set(handles.HT,'BackgroundColor',clr);
set(handles.SB,'BackgroundColor',clr);

set(handles.inoutput,'BackgroundColor',clr);

% Choose default command line output for mageegui
handles.output = hObject;

set(hObject, 'Name','FMPM for Windows (64-bit) V1.1');

im=imread('uncLogo.jpg');
handles.h1 = image(im,'Parent',handles.axes1);
set(handles.axes1,'Visible','off')  


im=imread('uncbiaslogo.jpg');
handles.h2 = image(im,'Parent',handles.axes2);
set(handles.axes2,'Visible','off')  

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FMPMgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FMPMgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in IODir.
function IODir_Callback(hObject, eventdata, handles)
% hObject    handle to IODir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 button_state = get(hObject,'Value');
    if button_state == get(hObject,'Max')
    IODir= uigetdir(' ','Select Input-Output folder');
        if IODir~=0
           set(handles.IODir,'UserData',IODir);
        end
    end
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of IODir


% --- Executes on button press in DataFn.
function DataFn_Callback(hObject, eventdata, handles)
% hObject    handle to DataFn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button_state = get(hObject,'Value');
    if button_state == get(hObject,'Max')
    DataFn_str = uigetfile('.mat','Select data file');
        if DataFn_str~=0
           set(handles.DataFn,'UserData',DataFn_str);
        end
    end
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of DataFn


% --- Executes on button press in CdesignFn.
function CdesignFn_Callback(hObject, eventdata, handles)
% hObject    handle to CdesignFn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CdesignFn
button_state = get(hObject,'Value');
    if button_state == get(hObject,'Max')
    CdesignFn_str = uigetfile('.mat','Select Cdesign file');
        if CdesignFn_str~=0
           set(handles.CdesignFn,'UserData',CdesignFn_str);
        end
    end
guidata(hObject,handles);


function expvarval_Callback(hObject, eventdata, handles)
% hObject    handle to expvarval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of expvarval as text
%        str2double(get(hObject,'String')) returns contents of expvarval as a double


% --- Executes during object creation, after setting all properties.
function expvarval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to expvarval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in HT.
function HT_Callback(hObject, eventdata, handles)
% hObject    handle to HT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    tic
    
    clc

    addpath('/');

    IODir= get(handles.IODir,'UserData');
       if isempty(IODir)
           warndlg('Please select the input-output folder by pressing "Input-Output" button')
          return;
       end   

    DataName= get(handles.DataFn,'UserData');
       if isempty(DataName)
           warndlg('Please select the Data file by pressing "Data" button')
          return;
       end   
       
    load(DataName)
    
    CdesignName = get(handles.CdesignFn,'UserData');
       if isempty(CdesignName)
           warndlg('Please select the Cdesign file by pressing "Cdesign" button')
          return;
       end 
    load(CdesignName)
       
    ExpVar_str = get(handles.expvarval,'String');
    ExpVar=str2double(ExpVar_str);
    
       
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      fprintf('Input-Output Directory:                   %s\n',IODir);
      fprintf('Data file name:                           %s\n',DataName);
      fprintf('Cdesign file name:                        %s\n',CdesignName);
      fprintf('ExpVar:                                   %s\n',ExpVar_str);   
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\nHypothesis Testing, please wait...\n');
 
        Xdesign=[Xdesign(:,1:2), log(Xdesign(:,3)), log(Xdesign(:,3)).^2];
        Zdesign= Xdesign(:,[1,3]);
        [NoSetup, Xdesign, Zdesign]=FMPM_ReadNormalize(Xdesign,Zdesign,Indicator, arclength);
        n=NoSetup(1);
        L0=NoSetup(2);
        p=NoSetup(3);
        sumri=NoSetup(4);
        pz=NoSetup(5);

        B0vector=zeros(size(Cdesign,1),L0);

        rng('default');
        [mh1,~,~] = FMPM_LL_wob(NoSetup,Indicator,arclength,Xdesign,Ydesign);
        [~,~,InvSigmats,~,~,ResYdesign,H0ResYdesign] = FMPM_LL_wb_H0(NoSetup,arclength,Xdesign,Ydesign,mh1,Cdesign,B0vector);
        [SigmaZeta,SigmaZetad,SigmaXi,SigmaXid] = FMPM_Cov_LeastSquare_e(NoSetup,Indicator,Zdesign,ResYdesign);
        [H0SigmaZeta,H0SigmaZetad,H0SigmaXi,H0SigmaXid] = FMPM_Cov_LeastSquare_e(NoSetup,Indicator,Zdesign,H0ResYdesign);
        [mh0] = BiSmooth_LL_bs(NoSetup,arclength,SigmaXi,pz,0);
        [SigmaXiP,SigmaXiE] = BiSmooth_LL(NoSetup,arclength,SigmaXi,SigmaXid,pz,mh0,0);
        [mh00] = BiSmooth_LL_bs(NoSetup,arclength,SigmaZeta,1,0);
        [SigmaZetaP,SigmaZetaE] = BiSmooth_LL(NoSetup,arclength,SigmaZeta, SigmaZetad,1,mh00,0);
        [mh0] = BiSmooth_LL_bs(NoSetup,arclength,H0SigmaXi,pz,0);
        [H0SigmaXiP,H0SigmaXiE] = BiSmooth_LL(NoSetup,arclength,H0SigmaXi,H0SigmaXid,pz,mh0,0);
        [mh00] = BiSmooth_LL_bs(NoSetup,arclength,H0SigmaZeta,1,0);
        [H0SigmaZetaP,H0SigmaZetaE] = BiSmooth_LL(NoSetup,arclength,H0SigmaZeta, H0SigmaZetad,1,mh00,0);
        SigmaXiE=kron(eye(L0),SigmaXiE);
        SigmaZetaE=SigmaZetaE*eye(L0);
        H0SigmaXiE=kron(eye(L0),H0SigmaXiE);
        H0SigmaZetaE=H0SigmaZetaE*eye(L0);

        [~,~,~,~,~,~,SigmaXiPCA, SigmaZetaPCA] = FMPM_FPCA_Cov(NoSetup,SigmaXiP,SigmaZetaP,SigmaXiE,SigmaZetaE,ExpVar);
        [~,~,~,~,~,~,H0SigmaXiPCA, H0SigmaZetaPCA] = FMPM_FPCA_Cov(NoSetup,H0SigmaXiP,H0SigmaZetaP,H0SigmaXiE,H0SigmaZetaE,ExpVar);
        [mh4,GCVs,vh,vYdesign,vYdesignH0, tempv, tempvH0] = FMPM_LL_Ref_wob_H0(NoSetup,Indicator,arclength,Xdesign,Zdesign,Ydesign, SigmaXiPCA,SigmaZetaPCA, H0SigmaXiPCA,H0SigmaZetaPCA);
        mh4=mh4/6;
        [efitBetas,efitBetas1,efitBetasH0,efitBetas1H0,efitYdesign,efitYdesignH0, ResYdesign, H0ResYdesign, InvXCovX]= FMPM_LL_Ref_wb_H0(NoSetup,Indicator, arclength,Xdesign,Ydesign,vYdesign,tempv,mh4,Cdesign,B0vector);
        efitBetas_2=efitBetas;

        [Lstat,Gstat] = FMPM_bstrp_stat(arclength,efitBetas,InvXCovX,Cdesign,B0vector);
        [ResEtasH0,efitEtasH0] = FMPM_LL_If_wob(NoSetup,arclength,H0ResYdesign);
        GG=1000;
        Smat=zeros(L0,GG);
        SAmat=zeros(GG,1);
        Pvalue=zeros(L0,1);
        PvalueG=0;
        for gii=1:GG
        [gYdesign] = FMPM_grs(NoSetup,Indicator,efitBetasH0,efitEtasH0,ResEtasH0,Xdesign);
        [gefitBetas,gefitBetas1]= FMPM_LL_Ref_wb(NoSetup,Indicator,arclength,Xdesign,gYdesign,Zdesign,SigmaXiPCA,SigmaZetaPCA,mh4);
        [gLstat,gGstat] = FMPM_bstrp_stat(arclength,gefitBetas,InvXCovX,Cdesign,B0vector);
        Smat(:,gii)=gLstat;
        SAmat(gii)=gGstat;
        MaxSn=max(gLstat);
        Pvalue=Pvalue+( (MaxSn*ones(L0,1))>Lstat);
        PvalueG=PvalueG+(gGstat>Gstat);
        end
        
        fprintf('P-values(L) are:\n');
        PvalL = Pvalue/GG
        
        fprintf('\n');
        fprintf('P-values(G) is:\n');
        PvalG = PvalueG/GG
        
    
        str=fullfile(IODir, 'PvalL.mat');
        save(str,'PvalL')
        
        str=fullfile(IODir, 'PvalG.mat');
        save(str,'PvalG')
                
  toc
  
  fprintf('+++++++Your job is finished!!+++++++\n');

  guidata(hObject,handles);


% --- Executes on button press in SB.
function SB_Callback(hObject, eventdata, handles)
% hObject    handle to SB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 clc
 tic
 
 addpath('/');

    IODir= get(handles.IODir,'UserData');
       if isempty(IODir)
           warndlg('Please select the input-output folder by pressing "Input-Output" button')
          return;
       end   

    DataName= get(handles.DataFn,'UserData');
       if isempty(DataName)
           warndlg('Please select the Data file by pressing "Data" button')
          return;
       end         
    load(DataName)
    
    ExpVar_str = get(handles.expvarval,'String');
    ExpVar=str2double(ExpVar_str);
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      fprintf('Input-Output Directory:                   %s\n',IODir);
      fprintf('Data file name:                           %s\n',DataName);;
      fprintf('ExpVar:                                   %s\n',ExpVar_str);   
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\nSimultaneous Confidence Interval, please wait...\n');
    
    Xdesign=[Xdesign(:,1:2), Xdesign(:,3), Xdesign(:,3).^2];
    Zdesign= Xdesign(:,[1,3]);
    [NoSetup, Xdesign, Zdesign]=FMPM_ReadNormalize(Xdesign,Zdesign,Indicator, arclength);
    n=NoSetup(1);
    L0=NoSetup(2);
    p=NoSetup(3);
    sumri=NoSetup(4);
    pz=NoSetup(5);
    GG=5000;
    rng('default');
    [mh1,GCVs,vh] = FMPM_LL_wob(NoSetup,Indicator,arclength,Xdesign,Ydesign);
    [efitBetas,efitBetas1,InvSigmats,efitYdesign,ResYdesign] = FMPM_LL_wb(NoSetup,arclength,Xdesign,Ydesign,mh1);

    [SigmaZeta,SigmaZetad,SigmaXi,SigmaXid] = FMPM_Cov_LeastSquare_e(NoSetup,Indicator,Zdesign,ResYdesign);
    [mh0] = BiSmooth_LL_bs(NoSetup,arclength,SigmaXi,pz,0);
    [SigmaXiP,SigmaXiE] = BiSmooth_LL(NoSetup,arclength,SigmaXi,SigmaXid,pz,mh0,0);
    [mh00] = BiSmooth_LL_bs(NoSetup,arclength,SigmaZeta,1,0);
    [SigmaZetaP,SigmaZetaE] = BiSmooth_LL(NoSetup,arclength,SigmaZeta, SigmaZetad,1,mh00,0);
    SigmaXiE=kron(eye(L0),SigmaXiE);
    SigmaZetaE=SigmaZetaE*eye(L0);
    [~,~,~,~,~,~,SigmaXiPCA, SigmaZetaPCA] = FMPM_FPCA_Cov(NoSetup,SigmaXiP,SigmaZetaP,SigmaXiE,SigmaZetaE,ExpVar);


    [mh4,GCVs,vh,vYdesign, tempv] = FMPM_LL_Ref_wob_SB_T(NoSetup,Indicator, arclength,Xdesign,Zdesign,Ydesign, SigmaXiPCA,SigmaZetaPCA);
    [efitBetas,efitBetas1,efitYdesign, ResYdesign,InvSigmats,Tempmat]= FMPM_LL_Ref_wb_SB_T(NoSetup,Indicator, arclength,Xdesign,Zdesign, Ydesign,vYdesign,tempv,SigmaXiPCA,SigmaZetaPCA, mh4);
    [Gvalue] = FMPM_Gval(NoSetup,Tempmat, InvSigmats, GG);

    LowerBd=zeros(p,L0,2);
    UpperBd=zeros(p,L0,2);
    alphas=[0.05 0.01];

    for pii=1:p;
    Gkl=squeeze(Gvalue(pii,:,:));
    Cklmat=zeros(1,GG);
    for gii=1:GG
    Cklmat(gii)=max(abs(squeeze(Gkl(:,gii))));
    end
    for ii=1:2
     Ckl=quantile(Cklmat,1-alphas(ii)); 
    LowerBd(pii,:,ii)=squeeze(efitBetas(pii,:))-Ckl/sqrt(n)*ones(1,L0);
    UpperBd(pii,:,ii)=squeeze(efitBetas(pii,:))+Ckl/sqrt(n)*ones(1,L0);
    end
    end
    
    str=fullfile(IODir, 'LowerBd.mat');
    save(str,'LowerBd')
    str=fullfile(IODir, 'UpperBd.mat');
    save(str,'UpperBd')
    
  toc
  
  fprintf('+++++++Your job is finished!!+++++++\n');

  guidata(hObject,handles);

