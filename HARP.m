function varargout = HARP(varargin)
% HARP M-file for HARP.fig
%      HARP, by itself, creates a new HARP or raises
%      the existing
%      singleton*.
%
%      H = HARP returns the handle to a new HARP or
%      the handle to
%      the existing singleton*.
%
%      HARP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HARP.M with the given input arguments.
%
%      HARP('Property','Value',...) creates a new HARP or raises the
%      existing singleton*.  Starting from the left, property value pairs
%      are
%      applied to the GUI before importHARP_4_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HARP_OpeningFcn via
%      varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Software:
% 		   Harmonic Phase (HARP)
%
% 	  Image Analysis and Communications Lab (IACL)
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage:
%
% This HARP software is intended for research and educational purposes
% only. You agree to use the HARP software according to the terms and
% conditions outlined in the HARP Disclaimer and End User License
% Agreement (HARP-EULA). In brief, this software has NOT been cleared by
% the Food and Drug Administration, and it may NOT be used for clinical
% or diagnostic purposes.  This software implements only the most basic
% functions of the HARP method and may not be useful for more advanced
% acquisition modes, unusual tag patterns, or advanced tag processing
% goals.  Users are permitted to modify this software for research and
% educational purposes provided that the research and core software is
% properly cited or acknowledged as provided herein and in the
% HARP-EULA.
%
% JHU owns the entire right, title, and interest in and to the "HARP"
% trademark and the patents related to the HARP software. JHU has
% exclusively licensed the HARP software to Diagnosoft, Inc.  A
% commercial, FDA cleared product is available from Diagnosoft,
% Inc. Permission for IACL and Johns Hopkins University to distribute
% this software for research educational purposes has been granted by
% Diagnosoft, Inc. However, the software may not be redistributed and no
% commercial use can be based on this software.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disclaimer:
%
% THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
% APPLICABLE LAW. THE COPYRIGHT HOLDERS AS WELL AS JOHNS HOPKINS
% UNIVERSITY AND DIAGNOSOFT, INC. PROVIDE THE PROGRAM AS IS WITHOUT
% WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED. THE ENTIRE RISK AS
% TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE
% PROGRAM PROVE DEFECTIVE, YOU ASSUME ALL THE COSTS AND RESPONSIBILITIES
% THAT MIGHT RESULT FROM ITS USE.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contributors:
%
% Fangxu Xing, fxing1@jhu.edu
% Sahar Soleimanifard
% Khaled Z. Abd-Elmoniem
% Harsh K. Agarwal
% Yuanming Suo
% Xiaofeng Liu
% Smita Sampath
% Vijay Parthasarathy
% Jonghye Woo
% Aaron Carass
% Nael F. Osman
% Jerry L. Prince, prince@jhu.edu
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acknowledgement:
%
% If you use this HARP software or methods derived from this software to
% process or produce any data for publication, please refer to the
% following publications:
%
% N. F. Osman, W. S. Kerwin, E. R. McVeigh, and J. L. Prince, "Cardiac
% Motion Tracking Using CINE Harmonic Phase (HARP) Magnetic Resonance
% Imaging", Mag. Res. Med., vol. 42, pp. 1048-1060, 1999.
%
% N. F. Osman, E. R. McVeigh, and J. L. Prince, "Imaging Heart Motion
% Using Harmonic Phase MRI", IEEE Trans. on Medical Imaging, vol. 19,
% No. 3, pp. 186-202, March 2000.
%
% X. Liu and J.L. Prince, "Shortest path refinement for motion
% estimation from tagged MR images," IEEE Trans Med Imaging, vol. 29,
% no. 8, pp.1560-1572, March 2010.
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% More information:
%
% IACL HARP page - http://www.iacl.ece.jhu.edu/static/harp/
%
% Wikipedia HARP page - http://en.wikipedia.org/wiki/HARP_algorithm
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright:
%
% Copyright (c) 1999-2013 Johns Hopkins University, Image Analysis
% and Communications Lab (IACL)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @HARP_OpeningFcn, ...
    'gui_OutputFcn',  @HARP_OutputFcn, ...
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


% --- Executes just before HARP is made visible.
function HARP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HARP (see VARARGIN)

% Optimize fft calls on first call.
% This makes subsequent fft calls much faster.
% HARP calls fft many times with the same length.
fftw('planner','patient');

% Choose default command line output for HARP
handles.output = hObject;
set(handles.PNL_SHOW_TRACKS,'SelectionChangeFcn',@PNL_SHOW_TRACKS_SelectionChangeFcn);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HARP wait for user response (see UIRESUME)
% uiwait(handles.FIG_IMPORTHARP);
SliceList = [];
AXIS_IMG_Pressed.IsTrue = false;
AXIS_IMG_Pressed.MouseDownStatus = 0;
setappdata(gcf,'AXIS_IMG_Pressed',AXIS_IMG_Pressed);

setappdata(gcf,'hROI',-1);
setappdata(gcf,'hBRC',-1);
setappdata(gcf,'hTLC',-1);
setappdata(gcf,'hCNTR',-1);
setappdata(gcf,'hCurPt',-1);

setappdata(gcf,'hfEllipse',-1);
setappdata(gcf,'hfCrcArm',-1);
setappdata(gcf,'hfSqrArm',-1);
setappdata(gcf,'hfCNTR',-1);
setappdata(0,'hImportHARPGui',gcf);
setappdata(0,'SliceList',SliceList);

hAXIS_CLR_LIMITS = findobj(gcf,'Tag','AXIS_CLR_LIMITS');
setappdata(0,'hAXIS_CLR_LIMITS',hAXIS_CLR_LIMITS);

hAXIS_IMG = findobj(gcf,'Tag','AXIS_IMG');
setappdata(0,'hAXIS_IMG',hAXIS_IMG);

ColorRange = InitializeImageColorRange;
setappdata(0,'ColorRange',ColorRange);

% If the user doesn't use the exit button, but closes it using the x at
% the top, then we should still go through the Exit button callback.
setappdata(gcf,'fhCloseRequestFunction',@BTN_EXT_Callback);
% Then use that function handle to define the Close request function for
% this figure.
set(gcf,'CloseRequestFcn',getappdata(gcf,'fhCloseRequestFunction'));

% Store the save dataset as a member function because sometimes you
% might need it from visualizeHARP
setappdata(gcf,'fhBTN_DTST_SAV_Callback',@BTN_DTST_SAV_Callback)

% Make GUI resizable
set(handles.FIG_IMPORTHARP,'resize','on');

% Compile HARP refinement c code
mex HarpTrackingRefinement_FastMarching.cpp;


% --- Outputs from this function are returned to the command line.
function varargout = HARP_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function ColorRange = InitializeImageColorRange
% Added by -SS- 01/04/10

% 'Ecc'
ColorRange.Ecc  = [-0.5 0.5];
% 'Err'
ColorRange.Err  = [-0.5 0.5];
% 'Exx'
ColorRange.Exx  = [-0.5 0.5];
% 'Eyy'
ColorRange.Eyy  = [-0.5 0.5];
% 'zEcc'
ColorRange.zEcc = [-0.5 0.5];
% 'zErr'
ColorRange.zErr = [-0.5 0.5];
% 'zExx'
ColorRange.zExx = [-0.5 0.5];
% 'zEyy'
ColorRange.zEyy = [-0.5 0.5];
% '2DEp1'
ColorRange.Ep1_2D = [-0.5 0.5];
% '2DEp2'
ColorRange.Ep2_2D = [-0.5 0.5];
% '2DCircAng1'
ColorRange.CircAng1_2D = [0 90];
% '2DCircAng2'
ColorRange.CircAng2_2D = [0 90];
% 'Ezz'
ColorRange.Ezz = [-0.5 0.5];
% 'Strain_S11_3D'
ColorRange.Strain_S11_3D = [-0.5 0.5];
% 'Strain_S12_3D'
ColorRange.Strain_S12_3D = [-0.5 0.5];
% '3D_Strain_S13'
ColorRange.Strain_S13_3D = [-0.5 0.5];
% '3D_Strain_S22'
ColorRange.Strain_S22_3D = [-0.5 0.5];
% 'Strain_S23_3D'
ColorRange.Strain_S23_3D = [-0.5 0.5];
% 'Strain_S33_3D'
ColorRange.Strain_S33_3D = [-0.5 0.5];
% 'Ep1_3D'
ColorRange.Ep1_3D = [-0.5 0.5];
% 'Ep2_3D'
ColorRange.Ep2_3D = [-0.5 0.5];
% 'Ep3_3D'
ColorRange.Ep3_3D = [-0.5 0.5];
% 'ESig1_3D'
ColorRange.ESig1_3D = [0.5 1.5];
% 'ESig2_3D'
ColorRange.ESig2_3D = [0.5 1.5];
% 'ESig3_3D'
ColorRange.ESig3_3D = [0.5 1.5];
% 'RadAng1_3D'
ColorRange.RadAng1_3D  = [0 90];
% 'CircAng1_3D'
ColorRange.CircAng1_3D = [0 90];
% 'LongAng1_3D'
ColorRange.LongAng1_3D = [0 90];
% 'RadAng2_3D'
ColorRange.RadAng2_3D  = [0 90];
% 'CircAng2_3D'
ColorRange.CircAng2_3D = [0 90];
% 'LongAng2_3D'
ColorRange.LongAng2_3D = [0 90];
% 'RadAng3_3D'
ColorRange.RadAng3_3D  = [0 90];
% 'CircAng3_3D'
ColorRange.CircAng3_3D = [0 90];
% 'LongAng3_3D'
ColorRange.LongAng3_3D = [0 90];
% 'FA_2D'
ColorRange.FA_2D   = [0 0.5];
% 'Indx_2D'
ColorRange.Indx_2D = [0 2];
%  'MD_2D'
ColorRange.MD_2D   = [0 2];
% 'FA'
ColorRange.FA    = [0 0.5];
% 'Indx1'
ColorRange.Indx1 = [0 2];
% 'Indx2'
ColorRange.Indx2 = [0 2];
%  'MD'
ColorRange.MD    = [0 2];


% --- Executes on button press in BTN_FILE_IMPRT.
function BTN_FILE_IMPRT_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_FILE_IMPRT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function executes when pressing "Add Slice" button to read in image
% data.
hImportHARPGui = getappdata(0,'hImportHARPGui');

% Find current path.
hEDT_PATH = findobj(hImportHARPGui,'Tag','EDT_PATH');
pathnameORG = get(hEDT_PATH,'string');
if isempty(pathnameORG)
    pathnameORG = pwd;
    set(hEDT_PATH,'string',pathnameORG);
end

% Detect the format type and scan type.
WhichFormatList = get(handles.LST_IMPORT_FILE_TYPE,'String');
FormatType = get(handles.LST_IMPORT_FILE_TYPE,'Value');
WhichScanList = get(handles.LST_IMPORT_SCAN_TYPE,'String');
ScanType = get(handles.LST_IMPORT_SCAN_TYPE,'Value');
if FormatType == 1
    errordlg('Please select file format!');
    return;
elseif ScanType == 1
    errordlg('Please select scan type!');
    return;
end

% Assign the number of import folders according to scan type.
if strcmpi(WhichFormatList{FormatType},'Philips PAR/REC Files')
    if sum(strcmpi(WhichScanList{ScanType},{'SPAMM-LINES','CSPAMM-LINES','MICSR'}))
        % One file for one (or two) tag direction, all slices, all frames.
        [filename,pathname,dummy] = uigetfile(...
            {'*.par;*.PAR','PHILIPS-header files (*.par)';...
            '*.rec;*.REC','PHILIPS-data files (*.rec)'},...
            'Load all data files',[pathnameORG filesep '*.par'],'MultiSelect','on');
    else
        errordlg('This Scan Type in This File Format is not supported for now!');
        return;
    end
    % if cancelled
    if isequal(filename,0) || isequal(pathname,0)
        return;
    end
elseif strcmpi(WhichFormatList{FormatType},'Siemens DICOM Files')
    % In dicom case image magnitude and phase are stored separately.
    % There should be two files for each image.
    Str =[{'Choose magnitude files for dynamic A1'},...
        {'Choose phase files for dynamic A1'},...
        {'Choose magnitude files for dynamic B1'},...
        {'Choose phase files for dynamic B1'},...
        {'Choose magnitude files for dynamic A2'},...
        {'Choose phase files for dynamic A2'},...
        {'Choose magnitude files for dynamic B2'},...
        {'Choose phase files for dynamic B2'}];
    
    if strcmpi(WhichScanList{ScanType},'SPAMM-LINES')

        % One file per tag direction, per slice, per frame.
        % DC peak is not eliminated.
        if get(findobj(hImportHARPGui,'Tag','CHK_TWODYNAM'),'value') % 2 dynamics
            NumberOfFolders = 2;
        else % only 1 dynamic
            NumberOfFolders = 1;
        end        
        pathname = cell(NumberOfFolders,1);
        filename = [];
        for FIndx = 1:NumberOfFolders
            [Allfiles, pathname{FIndx}, ~] = uigetfile({'*.dcm;*.DCM','Siemens Dicom files (*.dcm)'},...
                Str{FIndx*4-3},[pathnameORG filesep '*.dcm'],'MultiSelect','on');
            % if cancelled
            if isequal(Allfiles,0) || isequal(pathname{FIndx},0)
                return;
            end
            if iscell(Allfiles) % multi files
                filename = [filename;sort(Allfiles)];
            else % one file
                filename = [filename;Allfiles];
            end
        end
        
    elseif sum(strcmpi(WhichScanList{ScanType},'CSPAMM-LINES'))
        
        % Two images (four files) per tag direction, per slice, per frame.
        if get(findobj(hImportHARPGui,'Tag','CHK_TWODYNAM'),'value') % 2 dynamics
            NumberOfFolders = 8;
        else % only 1 dynamic
            NumberOfFolders = 4;
        end
        pathname = cell(NumberOfFolders,1);
        filename = [];
        for FIndx = 1:NumberOfFolders
            [Allfiles, pathname{FIndx}, ~] = uigetfile({'*.dcm;*.DCM','Siemens Dicom files (*.dcm)'},...
                Str{FIndx},[pathnameORG filesep '*.dcm'],'MultiSelect','on');
            % if cancelled
            if isequal(Allfiles,0) || isequal(pathname{FIndx},0)
                return;
            end
            if iscell(Allfiles) % multi files
                filename = [filename;sort(Allfiles)];
            else % one file
                filename = [filename;Allfiles];
            end
        end
        
    elseif sum(strcmpi(WhichScanList{ScanType},'MICSR'))
        
        % Two images (four files) per tag direction, per slice, per frame.
        if get(findobj(hImportHARPGui,'Tag','CHK_TWODYNAM'),'value') % 2 dynamics
            NumberOfFolders = 4;
        else % only 1 dynamic
            NumberOfFolders = 2;
        end
        pathname = cell(NumberOfFolders,1);
        filename = [];
        for FIndx = 1:NumberOfFolders
            [Allfiles, pathname{FIndx}, ~] = uigetfile({'*.dcm;*.DCM','Siemens Dicom files (*.dcm)'},...
                Str{2*FIndx-1},[pathnameORG filesep '*.dcm'],'MultiSelect','on');
            % if cancelled
            if isequal(Allfiles,0) || isequal(pathname{FIndx},0)
                return;
            end
            if iscell(Allfiles) % multi files
                filename = [filename;sort(Allfiles)];
            else % one file
                filename = [filename;Allfiles];
            end
        end
        
    else
        errordlg('This Scan Type in This File Format is not supported for now!');
        return;
    end
    
end

% Get SliceList in memory.
SliceList = getappdata(0,'SliceList');
OrigLength = length(SliceList);
ScanType = WhichScanList{ScanType};

% Get the file format.
if iscell(filename) % multiple files
    ext = filename{1,1}(end-3:end);
else % only 1 file
    ext = filename(end-3:end);
end

[DATA_IN_FILE, DATA_IN_FILE_INFO] = load_data(pathname,filename,ScanType);

celldisp(DATA_IN_FILE_INFO);

% Reading par file as well as DATA_IN_FILE_INFO
% to store Preparation Direction, This Direction as well as
% orientation is stored to verify if MPS coordinate system is
% Left handed or Right handed.

% DATA_IN_FILE_INFO stores the informtation from dicominfo file
% to retrieve plane information later on. It is assumed this
% information is similar for all time frames of one slice.

for NewSliceIndex = 1:length(DATA_IN_FILE.slice)
    
    SliceList{OrigLength+NewSliceIndex} = DATA_IN_FILE.slice{NewSliceIndex};
    SliceList{OrigLength+NewSliceIndex}.NXY = [70 70]; % NXY = Rectangular ROI width in first x and then y directions
    
    SliceList{OrigLength+NewSliceIndex}.ScanType = ScanType;
    SliceList{OrigLength+NewSliceIndex}.TAG_Spacing = 12;
    
    SliceList{OrigLength+NewSliceIndex}.NO_FRAMES = length(SliceList{OrigLength+NewSliceIndex}.dynamic{1}.phase);
    SliceList{OrigLength+NewSliceIndex}.ROI_CNTR  = [ones(1,SliceList{OrigLength+NewSliceIndex}.NO_FRAMES)*SliceList{OrigLength+NewSliceIndex}.NXY(1)/2;...
        ones(1,SliceList{OrigLength+NewSliceIndex}.NO_FRAMES)*SliceList{OrigLength+NewSliceIndex}.NXY(2)/2];
    SliceList{OrigLength+NewSliceIndex}.RefCNTR_Phase1 = 1;
    SliceList{OrigLength+NewSliceIndex}.RefCNTR_Phase2 = SliceList{OrigLength+NewSliceIndex}.NO_FRAMES;
    
    numPhase = length(SliceList{OrigLength+NewSliceIndex}.dynamic{1}.phase);
    SliceList{OrigLength+NewSliceIndex}.CurrentPoint = repmat([5; 5],1,numPhase); % Position of point relative to TLC
    Res =(SliceList{OrigLength+NewSliceIndex}.ReconResolution).*(SliceList{OrigLength+NewSliceIndex}.PxlSpacing);
    TAG = SliceList{OrigLength+NewSliceIndex}.TAG_Spacing;
    
    for DynIndex = 1:length(SliceList{OrigLength+NewSliceIndex}.dynamic)
        
        if strcmpi(ext,'.dcm')
            FindIndx = regexp(DATA_IN_FILE_INFO{NewSliceIndex}.Filename,filesep);
            filename = DATA_IN_FILE_INFO{NewSliceIndex}.Filename(FindIndx(end)+1:end);
        end
        SliceList{OrigLength+NewSliceIndex}.dynamic{DynIndex}.SliceInfo.FileName      = filename;
        SliceList{OrigLength+NewSliceIndex}.dynamic{DynIndex}.SliceInfo.SliceNumber   = NewSliceIndex;
        %-SS- 11/20/09
        %SliceName stores name of the file that contains
        %the slice, SliceNumber stores the number of slice
        %within that file (for files that contain more than one slice)
        SliceList{OrigLength+NewSliceIndex}.dynamic{DynIndex}.TLC            = [30;30];
        SliceList{OrigLength+NewSliceIndex}.dynamic{DynIndex}.TAG_Angle      = pi/2*(1-rem(DynIndex,2));
        SliceList{OrigLength+NewSliceIndex}.dynamic{DynIndex}.Omega          = [0;Res(1)/TAG]+[Res(1)/TAG;-Res(1)/TAG]*(1-rem(DynIndex,2));
        SliceList{OrigLength+NewSliceIndex}.dynamic{DynIndex}.Ratios         = [16; 16];%[1;1]*Res(1)/TAG*4/5;
        SliceList{OrigLength+NewSliceIndex}.dynamic{DynIndex}.Rotation       = 0;
        SliceList{OrigLength+NewSliceIndex}.dynamic{DynIndex}.Decay          = 0.05;
        SliceList{OrigLength+NewSliceIndex}.dynamic{DynIndex}.Filter_Shift   = [0;0];
        [SliceList{OrigLength+NewSliceIndex}.dynamic{DynIndex}.Data,...
            AvailDataList] = GenerateData(SliceList{OrigLength+NewSliceIndex}.dynamic{DynIndex}.phase);
        
        SliceList{OrigLength+NewSliceIndex}.dynamic{DynIndex}                = ...
            rmfield(SliceList{OrigLength+NewSliceIndex}.dynamic{DynIndex},'phase');
        
        SliceList{OrigLength+NewSliceIndex}.AvailResults = AvailDataList;
    end
    
    if strcmpi(ext,'.dcm')
        SliceList{OrigLength+NewSliceIndex}.Planes3D = Slice3d_info_dicom(DATA_IN_FILE_INFO{NewSliceIndex});
    else
        SliceList{OrigLength+NewSliceIndex}.PrpDir = DATA_IN_FILE_INFO.PrpDir;
        % from this PrpDir we can find if coordinate system is right-handed or left-handed
        SliceList{OrigLength+NewSliceIndex}.Planes3D = Slice3d_info(SliceList{OrigLength+NewSliceIndex});
    end
    SliceList{OrigLength+NewSliceIndex} = SaveOrgProcessInfo(SliceList{OrigLength+NewSliceIndex});
end

setappdata(0,'SliceList',SliceList);
UpdatePreviewROIGUI;
UpdateLST_SLICES;


function myslice = SaveOrgProcessInfo(myslice)

% This function is first added to the code on 12-14-09 by -SS-
% Stores the orginial info of "processinfo" data such as original
% sample numbers, filter status, total frames,...

% flag off, shows that never processed.
if ~isfield(myslice,'ProcessInfo')
    myslice.ProcessInfo = [];
end

% upsampling
if ~isfield(myslice.ProcessInfo,'UpsampleStatus')
    myslice.ProcessInfo.UpsampleStatus.flag = 'off';
    myslice.ProcessInfo.UpsampleStatus.OrgReconResolution = myslice.ReconResolution;
    myslice.ProcessInfo.UpsampleStatus.OrgPxlSpacing = myslice.PxlSpacing;
end

% Filtering
if ~isfield(myslice.ProcessInfo,'FilterStatus')
    myslice.ProcessInfo.FilterStatus.flag = 'dirty';
end

% 2D Strain Analysis
if ~isfield(myslice.ProcessInfo,'Strain2DStatus')
    myslice.ProcessInfo.Strain2DStatus.flag = 'dirty';
end

% Tracking
if ~isfield(myslice.ProcessInfo,'TrackStatus')
    myslice.ProcessInfo.TrackStatus.flag        = 'dirty';
    myslice.ProcessInfo.TrackStatus.Harp.flag   = 'dirty';
    myslice.ProcessInfo.TrackStatus.Refine.flag = 'dirty';
end


function  [mydata,AvailDataList] = GenerateData(myphase)

numFrames = length(myphase);

% First 4 types of data that are available to display
AvailDataList = {'Magnitude';'Spectrum';'Checkerboard';'ROI Magnitude'};

if(isfield(myphase{1},'recon'))
    for frame_ind=1:numFrames,
        mydata(frame_ind,:,:)=myphase{frame_ind}.recon;
    end
elseif(isfield(myphase{1},'magnitude'))
    mydata = zeros([numFrames size(myphase{1}.magnitude)]);
    
    if(isfield(myphase{1},'angle'))
        for frame_ind=1:numFrames,
            mydata(frame_ind,:,:)=myphase{frame_ind}.magnitude.*exp(sqrt(-1)*myphase{frame_ind}.angle);
        end
        AvailDataList{5} = 'Phase';
        AvailDataList{6} = 'Real';
        AvailDataList{7} = 'Imag';
    else
        for frame_ind=1:numFrames,
            mydata(frame_ind,:,:)=myphase{frame_ind}.magnitude;
        end
    end
    
elseif isfield(myphase{1},'real')
    mydata = zeros([numFrames size(myphase{1}.real)]);
    if(isfield(myphase{1},'imag'))
        for frame_ind=1:numFrames,
            mydata(frame_ind,:,:)=myphase{frame_ind}.real+sqrt(-1)*myphase{frame_ind}.imag;
        end
        AvailDataList{5} = 'Phase';
        AvailDataList{6} = 'Real';
        AvailDataList{7} = 'Imag';
    else
        for frame_ind=1:numFrames,
            mydata(frame_ind,:,:)=myphase{frame_ind}.real;
        end
    end
else
    mydata = zeros([numFrames size(myphase{1}.imag)]);
    for frame_ind=1:numFrames,
        mydata(frame_ind,:,:)=myphase{frame_ind}.imag;
    end
    AvailDataList = 'Imag';
end


% --- Executes on slider movement.
function SLDR_PHS_Callback(hObject, eventdata, handles)
% hObject    handle to SLDR_PHS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

UpdatePreviewImage('SLDR_PHS');


% --- Executes on slider movement.
function SLDR_IMG_MASK_Callback(hObject, eventdata, handles)
% hObject    handle to SLDR_IMG_MASK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

UpdateImageShow;


% --- Executes on slider movement.
function SLDR_IMG_BRIGHT_Callback(hObject, eventdata, handles)
% hObject    handle to SLDR_IMG_BRIGHT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of
%        slider
UpdateImageShow;


% --- Executes on selection change in LST_SLICES.
function LST_SLICES_Callback(hObject, eventdata, handles)
% hObject    handle to LST_SLICES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns LST_SLICES contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LST_SLICES
UpdateLST_DYNAMICS;
UpdatePreviewROIGUI;


function UpdateLST_DYNAMICS
SliceList = getappdata(0,'SliceList');
hImportHARPGui = getappdata(0,'hImportHARPGui');
if(isempty(hImportHARPGui))
    return;
end

hLST_SLICES  = findobj(hImportHARPGui,'Tag','LST_SLICES');
if(isempty(hLST_SLICES))
    return;end;
CurrentSlice = get(hLST_SLICES,'Value');

DynListLength=0;

if(~isempty(SliceList))
    DynListLength = length(SliceList{CurrentSlice}.dynamic);
end;

StringList=cell(1,DynListLength);
for DynIndex = 1:DynListLength ,
    StringList{DynIndex}=['Dynamic' num2str(DynIndex)];
end

hLST_DYNAMICS  = findobj(hImportHARPGui,'Tag','LST_DYNAMICS');
CurrentDynamic = get(hLST_DYNAMICS,'Value');
if CurrentDynamic > length(StringList)
    set(hLST_DYNAMICS,'Value',1);
end
set(hLST_DYNAMICS,'String',StringList);

UpdatePropertiesGUI;
UpdatePreviewGUI;


function UpdatePreviewGUI

SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

hImportHARPGui = getappdata(0,'hImportHARPGui');

hLST_SLICES  = findobj(hImportHARPGui,'Tag','LST_SLICES');
CurrentSlice = get(hLST_SLICES,'Value');

hLST_DYNAMICS = findobj(hImportHARPGui,'Tag','LST_DYNAMICS');
CurrentDynamic = get(hLST_DYNAMICS,'Value');

hLST_PHS = findobj(hImportHARPGui,'Tag','LST_PHS');
hLST_PRVWTYPE = findobj(hImportHARPGui,'Tag','LST_PRVWTYPE');
hSLDR_PHS = findobj(hImportHARPGui,'Tag','SLDR_PHS');

PhaseListIndex = 1:size(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Data,1);
set(hLST_PHS,'String',num2str(PhaseListIndex'),'Value',1);

if max(PhaseListIndex)>1
    set(hSLDR_PHS,'Min',1,'Max',max(PhaseListIndex));
    set(hSLDR_PHS,'SliderStep',[1/max(PhaseListIndex) 2/max(PhaseListIndex)]);
    set(hSLDR_PHS,'Value',1);
end

CurrentValue = get(hLST_PRVWTYPE,'Value');
MaxNewValue = length(SliceList{CurrentSlice}.AvailResults);
set(hLST_PRVWTYPE,'String',SliceList{CurrentSlice}.AvailResults);
if(CurrentValue>MaxNewValue) || (CurrentValue<1)
    set(hLST_PRVWTYPE,'Value',1);
end

UpdatePreviewImage('LST_PHS');


function UpdatePreviewROIGUI

SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

hImportHARPGui = getappdata(0,'hImportHARPGui');

hAXIS_IMG = getappdata(0,'hAXIS_IMG');
set(hImportHARPGui,'CurrentAxes',hAXIS_IMG);

hLST_SLICES     = findobj(hImportHARPGui,'Tag','LST_SLICES');
CurrentSlice    = get(hLST_SLICES,'Value');

hLST_DYNAMICS   = findobj(hImportHARPGui,'Tag','LST_DYNAMICS');
CurrentDynamic  = get(hLST_DYNAMICS,'Value');

hLST_PHS        = findobj(hImportHARPGui,'Tag','LST_PHS');
CurrentPhase    = get(hLST_PHS,'Value');

CKBUpSample = 1;

NXYT = SliceList{CurrentSlice}.ReconResolution;
NXY = SliceList{CurrentSlice}.NXY;
TLC = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TLC;
ROI_CNTR = SliceList{CurrentSlice}.ROI_CNTR;
CurPt = SliceList{CurrentSlice}.CurrentPoint(:,CurrentPhase)+TLC;

if(strcmp(SliceList{CurrentSlice}.ScanType,'d-zHARP'))
    OMEGA = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Filter_Shift;
else
    OMEGA = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Omega+SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Filter_Shift;
end

RATIOS = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Ratios;
ROTATION = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Rotation * pi / 180;

hLST_PRVWTYPE = findobj(hImportHARPGui,'Tag','LST_PRVWTYPE' );
PreviewType = get(hLST_PRVWTYPE,'Value');
AvailShowList = get(hLST_PRVWTYPE,'String');

hROI            = getappdata(hImportHARPGui,'hROI');
hTLC            = getappdata(hImportHARPGui,'hTLC');
hBRC            = getappdata(hImportHARPGui,'hBRC');
hCNTR           = getappdata(hImportHARPGui,'hCNTR');
hCurPt          = getappdata(hImportHARPGui,'hCurPt');

AXIS_IMG_Pressed= getappdata(hImportHARPGui,'AXIS_IMG_Pressed');

hfEllipse       = getappdata(hImportHARPGui,'hfEllipse'       );
hfCrcArm        = getappdata(hImportHARPGui,'hfCrcArm'        );
hfSqrArm        = getappdata(hImportHARPGui,'hfSqrArm'        );
hfCNTR          = getappdata(hImportHARPGui,'hfCNTR'          );

if(ishandle(hROI)),      delete(hROI) ;     end;
if(ishandle(hTLC)),      delete(hTLC) ;     end;
if(ishandle(hBRC)),      delete(hBRC) ;     end;
if(ishandle(hCNTR)),     delete(hCNTR);     end;
if(ishandle(hCurPt)),    delete(hCurPt);    end;

if(ishandle(hfEllipse)), delete(hfEllipse); end;
if(ishandle(hfCrcArm)),  delete(hfCrcArm) ; end;
if(ishandle(hfSqrArm)),  delete(hfSqrArm) ; end;
if(ishandle(hfCNTR)),    delete(hfCNTR)   ; end;

if strcmpi(AvailShowList{PreviewType},'Spectrum') %--> Draw Ellipse, ...
    nx    = NXYT(1);
    ny    = NXYT(2);
    ang   = pi*2/36;
    rot   = [cos(ROTATION) sin(ROTATION); -sin(ROTATION) cos(ROTATION)];
    
    for n = 0:35,
        p(1:2,n+1)=rot*[RATIOS(1) 0; 0 RATIOS(2)]*[cos(ang*n); sin(ang*n)]+OMEGA;
    end;
    F1_A      = rot*[RATIOS(1) 0; 0 RATIOS(2)]*[1;0]+OMEGA;
    F2_A      = rot*[RATIOS(1) 0; 0 RATIOS(2)]*[0;1]+OMEGA;
    
    hold on;
    hfEllipse = plot([p(2,1:36) p(2,1)]+ny/2+1,[p(1,1:36) p(1,1)]+nx/2+1,'y');
    hfCNTR    = plot(OMEGA(2)+ny/2+1,OMEGA(1)+nx/2+1, '+y');
    hfSqrArm  = plot(F1_A(2)+ny/2+1,F1_A(1)+nx/2+1  , 'sy');
    hfCrcArm  = plot(F2_A(2)+ny/2+1,F2_A(1)+nx/2+1  , 'oy');
    hold off;
    
    setappdata(hImportHARPGui,'hfEllipse',hfEllipse);
    setappdata(hImportHARPGui,'hfCNTR'   ,hfCNTR   );
    setappdata(hImportHARPGui,'hfSqrArm' ,hfSqrArm );
    setappdata(hImportHARPGui,'hfCrcArm' ,hfCrcArm );
    
else
    if CurPt(1)-TLC(1)<=1 || CurPt(1)-TLC(1)>NXY(1)*CKBUpSample || ...
            CurPt(2)-TLC(2)<=1 || CurPt(2)-TLC(2)>NXY(2)*CKBUpSample
        CurPt(1) = TLC(1) + round(NXY(1)*CKBUpSample/4);
        CurPt(2) = TLC(2) + round(NXY(2)*CKBUpSample/4);
        SliceList{CurrentSlice}.CurrentPoint(:, CurrentPhase) = CurPt-TLC;
        setappdata(0,'SliceList',SliceList);
    end
    
    if sum(strcmpi(AvailShowList{PreviewType},{'Magnitude','Phase','Checkerboard','Real','Imag'})) %Special--> Draw ROI for all others
        p(1:2,1) = TLC;
        p(1:2,2) = TLC+[NXY(1)-1;0];
        p(1:2,3) = TLC+[NXY(1)-1;NXY(2)-1];
        p(1:2,4) = TLC+[0;NXY(2)-1];
        
        hold on;
        
        hROI=plot([p(2,1:4) p(2,1)],[p(1,1:4) p(1,1)],'y');
        
        R1_A=p(1:2,1); % top left corner, motion
        R2_A=p(1:2,3); % bottom right corner, resizing
        R3_A=R1_A+ROI_CNTR(:,CurrentPhase); %Center
        
        if ((AXIS_IMG_Pressed.MouseDownStatus == 1)||(AXIS_IMG_Pressed.MouseDownStatus ==11))
            hTLC=plot(R1_A(2),R1_A(1),'sg');%Top Left Corner
        else
            hTLC=plot(R1_A(2),R1_A(1),'sy');
        end;
        if (AXIS_IMG_Pressed.MouseDownStatus == 2)
            hBRC=plot(R2_A(2),R2_A(1),'og');%Bottom Right Corner
        else
            hBRC=plot(R2_A(2),R2_A(1),'oy');
        end;
        if ((AXIS_IMG_Pressed.MouseDownStatus == 3)||(AXIS_IMG_Pressed.MouseDownStatus ==13))
            hCNTR=plot(R3_A(2),R3_A(1),'xg'); %center of ROI
        else
            hCNTR=plot(R3_A(2),R3_A(1),'xy');
        end;
    end
    hold on
    % Draw current point
    if sum(strcmpi(AvailShowList{PreviewType},{'Magnitude','Phase','Checkerboard','Real','Imag'}))
        hCurPt = plot(CurPt(2), CurPt(1), 'yo', 'markersize', 4, 'linewidth', 1.5);
    else
        hCurPt = plot(CurPt(2)-TLC(2), CurPt(1)-TLC(1), 'yo', 'markersize', 4, 'linewidth', 1.5);
    end
    
    hold off;
    setappdata(hImportHARPGui,'hROI'  , hROI  );
    setappdata(hImportHARPGui,'hBRC'  , hBRC  );
    setappdata(hImportHARPGui,'hTLC'  , hTLC  );
    setappdata(hImportHARPGui,'hCNTR' , hCNTR );
    setappdata(hImportHARPGui,'hCurPt', hCurPt);
end


function UpdateLST_SLICES

SliceList = getappdata(0,'SliceList');
hImportHARPGui = getappdata(0,'hImportHARPGui');

SliceListLength = length(SliceList);

StringList = cell(SliceListLength,1);
for SliceIndex=1:SliceListLength
    SlcNo = SliceList{SliceIndex}.dynamic{1}.SliceInfo.SliceNumber;
    StringList{SliceIndex} = ['Slice' num2str(SlcNo) '-' SliceList{SliceIndex}.ScanType];
end

hLST_SLICES=findobj(hImportHARPGui,'Tag','LST_SLICES');
set(hLST_SLICES,'String',StringList);
set(hLST_SLICES,'Value',1);

UpdateLST_DYNAMICS;


% --- Executes on selection change in LST_DYNAMICS.
function LST_DYNAMICS_Callback(hObject, eventdata, handles)
% hObject    handle to LST_DYNAMICS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns LST_DYNAMICS contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LST_DYNAMICS

UpdatePropertiesGUI;

if (getappdata(0,'hAXIS_IMG') == hObject)
    UpdatePreviewImage('ImageAxes');
else
    UpdatePreviewGUI;
end

hImportHARPGui  =getappdata(0,'hImportHARPGui');

if(findobj(hImportHARPGui,'Tag','LST_DYNAMICS') == hObject)
    %BTN_ANIM_Callback(0, 0, 0);%-SS-
end


function UpdatePropertiesGUI

SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

hImportHARPGui = getappdata(0,'hImportHARPGui');
ColorRange = getappdata(0,'ColorRange');

hLST_SLICES     = findobj(hImportHARPGui,'Tag','LST_SLICES');
CurrentSlice    = get(hLST_SLICES,'Value');
hLST_DYNAMICS   = findobj(hImportHARPGui,'Tag','LST_DYNAMICS');
CurrentDynamic  = get(hLST_DYNAMICS,'Value');

hLST_PRVWTYPE   = findobj(hImportHARPGui,'Tag','LST_PRVWTYPE');
ReqPreviewType  = get(hLST_PRVWTYPE,'Value');
AvailShowList   = get(hLST_PRVWTYPE,'String');

hEDT_FOV             = findobj(hImportHARPGui,'Tag','EDT_FOV');
hEDT_RES             = findobj(hImportHARPGui,'Tag','EDT_RES');
hEDT_EDT_TAG_SPACING = findobj(hImportHARPGui,'Tag','EDT_TAG_SPACING');
hEDT_EDT_TAG_ANGLE   = findobj(hImportHARPGui,'Tag','EDT_TAG_ANGLE');
hEDT_SLC_LOCATION    = findobj(hImportHARPGui,'Tag','EDT_SLC_LOCATION');
hEDT_SLC_NUMBER      = findobj(hImportHARPGui,'Tag','EDT_SLC_NUMBER');
hEDT_DATSET_NAME     = findobj(hImportHARPGui,'Tag','EDT_DATSET_NAME');
hEDT_TLC             = findobj(hImportHARPGui,'Tag','EDT_TLC');
hEDT_NXY             = findobj(hImportHARPGui,'Tag','EDT_NXY');
hEDT_OMEGA           = findobj(hImportHARPGui,'Tag','EDT_OMEGA');
hEDT_FLTR_SHFT       = findobj(hImportHARPGui,'Tag','EDT_FLTR_SHFT');
hEDT_RATIOS          = findobj(hImportHARPGui,'Tag','EDT_RATIOS');
hEDT_DECAY           = findobj(hImportHARPGui,'Tag','EDT_DECAY');
hEDT_ROTATION        = findobj(hImportHARPGui,'Tag','EDT_ROTATION');
hEDT_CLR_RANGE       = findobj(hImportHARPGui,'Tag','EDT_CLR_RANGE');
hSLDR_IMG_MASK       = findobj(hImportHARPGui,'Tag','SLDR_IMG_MASK');
hSLDR_IMG_BRIGHT     = findobj(hImportHARPGui,'Tag','SLDR_IMG_BRIGHT');
hBTN_TRACK_PT = findobj(hImportHARPGui,'Tag','BTN_TRACK_PT');

if ~strcmpi(AvailShowList{ReqPreviewType},'Spectrum')
    set(hSLDR_IMG_BRIGHT,'enable','on');
else
    set(hSLDR_IMG_BRIGHT,'enable','off');
end

HideMask = sum(strcmpi(AvailShowList{ReqPreviewType},[{'Magnitude'},{'Spectrum'},{'Checkerboard'},{'Phase'},{'Real'},{'Imag'}]));
if isfield(SliceList{CurrentSlice},'HarpMagnitude') && ~HideMask
    set(hBTN_TRACK_PT,'enable','on');
    set(hSLDR_IMG_MASK,'enable','on');
else
    set(hBTN_TRACK_PT,'enable','off');
    set(hSLDR_IMG_MASK,'enable','off');
end

set(hEDT_EDT_TAG_SPACING,'string',num2str(SliceList{CurrentSlice}.TAG_Spacing));
set(hEDT_EDT_TAG_ANGLE  ,'string',num2str(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TAG_Angle/pi*180));
set(hEDT_FOV            ,'string',num2str(int32([(SliceList{CurrentSlice}.PxlSpacing).*(SliceList{CurrentSlice}.ReconResolution) ...
    SliceList{CurrentSlice}.Thickness])));
set(hEDT_RES            ,'string',num2str(int32(SliceList{CurrentSlice}.ReconResolution)));

if iscell(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.SliceInfo.FileName)
    FileName = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.SliceInfo.FileName{1};
else
    FileName = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.SliceInfo.FileName;
end
set(hEDT_SLC_LOCATION,'string',FileName);
set(hEDT_SLC_NUMBER,'string',SliceList{CurrentSlice}.dynamic{CurrentDynamic}.SliceInfo.SliceNumber);

if isfield(SliceList{CurrentSlice},'DatasetInfo')
    set(hEDT_DATSET_NAME,'string',SliceList{CurrentSlice}.DatasetInfo.FileName);
else
    set(hEDT_DATSET_NAME,'string',' ');
end
set(hEDT_TLC            ,'string',num2str(int32(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TLC)'));
set(hEDT_NXY            ,'string',num2str(int32(SliceList{CurrentSlice}.NXY)));
set(hEDT_OMEGA          ,'string',num2str((SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Omega)'));
set(hEDT_FLTR_SHFT      ,'string',num2str((SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Filter_Shift)'));
set(hEDT_RATIOS         ,'string',num2str((SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Ratios)'));
set(hEDT_DECAY          ,'string',num2str(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Decay));
set(hEDT_ROTATION       ,'string',num2str(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Rotation));

if isfield(ColorRange,AvailShowList{ReqPreviewType})
    ImRange = eval(['ColorRange.' AvailShowList{ReqPreviewType}]);
    set(findobj(hImportHARPGui,'Tag','AXIS_CLR_LIMITS'),'Visible','on');
    set(hEDT_CLR_RANGE,'Visible','on');
    set(hEDT_CLR_RANGE,'string',num2str(ImRange));
    
else
    set(hEDT_CLR_RANGE,'Visible','off');
    set(findobj(hImportHARPGui,'Tag','AXIS_CLR_LIMITS'),'Visible','off');
end

hRAD_SHOW_TRK_HARP = findobj(hImportHARPGui,'Tag','RAD_SHOW_TRK_HARP');
hRAD_SHOW_TRK_REFINE = findobj(hImportHARPGui,'Tag','RAD_SHOW_TRK_REFINE');
if strcmpi(SliceList{CurrentSlice}.ProcessInfo.TrackStatus.Harp.flag,'clean')
    set(hRAD_SHOW_TRK_HARP,'Enable','on');
else
    set(hRAD_SHOW_TRK_HARP,'Enable','off');
end
if strcmp(SliceList{CurrentSlice}.ProcessInfo.TrackStatus.Refine.flag,'clean')
    set(hRAD_SHOW_TRK_REFINE,'Enable','on');
else
    set(hRAD_SHOW_TRK_REFINE,'Enable','off');
end


% --- Executes on button press in BTN_SLC_DLT.
function BTN_SLC_DLT_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_SLC_DLT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

hImportHARPGui  = getappdata(0,'hImportHARPGui');
hLST_SLICES = findobj(hImportHARPGui,'Tag','LST_SLICES');
CurrentSlice = get(hLST_SLICES,'Value');

TempInd = 1:length(SliceList);
SliceList = {SliceList{TempInd(TempInd~=CurrentSlice)}};
setappdata(0,'SliceList',SliceList);

UpdateLST_SLICES;

if isempty(SliceList)
    axes(handles.AXIS_IMG)
    cla;
end


% --- Executes on button press in BTN_SET_PATH.
function BTN_SET_PATH_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_SET_PATH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Sets current data path to any folder.

% Get current data folder.
hImportHARPGui  =getappdata(0,'hImportHARPGui');
hEDT_PATH = findobj(hImportHARPGui,'Tag','EDT_PATH');
pathnameORG = get(hEDT_PATH,'string');
if isempty(pathnameORG)
    pathnameORG = pwd;
end
PathName = uigetdir(pathnameORG,'Select the data folder...');

% Set current data folder.
if PathName ~= 0
    set(hEDT_PATH,'string',PathName);
end


% --- Executes on button press in BTN_DTST_SAV.
function BTN_DTST_SAV_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_DTST_SAV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

hImportHARPGui = getappdata(0,'hImportHARPGui');
hEDT_PATH = findobj(hImportHARPGui,'Tag','EDT_PATH');
pathname = get(hEDT_PATH,'string');
filename = 'new';

if isfield(SliceList{1},'DatasetInfo')
    % check if there is a date extension
    % date format: 'yyyy-mm-dd'
    filename = SliceList{1}.DatasetInfo.FileName;
    [match start]=regexpi(filename,'-(\d{4})-(\d{2})-(\d{2})','match','start');
    if ~isempty(match)
        filename = filename(1:start-1);
    end
end

if isempty(pathname)
    pathname = pwd; %current directory
end

% The program gives suggestions for dataset names based on the information
% saved in first slice.

[filename, pathname, dummy]    = uiputfile({'*.mat'...
    ,'Matlab files (*.mat)'},'Please select a folder...',...
    [pathname filesep filename '-' datestr(clock,'yyyy-mm-dd')]);
if isequal(filename,0) || isequal(pathname,0)
    % User pressed cancel
    return;
end

filename = filename(1:end-4); %excluding .mat extension
pathname = pathname(1:end-1); % excluding last filesep to be consistent

% This is taking a bit long, so add waitbar
hWaitbar = waitbar(0,'Saving dataset. Please wait...');
waitbar(0.25,hWaitbar);

for ind = 1:length(SliceList)
    SliceList{ind}.DatasetInfo.FileName = filename;
    SliceList{ind}.DatasetInfo.PathName = pathname;
    SliceList{ind}.DatasetInfo.Time = datestr(clock,0);
end

waitbar(0.5,hWaitbar);

save([pathname filesep filename '.mat'],'SliceList');
setappdata(0,'SliceList',SliceList);
UpdatePropertiesGUI;

waitbar(1,hWaitbar);
% Close the waitbar;
close(hWaitbar);


% --- Executes on button press in BTN_DTST_LOD.
function BTN_DTST_LOD_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_DTST_LOD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SliceList      = getappdata(0,'SliceList');
hImportHARPGui = getappdata(0,'hImportHARPGui');
hEDT_PATH = findobj(hImportHARPGui,'Tag','EDT_PATH');
CurrentPath = get(hEDT_PATH,'string');
if isempty(CurrentPath)
    CurrentPath = pwd;
    set(hEDT_PATH,'string',CurrentPath);
end

[filename, pathname, dummy] = uigetfile(...
    {'*.mat','Matlab files (*.mat)'},'Choose a file',...
    [CurrentPath filesep '*.mat']);

if isequal(filename,0) || isequal(pathname,0)
    % User pressed cancel
    return;
end

hWaitbar = waitbar(0,'Loading dataset. Please wait...');
waitbar(0.5,hWaitbar);

dataset = load([pathname filename],'SliceList');
waitbar(0.75,hWaitbar);

% This is useful when loading older datasets
for slice = 1:length(dataset.SliceList)
    dataset.SliceList{slice} = SaveOrgProcessInfo(dataset.SliceList{slice});
    
    if size(dataset.SliceList{slice}.NXY) == 1
        dataset.SliceList{slice}.NXY(2) = dataset.SliceList{slice}.NXY(1);
    end
end

if(isempty(SliceList))
    SliceList = dataset.SliceList;
else
    OrigLength = length(SliceList);
    for NewSliceIndex = 1:length(dataset.SliceList)
        SliceList{OrigLength+NewSliceIndex} = dataset.SliceList{NewSliceIndex};
    end
end

waitbar(1,hWaitbar);

setappdata(0,'SliceList',SliceList);

close(hWaitbar);
UpdateLST_SLICES;
UpdatePreviewROIGUI;


% --- Executes on button press in BTN_EXT.
function BTN_EXT_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_EXT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FIG_IMPORTHARP_CloseRequestFcn(getappdata(0,'hImportHARPGui'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes when user attempts to close FIG_IMPORTHARP.
function FIG_IMPORTHARP_CloseRequestFcn(hObject)
% hObject    handle to FIG_IMPORTHARP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isappdata(0,'SliceList')
    rmappdata(0,'SliceList');
end

if isappdata(0,'hThumbnailAxes')
    rmappdata(0,'hThumbnailAxes');
end

if isappdata(0,'hCurrentImageAxes')
    rmappdata(0,'hCurrentImageAxes');
end

if isappdata(0,'hPlotAxes')
    rmappdata(0,'hPlotAxes');
end

% If Visualize module exists, then kill it, and remove the handle from
% the root. This is good house-keeping.
if isappdata(0,'hVisualizeHARPGui')
    hVisualizeHARPGui = getappdata(0,'hVisualizeHARPGui');
    if ishandle(hVisualizeHARPGui)
        close(hVisualizeHARPGui);
    end
    rmappdata(0,'hVisualizeHARPGui');
    % Note that there is separate close request function for each
    % separate GUI, which will then clean up the root handles as well.
end

if(isappdata(0,'hImportHARPGui'))
    rmappdata(0,'hImportHARPGui');
end
% Finally destroy the current figure also.
delete(hObject);


function UpdateSliceAndDynamicProperties
SliceList       =getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

hImportHARPGui  =getappdata(0,'hImportHARPGui');

hLST_SLICES     = findobj(hImportHARPGui,'Tag','LST_SLICES');
if(isempty(hLST_SLICES))
    return;
end;

CurrentSlice    = get(hLST_SLICES,'Value');
hLST_DYNAMICS   = findobj(hImportHARPGui,'Tag','LST_DYNAMICS');
CurrentDynamic  = get(hLST_DYNAMICS,'Value');
%hLST_PHS        = findobj(hImportHARPGui,'Tag','LST_PHS');
%CurrentPhase    = get(hLST_PHS,'Value');

hEDT_EDT_TAG_SPACING = findobj(hImportHARPGui,'Tag','EDT_TAG_SPACING');
hEDT_EDT_TAG_ANGLE   = findobj(hImportHARPGui,'Tag','EDT_TAG_ANGLE');
hLST_SCAN_TYPE       = findobj(hImportHARPGui,'Tag','LST_SCAN_TYPE');
hEDT_TLC             = findobj(hImportHARPGui,'Tag','EDT_TLC');
hEDT_NXY             = findobj(hImportHARPGui,'Tag','EDT_NXY');
hEDT_OMEGA           = findobj(hImportHARPGui,'Tag','EDT_OMEGA');
hEDT_RATIOS          = findobj(hImportHARPGui,'Tag','EDT_RATIOS');
hEDT_FLTER_SHFT      = findobj(hImportHARPGui,'Tag','EDT_FLTR_SHFT');
hEDT_DECAY           = findobj(hImportHARPGui,'Tag','EDT_DECAY');
hEDT_ROTATION        = findobj(hImportHARPGui,'Tag','EDT_ROTATION');


x = str2double(get(hEDT_EDT_TAG_ANGLE,'string'));
if(~isnan(x) && (length(x)==1))
    if (x~= SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TAG_Angle)
        
        SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TAG_Angle = 2*pi/360*x;
        FOV = SliceList{CurrentSlice}.ReconResolution .* SliceList{CurrentSlice}.PxlSpacing ;
        SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Omega = ...
            FOV(1)/SliceList{CurrentSlice}.TAG_Spacing * ([sin(2*pi/360*x); cos(2*pi/360*x)]);
        SliceList{CurrentSlice}.ProcessInfo.FilterStatus.flag = 'dirty';
    end
else
    disp('Tag angle must be a 1x1 number. Input ignored');
    set(hEDT_EDT_TAG_ANGLE,'string', num2str(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TAG_Angle));
end;

[x,ok] = str2num(get(hEDT_EDT_TAG_SPACING,'string'));
if(ok && (length(x)==1))
    if (x ~= SliceList{CurrentSlice}.TAG_Spacing)
        SliceList{CurrentSlice}.TAG_Spacing = x;
        tag_angle = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TAG_Angle;
        FOV = SliceList{CurrentSlice}.ReconResolution .* SliceList{CurrentSlice}.PxlSpacing ;
        SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Omega = FOV(1)/x * ([sin(tag_angle); cos(tag_angle)]);
        set(hEDT_OMEGA,'string',num2str((SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Omega)'));
        SliceList{CurrentSlice}.ProcessInfo.FilterStatus.flag = 'dirty';
    end
else
    disp('tag spacing must be a 1x1 number. Input ignored');
    set(hEDT_EDT_TAG_SPACING,'string', num2str(SliceList{CurrentSlice}.TAG_Spacing));
end;

[x,ok] = str2num(get(hEDT_TLC,'string'));
if(ok && (length(x)==2))
    if ~isequal(x,SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TLC')
        SliceList{CurrentSlice}.dynamic{1}.TLC = x';
        SliceList{CurrentSlice}.dynamic{2}.TLC = x';
        SliceList{CurrentSlice}.ProcessInfo.FilterStatus.flag = 'dirty';
    end
else
    disp('TLC must be a 1x2 number. Input ignored');
    set(hEDT_TLC,'string',num2str((SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TLC)'));
end;

[x,ok] = str2num(get(hEDT_NXY,'string'));
if(ok && (length(x)==2))
    if ~isequal(x,SliceList{CurrentSlice}.NXY)
        SliceList{CurrentSlice}.NXY = x;
        SliceList{CurrentSlice}.ProcessInfo.FilterStatus.flag = 'dirty';
    end
else
    disp('NXY must be a 1x2 number. Input ignored');
    set(hEDT_NXY,'string',num2str(SliceList{CurrentSlice}.NXY));
end;

[dummy,ok] = str2num(get(hEDT_OMEGA,'string'));
if(ok)
    %SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Omega = x';
else
    disp('OMEGA must be a number. Input ignored');
end;

[x,ok] = str2num(get(hEDT_RATIOS,'string'));
if(ok && (length(x)==2))
    if (x ~= SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Ratios')
        SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Ratios = x';
        SliceList{CurrentSlice}.ProcessInfo.FilterStatus.flag = 'dirty';
    end
else
    disp('Filter ratios must be a 1x2 number. Input ignored');
    set(hEDT_RATIOS,'string',num2str((SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Ratios)'));
end;

x = str2double(get(hEDT_DECAY,'string'));
if(~isnan(x) && (length(x)==1))
    if (x ~= SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Decay)
        SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Decay = x;
        SliceList{CurrentSlice}.ProcessInfo.FilterStatus.flag = 'dirty';
    end
else
    disp('Filter decay must be a 1x1 number. Input ignored');
    set(hEDT_DECAY,'string',num2str(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Decay));
end;

[x,ok] = str2num(get(hEDT_ROTATION,'string'));
if(ok && (length(x)==1))
    if (x ~= SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Rotation)
        SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Rotation = x;
        SliceList{CurrentSlice}.ProcessInfo.FilterStatus.flag = 'dirty';
    end
else
    disp('rotation must be a 1x1 number. Input ignored');
    set(hEDT_ROTATION,'string',num2str(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Rotation));
end;

[x,ok] = str2num(get(hEDT_FLTER_SHFT,'string'));
if(ok && (length(x)==2))
    if (x ~= SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Filter_Shift')
        SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Filter_Shift = x';
        SliceList{CurrentSlice}.ProcessInfo.FilterStatus.flag = 'dirty';
    end
else
    disp('Filter shift must be a 1x2 number. Input ignored');
    set(hEDT_FLTER_SHFT,'string',num2str((SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Filter_Shift)'));
end;

ScanTypes = get(hLST_SCAN_TYPE,'String');
NewType   = get(hLST_SCAN_TYPE,'Value');
if NewType ~= 1
    if ~strcmp(SliceList{CurrentSlice}.ScanType,ScanTypes{NewType})
        SliceList{CurrentSlice}.ScanType = ScanTypes{NewType};
        SliceList{CurrentSlice}.ProcessInfo.FilterStatus.flag = 'dirty';
    end
end


setappdata(0,'SliceList',SliceList);
UpdatePreviewROIGUI;


% --- Executes on selection change in LST_IMPORT_SCAN_TYPE.
function LST_IMPORT_SCAN_TYPE_Callback(hObject, eventdata, handles)
% hObject    handle to LST_IMPORT_SCAN_TYPE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LST_IMPORT_SCAN_TYPE contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LST_IMPORT_SCAN_TYPE


% --- Executes during object creation, after setting all properties.
function LST_IMPORT_SCAN_TYPE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LST_IMPORT_SCAN_TYPE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in LST_IMPORT_FILE_TYPE.
function LST_IMPORT_FILE_TYPE_Callback(hObject, eventdata, handles)
% hObject    handle to LST_IMPORT_FILE_TYPE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LST_IMPORT_FILE_TYPE contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LST_IMPORT_FILE_TYPE


% --- Executes during object creation, after setting all properties.
function LST_IMPORT_FILE_TYPE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LST_IMPORT_FILE_TYPE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in LST_PHS.
function LST_PHS_Callback(hObject, eventdata, handles)
% hObject    handle to LST_PHS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns LST_PHS contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LST_PHS

UpdatePreviewImage('LST_PHS');


function UpdatePreviewImage(FromWhich)
% This function updates the time frame slider and list, then the image.

hImportHARPGui = getappdata(0,'hImportHARPGui');
hSLDR_PHS = findobj(hImportHARPGui,'Tag','SLDR_PHS');
hLST_PHS = findobj(hImportHARPGui,'Tag','LST_PHS');

% Update frame slider or list according to each other.
if(strcmp(FromWhich,'SLDR_PHS'))
    set(hLST_PHS,'Value',int32(get(hSLDR_PHS,'Value')));
else
    set(hSLDR_PHS,'Value',get(hLST_PHS,'Value'));
end

UpdateImageShow;


function UpdateImageShow(WhichAxes)

SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

hImportHARPGui = getappdata(0,'hImportHARPGui');
hAXIS_IMG = getappdata(0,'hAXIS_IMG');
ColorRange = getappdata(0,'ColorRange');

hLST_SLICES = findobj(hImportHARPGui,'Tag','LST_SLICES');
CurrentSlice = get(hLST_SLICES,'Value');

hLST_DYNAMICS = findobj(hImportHARPGui,'Tag','LST_DYNAMICS');
CurrentDynamic = get(hLST_DYNAMICS,'Value');

hLST_PHS = findobj(hImportHARPGui,'Tag','LST_PHS');
CurrentPhase = get(hLST_PHS,'Value');

hLST_PRVWTYPE = findobj(hImportHARPGui,'Tag','LST_PRVWTYPE');
ReqPreviewType = get(hLST_PRVWTYPE,'Value');
AvailShowList = get(hLST_PRVWTYPE,'String');

hSLDR_IMG_BRIGHT = findobj(hImportHARPGui,'Tag','SLDR_IMG_BRIGHT');
IntensityThresh = 100 - get(hSLDR_IMG_BRIGHT,'Value');

MyData = squeeze(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Data(CurrentPhase,:,:));

TLC = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TLC;
NXY = SliceList{CurrentSlice}.NXY;

figure(hImportHARPGui);

switch AvailShowList{ReqPreviewType}
    case 'Magnitude'
        if isreal(MyData) % To allow negative values for CSPAMM %-SS- 09/28/10
            MyImage = (MyData);
        else
            MyImage = abs(MyData);
        end
        ImClrMap=gray;
    case 'ROI Magnitude'
        MyImage = squeeze(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Data(CurrentPhase,TLC(1)+(1:NXY(1)),TLC(2)+(1:NXY(2))));
        if isreal(MyImage) % To allow negative values for CSPAMM %-SS- 09/28/10
            MyImage = (MyImage);
        else
            MyImage = abs(MyImage);
        end
        ImClrMap=gray;%ImRange=[ ];
    case 'Spectrum'
        MyImage=(abs(fftshift(fft2((MyData)))));ImClrMap=gray;%ImRange=[];
    case 'Checkerboard' % Working on adding this -JLP- 12/7/11
        if CurrentDynamic == 1
            OtherDynamic = 2;
        else
            OtherDynamic = 1;
        end
        MyData2=squeeze(SliceList{CurrentSlice}.dynamic{OtherDynamic}.Data(CurrentPhase,:,:));
        if isreal(MyData)
            MyImage = (MyData.*MyData2);
        else
            MyImage = abs(MyData.*MyData2);
        end
        ImClrMap=gray;
    case 'Phase'
        MyImage=angle(MyData);ImClrMap=gray;%ImRange=[];
    case 'Real'
        MyImage=real(MyData);ImClrMap=gray;%ImRange=[];
    case 'Imag'
        MyImage=imag(MyData);ImClrMap=gray;%ImRange=[];
    case 'HarpPhase1'
        %         MyImage=cos(squeeze(SliceList{CurrentSlice}.HarpPhase1(CurrentPhase,:,:)));ImClrMap=gray;%ImRange=[];
        MyImage=wrap(0.5*wrap(squeeze(SliceList{CurrentSlice}.HarpPhase1(CurrentPhase,:,:))));ImClrMap=gray;%ImRange=[];
    case 'HarpPhase2'
        %         MyImage=cos(squeeze(SliceList{CurrentSlice}.HarpPhase2(CurrentPhase,:,:)));ImClrMap=gray;%ImRange=[];
        MyImage=wrap(0.5*wrap(squeeze(SliceList{CurrentSlice}.HarpPhase2(CurrentPhase,:,:))));ImClrMap=gray;%ImRange=[];
    case 'HarpMagnitude'
        MyImage=squeeze(SliceList{CurrentSlice}.HarpMagnitude(CurrentPhase,:,:));ImClrMap=gray;%ImRange=[];
    case 'Ecc'
        MyImage=squeeze(SliceList{CurrentSlice}.Strain2DInfoSet.EC_STRAIN(CurrentPhase,:,:));ImClrMap=jet;%ImRange=[-0.5 0.5];
    case 'Err'
        MyImage=squeeze(SliceList{CurrentSlice}.Strain2DInfoSet.ER_STRAIN(CurrentPhase,:,:));ImClrMap=jet;%ImRange=[-0.5 0.5];
    case 'Exx'
        MyImage=squeeze(SliceList{CurrentSlice}.Strain2DInfoSet.EX_STRAIN(CurrentPhase,:,:));ImClrMap=jet;%ImRange=[-0.5 0.5];
    case 'Eyy'
        MyImage=squeeze(SliceList{CurrentSlice}.Strain2DInfoSet.EY_STRAIN(CurrentPhase,:,:));ImClrMap=jet;%ImRange=[-0.5 0.5];
    case 'Ep1_2D'
        MyImage=squeeze(SliceList{CurrentSlice}.Strain2DInfoSet.lambda(CurrentPhase,:,:,1));ImClrMap=jet;%ImRange=[-0.5 0.5];
    case 'Ep2_2D'
        MyImage=squeeze(SliceList{CurrentSlice}.Strain2DInfoSet.lambda(CurrentPhase,:,:,2));ImClrMap=jet;%ImRange=[-0.5 0.5];
    case 'CircAng1_2D'
        MyImage=squeeze(SliceList{CurrentSlice}.Strain2DInfoSet.angle(CurrentPhase,:,:,1))/pi*180;ImClrMap=jet;%ImRange=[0 90];
    case 'CircAng2_2D'
        MyImage=squeeze(SliceList{CurrentSlice}.Strain2DInfoSet.angle(CurrentPhase,:,:,2))/pi*180;ImClrMap=jet;%ImRange=[0 90];
end

CKBUpSample = 1;
% -SS- 01/25/10
% Creating a mask for visualization
% Based on the slider value (user have control over this value)
Mask = ones(size(MyImage));
HideMask = sum(strcmpi(AvailShowList{ReqPreviewType},[{'Magnitude'},{'Spectrum'},{'Checkerboard'},{'Phase'},{'Real'},{'Imag'}]));

if isfield(SliceList{CurrentSlice},'HarpMagnitude') && ~HideMask
    
    hSLDR_IMG_MASK = findobj(hImportHARPGui,'Tag','SLDR_IMG_MASK');
    
    MImg = squeeze(SliceList{CurrentSlice}.HarpMagnitude(CurrentPhase,:,:));
    TH = (max(max(MImg)) - min(min(MImg)))*(get(hSLDR_IMG_MASK,'Value'));
    Mask(logical(MImg < TH)) = NaN;
end

Mask = repmat(Mask,[1,1,3]);

% Updating ImRange
% -MVM- 07/30/10
% Brightness bar added.
if isfield(ColorRange,AvailShowList{ReqPreviewType})
    ImRange = eval(['ColorRange.' AvailShowList{ReqPreviewType}]);
    MyImage(logical(MyImage<ImRange(1))) = ImRange(1);
    MyImage(logical(MyImage>ImRange(2))) = ImRange(2);
else
    ImRange = [min(MyImage(:)) max(MyImage(:))];
    ImRange = ImRange*(IntensityThresh/100);
end
TempImageRGB=ind2rgb(ceil((MyImage-ImRange(1))/...
    (ImRange(2)-ImRange(1)+eps)*size(ImClrMap,1)),ImClrMap);

TempImageRGB = TempImageRGB.*Mask;

axes(hAXIS_IMG); %#ok<MAXES>
if nargin ~= 0
    axes(WhichAxes);hold off; %#ok<MAXES>
    axis off;
end

switch AvailShowList{ReqPreviewType}
    case {'Magnitude';'Spectrum';'Checkerboard';'Phase';'Real';'Imag'}
        imagesc(TempImageRGB);
        colormap(ImClrMap);
    otherwise
        imagesc(TempImageRGB);
        axis equal; axis tight;
        colormap(ImClrMap);
        hPNL_SHOW_TRACKS = findobj(hImportHARPGui,'Tag','PNL_SHOW_TRACKS');
        if strcmp(SliceList{CurrentSlice}.ProcessInfo.TrackStatus.flag,'clean')
            hold on;
            hObject = get(hPNL_SHOW_TRACKS,'SelectedObject');
            Selection = get(hObject, 'String');
            if strcmp(Selection,'HARP') && ...
                    strcmp(SliceList{CurrentSlice}.ProcessInfo.TrackStatus.Harp.flag,'clean')
                %using traditional tracking
                TrackData = squeeze(SliceList{CurrentSlice}.HarpTrackInfoSet.TrackForwardData(CurrentPhase,:,:,:));
            elseif strcmp(Selection,'HARP Refine') && ...
                    strcmp(SliceList{CurrentSlice}.ProcessInfo.TrackStatus.Refine.flag,'clean')
                %using refinement
                MotionForwardX = SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardX;
                MotionForwardY = SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardY;
                
                TrackData(:,:,1) = meshgrid(1:NXY(1),1:NXY(2))' + squeeze(sum(MotionForwardX(1:CurrentPhase-1,:,:),1)); %-SYM- 07/12/2010
                TrackData(:,:,2) = meshgrid(1:NXY(2),1:NXY(1))  + squeeze(sum(MotionForwardY(1:CurrentPhase-1,:,:),1));
            else
                TrackData = NaN(NXY(1),NXY(2),2);
            end
            
            % creating the mask
            % NOTE: --important-- mask should always be created on first phase for tracking
            mask = ones(NXY(1),NXY(2));
            MImg = squeeze(SliceList{CurrentSlice}.HarpMagnitude(1,:,:));
            TH = (max(max(MImg)) - min(min(MImg)))*(get(hSLDR_IMG_MASK,'Value'));
            mask(logical(MImg < TH)) = NaN;
            
            TrackData(:,:,1) = TrackData(:,:,1).*mask;
            TrackData(:,:,2) = TrackData(:,:,2).*mask;
            TrackData_coordX = TrackData(1:2:end,1:2:end,2);
            TrackData_coordY = TrackData(1:2:end,1:2:end,1);
            indexTrackData_coord = find(0<TrackData_coordX <= NXY(2) & 0<TrackData_coordY <= NXY(1));
            plot(TrackData_coordX(indexTrackData_coord)*CKBUpSample,TrackData_coordY(indexTrackData_coord)*CKBUpSample,...
                'r.','MarkerSize',5);
            axis equal;
            xlim([0.5,NXY(1)+0.5]); ylim([0.5,NXY(2)+0.5]);
            hold off;
        end
end

UpdatePreviewROIGUI;


function EDT_CLR_RANGE_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_CLR_MAX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_CLR_MAX as text
%        str2double(get(hObject,'String')) returns contents of EDT_CLR_MAX as a double
hImportHARPGui  = getappdata(0,'hImportHARPGui');
ColorRange      = getappdata(0,'ColorRange');

hLST_PRVWTYPE   = findobj(hImportHARPGui,'Tag','LST_PRVWTYPE');
ReqPreviewType  = get(hLST_PRVWTYPE,'Value');
AvailShowList   = get(hLST_PRVWTYPE,'String');

hEDT_CLR_RANGE  = findobj(hImportHARPGui,'Tag','EDT_CLR_RANGE');

[x,ok] = str2num(get(hEDT_CLR_RANGE,'string'));
if(ok && (length(x)==2)) && (x(1) < x(2))
    evalc(['ColorRange.' AvailShowList{ReqPreviewType} ' = x']);
    
else
    errordlg('Please input a valid range!');
    return;
end

setappdata(0,'ColorRange',ColorRange);

Update_Colorbar;
UpdatePropertiesGUI;
UpdateImageShow;


% --- Executes on button press in BTN_EXPRT_ALL_IMG.
function BTN_EXPRT_ALL_IMG_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_EXPRT_ALL_IMG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SliceList = getappdata(0,'SliceList');
if isempty(SliceList)
    return;
end

hImportHARPGui = getappdata(0,'hImportHARPGui');
hAXIS_IMG = getappdata(0,'hAXIS_IMG');

hLST_SLICES = findobj(hImportHARPGui,'Tag','LST_SLICES');
CurrentSlice = get(hLST_SLICES,'Value');

hLST_PHS = findobj(hImportHARPGui,'Tag','LST_PHS');

hEDT_PATH = findobj(hImportHARPGui,'Tag','EDT_PATH');
pathname = get(hEDT_PATH,'string');
filename = 'new';

hLST_PRVWTYPE = findobj(hImportHARPGui,'Tag','LST_PRVWTYPE');
ReqPreviewType  = get(hLST_PRVWTYPE,'Value');
AvailShowList = get(hLST_PRVWTYPE,'String');

if isfield(SliceList{CurrentSlice},'DatasetInfo')
    % check if there is a date extension
    % date format: 'yyyy-mm-dd'
    filename = SliceList{CurrentSlice}.DatasetInfo.FileName;
    [match,start]=regexpi(filename,'-(\d{4})-(\d{2})-(\d{2})','match','start');
    if ~isempty(match)
        filename = filename(1:start-1);
    end
end

if isempty(pathname)
    pathname = pwd;
end

filename = [filename '-' AvailShowList{ReqPreviewType} '-Slc' num2str(CurrentSlice) '.png'];

[filename, pathname] = uiputfile(...
    {'*.png','PNG Files (*.png)';...
    '*.jpg','JPEG Files (*.jpg)';...
    '*.fig','Matlab Figures (*.fig)'},...
    'Please select a folder... (empty folder recommended)',[pathname filesep filename]);
if isequal(filename,0) || isequal(pathname,0)
    % User pressed cancel
    return;
end

NoPhases = SliceList{CurrentSlice}.NO_FRAMES;
WaitBar = waitbar(0,'Saving images. Please wait...');

for my_index=1:NoPhases
    set(hLST_PHS,'Value',my_index);
    UpdatePreviewImage('LST_PHS');
    %     UpdateImageShow(gca(h));
    %     axis image;
    %     axis off;
    %[indImage,Map] = rgb2ind(frame2im(getframe(hAXIS_IMG)),1024);
    movegui(hAXIS_IMG); %-SYM- 07/14/2010 Use this if user uses dual screen
    % display and figure is on 2nd screen
    [indImage,Map] = rgb2ind(frame2im(getframe),1024);
    fullName = [pathname filename(1:end-4) '-' num2str(my_index) filename(end-3:end)];
    waitbar(my_index/NoPhases,WaitBar);
    imwrite(indImage,Map,fullName);
end
close(WaitBar);


% --- Executes on button press in BTN_EXPRT_MOV.
function BTN_EXPRT_MOV_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_EXPRT_MOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SliceList = getappdata(0,'SliceList');
if isempty(SliceList)
    return;
end

hImportHARPGui = getappdata(0,'hImportHARPGui');
hAXIS_IMG = getappdata(0,'hAXIS_IMG');

hLST_SLICES = findobj(hImportHARPGui,'Tag','LST_SLICES');
CurrentSlice = get(hLST_SLICES,'Value');

hLST_PHS = findobj(hImportHARPGui,'Tag','LST_PHS');

filename = 'new';
hEDT_PATH = findobj(hImportHARPGui,'Tag','EDT_PATH');
pathname = get(hEDT_PATH,'string');

hLST_PRVWTYPE   = findobj(hImportHARPGui,'Tag','LST_PRVWTYPE');
ReqPreviewType  = get(hLST_PRVWTYPE,'Value');
AvailShowList   = get(hLST_PRVWTYPE,'String');

if isfield(SliceList{CurrentSlice},'DatasetInfo')
    % check if there is a date extension
    % date format: 'yyyy-mm-dd'
    filename = SliceList{CurrentSlice}.DatasetInfo.FileName;
    [match,start]=regexpi(filename,'-(\d{4})-(\d{2})-(\d{2})','match','start');
    if ~isempty(match)
        filename = filename(1:start-1);
    end
end

if isempty(pathname)
    pathname = pwd;
end

filename = [filename '-' AvailShowList{ReqPreviewType} '-Slc' num2str(CurrentSlice) '.gif'];

[filename, pathname] = uiputfile(...
    {'*.gif','GIF Files (*.gif)';...
    '*.avi','AVI Files (*.avi)'},...
    'Please select a folder...',[pathname filesep filename]);
if isequal(filename,0) || isequal(pathname,0)
    % User pressed cancel
    return;
end
[~,ext] = strtok(filename,'.');

NoPhases = SliceList{CurrentSlice}.NO_FRAMES;
WaitBar = waitbar(0,'Saving video. Please wait...');

if isequal(ext,'.gif') == 1
    for my_index = [1:NoPhases NoPhases-1:-1:2]
        
        set(hLST_PHS,'Value',my_index);
        UpdatePreviewImage('LST_PHS');
        movegui(hAXIS_IMG); %-SYM- 07/14/2010 Use this if user uses dual screen
        % display and figure is on 2nd screen
        [indImage,Map] = rgb2ind(frame2im(getframe(hAXIS_IMG)),1024);
        fullName = [pathname filename];
        
        if my_index == 1
            imwrite(indImage,Map,fullName,'DelayTime',0.05);
        else
            imwrite(indImage,Map,fullName,'WriteMode','append','DelayTime',0.05);
        end
        waitbar(my_index/NoPhases,WaitBar);
    end
    close(WaitBar);
    
elseif isequal(ext,'.avi') == 1
    frameIDX = 1;
    for my_index = [1:NoPhases NoPhases-1:-1:2]
        
        set(hLST_PHS,'Value',my_index);
        UpdatePreviewImage('LST_PHS');
        movegui(hAXIS_IMG); %-SYM- 07/14/2010 Use this if user uses dual screen
        % display and figure is on 2nd screen
        [indImage,Map] = rgb2ind(frame2im(getframe(hAXIS_IMG)),1024);
        fullName = [pathname filename];
        
        if my_index == 1
            %imshow(indImage,Map);
            axis off;
            FILM(frameIDX) = getframe;
        else
            %imshow(indImage,Map);
            axis off;
            FILM(frameIDX) = getframe;
        end
        frameIDX = frameIDX + 1;
        waitbar(my_index/NoPhases,WaitBar);
    end
    movie2avi(FILM,fullName,'compression','None','fps',15);
    close(WaitBar);
else
    errordlg('Please choose either GIF or AVI format.');
    close(WaitBar);
end


% --- Executes on selection change in LST_PRVWTYPE.
function LST_PRVWTYPE_Callback(hObject, eventdata, handles)
% hObject    handle to LST_PRVWTYPE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns LST_PRVWTYPE contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LST_PRVWTYPE

Update_Colorbar;
UpdateImageShow;
UpdatePropertiesGUI;


function Update_Colorbar

% This function updates the colorbar and color ranges of main figure
hImportHARPGui   = getappdata(0,'hImportHARPGui');
hAXIS_CLR_LIMITS = getappdata(0,'hAXIS_CLR_LIMITS');

ColorRange       = getappdata(0,'ColorRange');

hLST_PRVWTYPE   = findobj(hImportHARPGui,'Tag','LST_PRVWTYPE');
ReqPreviewType  = get(hLST_PRVWTYPE,'Value');
AvailShowList   = get(hLST_PRVWTYPE,'String');

if isfield(ColorRange,AvailShowList{ReqPreviewType})
    ImRange = eval(['ColorRange.' AvailShowList{ReqPreviewType}]);
    set(findobj(hImportHARPGui,'Tag','EDT_CLR_RANGE'),'Visible','on');
    set(findobj(hImportHARPGui,'Tag','AXIS_CLR_LIMITS'),'Visible','on');
    
    %updating colorbar
    axes(hAXIS_CLR_LIMITS) %#ok<MAXES>
    
    MIN_CLR = ImRange(1);
    MAX_CLR = ImRange(2);
    inc     = (MAX_CLR - MIN_CLR)/50;
    imagesc(ones(1,size((MIN_CLR:inc:MAX_CLR),2)),(MIN_CLR:inc:MAX_CLR)',(MIN_CLR:inc:MAX_CLR)');
    set(gca,'Clim',[MIN_CLR,MAX_CLR]);
    set(gca,'XTick',[],'XTickLabel',[],'Ylim',[MIN_CLR MAX_CLR],...
        'YAxisLocation','right','YDir','normal');
    %axes(hAXIS_IMG);
else
    axis(hAXIS_CLR_LIMITS,'off');
    set(findobj(hImportHARPGui,'Tag','EDT_CLR_RANGE'),'Visible','off');
    set(findobj(hImportHARPGui,'Tag','AXIS_CLR_LIMITS'),'Visible','off');
end


function EDT_NXY_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_NXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_NXY as text
%        str2double(get(hObject,'String')) returns contents of EDT_NXY as a double

SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end
UpdateSliceAndDynamicProperties;
UpdatePreviewROIGUI;


function EDT_TLC_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_TLC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_TLC as text
%        str2double(get(hObject,'String')) returns contents of EDT_TLC as a double

SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end
UpdateSliceAndDynamicProperties;
UpdatePreviewROIGUI;


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function FIG_IMPORTHARP_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to FIG_IMPORTHARP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SliceList       = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

hImportHARPGui  = getappdata(0,'hImportHARPGui');
hLST_SLICES     = findobj(hImportHARPGui,'Tag','LST_SLICES');
CurrentSlice    = get(hLST_SLICES,'Value');
hLST_DYNAMICS   = findobj(hImportHARPGui,'Tag','LST_DYNAMICS');
CurrentDynamic  = get(hLST_DYNAMICS,'Value');
hLST_PRVWTYPE   = findobj(hImportHARPGui,'Tag','LST_PRVWTYPE');
PreviewType       = get(hLST_PRVWTYPE,'Value');
AvailShowList   = get(hLST_PRVWTYPE,'String');

hLST_PHS        = findobj(hImportHARPGui,'Tag','LST_PHS');
CurrentPhase    = get(hLST_PHS,'Value');

nx       = SliceList{CurrentSlice}.ReconResolution(1);
ny       = SliceList{CurrentSlice}.ReconResolution(2);

numPhase = size(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Data, 1);
BTN_STATUS=get(hImportHARPGui,'selectiontype');

cp=get(handles.AXIS_IMG,'CurrentPoint');
cx=cp(1); cy=cp(3);

AllImage = sum(strcmpi(AvailShowList{PreviewType},{'magnitude','phase','real','imag'}));
Spectrum = strcmpi(AvailShowList{PreviewType},'Spectrum');
FOVImage = ~ (AllImage || Spectrum);

if AllImage && ... % showing image
        (cx<=0 || cx>=nx || cy<=0 || cy>=ny) % if cp is outside the ROI
    
    AXIS_IMG_Pressed.IsTrue = false;
    AXIS_IMG_Pressed.MouseDownStatus     = 0;
    setappdata(hImportHARPGui,'AXIS_IMG_Pressed',AXIS_IMG_Pressed);
    return;
end

NXY = SliceList{CurrentSlice}.NXY;
CKBUpSample = 1;
if FOVImage &&... % showing image
        (cx<=0 || cx>=NXY(2)*CKBUpSample || cy<=0 || cy>=NXY(1)*CKBUpSample) % if cp is outside the ROI
    AXIS_IMG_Pressed.IsTrue = false;
    AXIS_IMG_Pressed.MouseDownStatus     = 0;
    setappdata(hImportHARPGui,'AXIS_IMG_Pressed',AXIS_IMG_Pressed);
    return;
end

if Spectrum %Spectrum --> Draw Ellipse, ...
    OMEGA    = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Omega + SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Filter_Shift;
    RATIOS   = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Ratios;
    ROTATION = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Rotation * pi / 180;
    
    rot      = [cos(ROTATION) sin(ROTATION); -sin(ROTATION) cos(ROTATION)];
    F1_A     = rot*[RATIOS(1) 0; 0 RATIOS(2)]*[1;0]+OMEGA;
    F2_A     = rot*[RATIOS(1) 0; 0 RATIOS(2)]*[0;1]+OMEGA;
    
    hpos=get(handles.AXIS_IMG,'Position');
    POINT_SIZEY=1; %hpos(4)/ny;
    POINT_SIZEX=1; %hpos(3)/nx;
    cp=get(handles.AXIS_IMG,'CurrentPoint');
    cx=cp(1); cy=cp(3);
    cx=cx-nx/2-1;
    cy=cy-ny/2-1;
    
    AXIS_IMG_Pressed.IsTrue              = true;
    AXIS_IMG_Pressed.MouseDownStatus     = 0;
    
    if (norm([cy;cx]-F1_A)<4)
        AXIS_IMG_Pressed.MouseDownStatus = 4;     % Capture the square
    elseif (norm([cy;cx]-F2_A)<4)
        AXIS_IMG_Pressed.MouseDownStatus = 5;     % Capture the circle
    elseif (norm([cy;cx]-OMEGA)<4)
        AXIS_IMG_Pressed.MouseDownStatus = 6;     % Capture the origin
    end;
    
else%A whole/reduced FOV spatial image is shown --> Draw ROI for all others
    NXY         = SliceList{CurrentSlice}.NXY;
    TLC         = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TLC;
    ROI_CNTR    = SliceList{CurrentSlice}.ROI_CNTR;
    REF_NXY     = NXY;
    REF_TLC     = TLC;
    REF_ROI_CNTR= ROI_CNTR(:,CurrentPhase);
    
    p(1:2,1) = TLC;
    p(1:2,2) = TLC+[NXY(1)-1;0];
    p(1:2,3) = TLC+[NXY(1)-1;NXY(2)-1];
    p(1:2,4) = TLC+[0;NXY(2)-1];
    
    R1_A     = p(1:2,1); % top left corner, motion
    R2_A     = p(1:2,3); % bottom right corner, resizing
    R3_A     = R1_A+REF_ROI_CNTR; %Center
    
    AXIS_IMG_Pressed.IsTrue          = true;
    AXIS_IMG_Pressed.MouseDownStatus = 0;
    AXIS_IMG_Pressed.REF_NXY         = REF_NXY;
    AXIS_IMG_Pressed.REF_TLC         = REF_TLC;
    AXIS_IMG_Pressed.REF_ROI_CNTR    = REF_ROI_CNTR;
    
    hpos=get(handles.AXIS_IMG,'Position');
    POINT_SIZEY=1;%hpos(4)/ny;
    POINT_SIZEX=1;%hpos(3)/nx;
    cp=get(handles.AXIS_IMG,'CurrentPoint');
    cx=cp(1); cy=cp(3);
    %cx=(cx-1);%/POINT_SIZEX;
    %cy=(ny-cy);%/POINT_SIZEY;
    
    AXIS_IMG_Pressed.REF_POINT=[cy;cx];
    
    if (norm([cy;cx]-R1_A)<4)
        if(strcmp(BTN_STATUS,'normal'))
            AXIS_IMG_Pressed.MouseDownStatus = 1;     % Capture the square in the top left corner
        else
            AXIS_IMG_Pressed.MouseDownStatus =11;     % Capture the square in the top left corner with 'SHFT' or 'CTRL'
        end;
    elseif (norm([cy;cx]-R2_A)<4)
        AXIS_IMG_Pressed.MouseDownStatus = 2;     % Capture the circle bottom right corner
    elseif (norm([cy;cx]-R3_A)<4)
        AXIS_IMG_Pressed.MouseDownStatus = 3;     % Capture the center entrior of the square
    end;
end;

% if the button is not clicked on one of the above markers, it is taken
% as selecting the current point (seed point)
if AXIS_IMG_Pressed.MouseDownStatus == 0
    if AllImage % if the displayed is the whole image
        cx = cx-TLC(2);
        cy = cy-TLC(1);
    end
    if Spectrum ||... % if spectrum, do nothing
            cx<=0 || cx>=NXY(2)*CKBUpSample || cy<=0 || cy>=NXY(1)*CKBUpSample % if cp is outside the ROI
        % do nothing
    else
        pt = [cy;cx];% + TLC;%We save the currentPoint relative to the TLC. This is safer in case the ROI is not the
        %same in all dynamics
        %We reset all the track of the point to that selected point.
        SliceList{CurrentSlice}.CurrentPoint = repmat(pt,1,size(SliceList{CurrentSlice}.CurrentPoint,2));
        
        % track current point over all time phases
        if isfield(SliceList{CurrentSlice}, 'HarpPhase1')
            sHarp1 = size(SliceList{CurrentSlice}.HarpPhase1);
        else
            sHarp1 = zeros(0);
        end
        if isfield(SliceList{CurrentSlice}, 'HarpPhase2')
            sHarp2 = size(SliceList{CurrentSlice}.HarpPhase2);
        else
            sHarp2 = zeros(0);
        end
        
        if length(sHarp1) ~= 3 || sHarp1(1) ~= numPhase || sHarp1(2) ~= SliceList{CurrentSlice}.NXY(1) || ...
                sHarp1(3) ~= SliceList{CurrentSlice}.NXY(2)
            % do not track currentPoint
            
        else
            % track currentPoint using traditional HARP tracking
            dim = 2;
            if length(sHarp2) ~= 3 || (sHarp2(1) ~= sHarp1(1) || sHarp2(2) ~= sHarp1(2) || sHarp2(3) ~= sHarp1(3))
                dim = 1;
            end
            setappdata(0, 'SliceList', SliceList);
            disp('Point coordinates:'); disp(pt);
            curPt_track = TrackPoint(pt/CKBUpSample, CurrentSlice,CurrentPhase, dim, 0);
            SliceList   = getappdata(0, 'SliceList');
            SliceList{CurrentSlice}.CurrentPoint = curPt_track*CKBUpSample;
        end
        
        setappdata(0, 'SliceList', SliceList);
        UpdatePreviewROIGUI;
    end
end

setappdata(hImportHARPGui,'AXIS_IMG_Pressed',AXIS_IMG_Pressed);


% --- Executes on mouse motion over figure - except title and menu.
function FIG_IMPORTHARP_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to FIG_IMPORTHARP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

hImportHARPGui = getappdata(0,'hImportHARPGui');

hLST_SLICES     = findobj(hImportHARPGui,'Tag','LST_SLICES');
CurrentSlice    = get(hLST_SLICES,'Value');

hLST_DYNAMICS   = findobj(hImportHARPGui,'Tag','LST_DYNAMICS');
CurrentDynamic  = get(hLST_DYNAMICS,'Value');

hLST_PHS        = findobj(hImportHARPGui,'Tag','LST_PHS');
CurrentPhase    = get(hLST_PHS,'Value');

hLST_PRVWTYPE   = findobj(hImportHARPGui,'Tag','LST_PRVWTYPE');

if(isempty(hImportHARPGui)||~ishandle(hImportHARPGui))
    return;
end

AXIS_IMG_Pressed= getappdata(hImportHARPGui,'AXIS_IMG_Pressed');
NXY             =  SliceList{CurrentSlice}.NXY;
if(~(AXIS_IMG_Pressed.IsTrue))
    return;
end

% SliceList{CurrentSlice}.FilterResultStatus = 'dirty';
% SliceList{CurrentSlice}.StrainResultStatus = 'dirty';
% SliceList{CurrentSlice}.TracknResultStatus = 'dirty';

nx              = SliceList{CurrentSlice}.ReconResolution(1);
ny              = SliceList{CurrentSlice}.ReconResolution(2);
%hAXIS_IMG       = findobj(hImportHARPGui,'Type','axes'); %comment by -SS-
%01/05/10 because of colorbar axes
hAXIS_IMG       = getappdata(0,'hAXIS_IMG');

hpos            = get(hAXIS_IMG, 'Position');
POINT_SIZEY     = 1;%hpos(4)/ny;
POINT_SIZEX     = 1;%hpos(3)/nx;
cp              = get(hAXIS_IMG, 'CurrentPoint');

NO_FRAMES       = size(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Data,1);

if(get(hLST_PRVWTYPE,'Value')==2)%Spectrum--> Draw Ellipse, ...
    OMEGA    = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Omega + SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Filter_Shift;
    RATIOS   = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Ratios;
    ROTATION = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Rotation * pi / 180;
    
    cx              = cp(1);
    cy              = cp(3);
    cx=cx-nx/2-1;
    cy=cy-ny/2-1;
    
    if(AXIS_IMG_Pressed.MouseDownStatus == 4) % Reshaping and rotation of the ellipse
        ROTATION =-atan2(cx-OMEGA(2),cy-OMEGA(1))*180/pi;
        RATIOS(1)=abs(norm([cy;cx]-OMEGA));
        SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Ratios   = RATIOS;
        SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Rotation = ROTATION;
    elseif(AXIS_IMG_Pressed.MouseDownStatus == 5) % Resizing and rotation of the ellipse
        ROTATION = -atan2(cx-OMEGA(2),cy-OMEGA(1))*180/pi+90;
        newscale = abs(norm([cy;cx]-OMEGA))/RATIOS(2);
        RATIOS   = RATIOS*newscale;
        SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Ratios   = RATIOS;
        SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Rotation = ROTATION;
    elseif(AXIS_IMG_Pressed.MouseDownStatus == 6) % Moving the ellipse
        Pos    = [cy; cx];
        tOmega = SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Omega;
        SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Filter_Shift = Pos - tOmega;
    end;
else %WhichPreview --> Magnitude/phase
    REF_POINT = AXIS_IMG_Pressed.REF_POINT;
    cx              = cp(1);
    cy              = cp(3);
    %cx              = (cx-1)/POINT_SIZEX;
    %cy              = (ny*POINT_SIZEY-cy)/POINT_SIZEY;
    dp              = [cy;cx]-REF_POINT;
    
    REF_ROI_CNTR    = AXIS_IMG_Pressed.REF_ROI_CNTR;
    REF_TLC         = AXIS_IMG_Pressed.REF_TLC;
    REF_NXY         = AXIS_IMG_Pressed.REF_NXY;
    REF_POINT       = AXIS_IMG_Pressed.REF_POINT;
    
    if ((AXIS_IMG_Pressed.MouseDownStatus == 1)||(AXIS_IMG_Pressed.MouseDownStatus ==11)) %draging square in tlc moves TLC of ally dynamics
        tlc=round(REF_TLC+dp);
        if ((tlc+NXY(2)-1)>ny)
            tlc=ny-NXY(2)+1;
        end;
        if ((tlc+NXY(1)-1)>nx)
            tlc=nx-NXY(1)+1;
        end;
        
        TLC=tlc;
        tempTLC=SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TLC;
        
        iIndex=CurrentDynamic;%If SHFT/CTRL are pressed, TLC is changed for the current dyn. only
        if((AXIS_IMG_Pressed.MouseDownStatus == 1))%Otherwise TLC is displaced in all dynamics.
            iIndex=1:length(SliceList{CurrentSlice}.dynamic);
        end;
        
        
        for iIndex=iIndex
            SliceList{CurrentSlice}.dynamic{iIndex}.TLC=...
                SliceList{CurrentSlice}.dynamic{iIndex}.TLC+TLC-tempTLC;
        end;
        
    elseif (AXIS_IMG_Pressed.MouseDownStatus == 2)%This is the bottom right circle (to change NXY)
        center_p = REF_TLC+REF_NXY'/2;
        NXY1      = 2*max([max(abs(round(center_p(1)+REF_NXY(1)/2+dp(1))-center_p(1))),1]);
        NXY2      = 2*max([max(abs(round(center_p(2)+REF_NXY(2)/2+dp(2))-center_p(2))),1]);
        NXY      =[max([REF_ROI_CNTR(1) NXY1]) max([REF_ROI_CNTR(2) NXY2])];
        %         tlc      = (center_p - NXY/2);
        %
        %         if ((tlc+NXY-1)>ny)
        %             tlc=ny-NXY+1;
        %         end;
        %         if ((tlc+NXY-1)>nx)
        %             tlc=nx-NXY+1;
        %         end;
        %         TLC = tlc;
        %         tempTLC=squeeze(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TLC(:,CurrentPhase));
        %
        %         for iIndex=1:length(SliceList{CurrentSlice}.dynamic)
        %             SliceList{CurrentSlice}.dynamic{iIndex}.TLC=...
        %                 SliceList{CurrentSlice}.dynamic{iIndex}.TLC+TLC-tempTLC;
        %         end;
        
        %         SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TLC = TLC';
        SliceList{CurrentSlice}.NXY = NXY;
        
    elseif (AXIS_IMG_Pressed.MouseDownStatus == 3)%This is the center cross
        roi_cntr=round(REF_ROI_CNTR+dp);
        tmpIndex = roi_cntr>NXY';
        roi_cntr(tmpIndex)=NXY(tmpIndex')';
        roi_cntr(roi_cntr<0)=0;
        if ((roi_cntr(1)+NXY(2)-1)>ny)
            roi_cntr(1)=ny-NXY(2)+1;
        end;
        if ((roi_cntr(2)+NXY(1)-1)>nx)
            roi_cntr(2)=nx-NXY(1)+1;
        end;
        ROI_CNTR = roi_cntr;
        tempROI_CNTR=squeeze(SliceList{CurrentSlice}.ROI_CNTR(:,CurrentPhase));
        
        SliceList{CurrentSlice}.ROI_CNTR(:,CurrentPhase)=...
            SliceList{CurrentSlice}.ROI_CNTR(:,CurrentPhase)+ROI_CNTR-tempROI_CNTR;
        if(abs(SliceList{CurrentSlice}.RefCNTR_Phase1-CurrentPhase)<=...
                abs(SliceList{CurrentSlice}.RefCNTR_Phase2-CurrentPhase))
            SliceList{CurrentSlice}.RefCNTR_Phase1=CurrentPhase;
        else
            SliceList{CurrentSlice}.RefCNTR_Phase2=CurrentPhase;
        end;
        
        X=[SliceList{CurrentSlice}.RefCNTR_Phase1 SliceList{CurrentSlice}.RefCNTR_Phase2];
        if (SliceList{CurrentSlice}.NO_FRAMES > 1) %interpolation for differenet phases
            Y=[SliceList{CurrentSlice}.ROI_CNTR(1,X(1)) SliceList{CurrentSlice}.ROI_CNTR(1,X(2))];
            SliceList{CurrentSlice}.ROI_CNTR(1,:)=interp1(X,Y,(1:NO_FRAMES),'linear','extrap');
            
            Y=[SliceList{CurrentSlice}.ROI_CNTR(2,X(1)) SliceList{CurrentSlice}.ROI_CNTR(2,X(2))];
            SliceList{CurrentSlice}.ROI_CNTR(2,:)=interp1(X,Y,(1:NO_FRAMES),'linear','extrap');
        end
    end;
end %WhichPreview
setappdata(0,'SliceList',SliceList);
LST_DYNAMICS_Callback(hAXIS_IMG, eventdata, handles);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function FIG_IMPORTHARP_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to FIG_IMPORTHARP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SliceList = getappdata(0,'SliceList');
if isempty(SliceList)
    return;
end

hImportHARPGui = getappdata(0,'hImportHARPGui');
AXIS_IMG_Pressed = getappdata(hImportHARPGui,'AXIS_IMG_Pressed');

if AXIS_IMG_Pressed.IsTrue
    AXIS_IMG_Pressed.IsTrue = false;
    AXIS_IMG_Pressed.MouseDownStatus = 0;
    setappdata(0,'SliceList',SliceList);
end

UpdatePreviewROIGUI;
setappdata(hImportHARPGui,'AXIS_IMG_Pressed',AXIS_IMG_Pressed);


% --- Executes on selection change in LST_PRPGT_STNGS.
function LST_PRPGT_STNGS_Callback(hObject, eventdata, handles)
% hObject    handle to LST_PRPGT_STNGS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns LST_PRPGT_STNGS contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LST_PRPGT_STNGS
if ~get(hObject,'Value')
    set(hObject,'Value',1);
    return;
end

SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    set(hObject,'Value',1);
    return;
end

hImportHARPGui = getappdata(0,'hImportHARPGui');
hLST_SLICES = findobj(hImportHARPGui,'Tag','LST_SLICES');
CurrentSlice = get(hLST_SLICES,'Value');
nSlices = length(SliceList);
PropagateWhich = get(hObject,'Value') - 1;

switch PropagateWhich
    case 0
        % Do nothing
    case 1 % ROI -->> All of same scan-type & parallel
        for iIndex=1:nSlices,
            
            Norm1 = SliceList{iIndex}.Planes3D.norm;
            Norm2 = SliceList{CurrentSlice}.Planes3D.norm;
            RelAngle = acos(dot(Norm1,Norm2));
            
            if strcmp(SliceList{iIndex}.ScanType,SliceList{CurrentSlice}.ScanType) && RelAngle<(0.001*pi)
                nDynamics = length(SliceList{iIndex}.dynamic);
                SliceList{iIndex}.NXY = SliceList{CurrentSlice}.NXY;
                SliceList{iIndex}.ROI_CNTR = SliceList{CurrentSlice}.ROI_CNTR;
                for jIndex=1:nDynamics,
                    SliceList{iIndex}.dynamic{jIndex}.TLC = SliceList{CurrentSlice}.dynamic{jIndex}.TLC;
                end
                SliceList{iIndex}.ProcessInfo.FilterStatus.flag = 'dirty';
            end
            
        end
    case 2 % Coarse Fix for Center of Filter --> all
        % Added by -HA-
        for iIndex=1:nSlices
            NXY = SliceList{iIndex}.ReconResolution;
            FILTER_T = fftshift(design_filter([0 0], ...
                SliceList{iIndex}.dynamic{1}.Ratios,SliceList{iIndex}.dynamic{1}.Decay, ...
                SliceList{iIndex}.dynamic{1}.Rotation,NXY(1),NXY(2)));
            [Y,X] = meshgrid(0:NXY(2)-1,0:NXY(1)-1);
            Womega =  SliceList{iIndex}.ReconResolution .* ...
                SliceList{iIndex}.PxlSpacing/SliceList{iIndex}.TAG_Spacing;
            
            nDynamics = length(SliceList{iIndex}.dynamic);
            for jIndex=1:nDynamics,
                omega = Womega';
                testImage = squeeze(SliceList{iIndex}.dynamic{jIndex}.Data(1,:,:));
                % Check 4 omega
                bestMag = 0;
                IX = [1 0 -1]; if strcmp(SliceList{iIndex}.ScanType,'zHARP'), IX = [1 0]; end
                for ix = IX
                    for iy = IX
                        EXP = exp(-1i*2*pi*(X*omega(1)*ix/NXY(1) + Y*omega(2)*iy/NXY(2)));
                        testMag = sum(sum(abs(fft2(testImage.*EXP).*FILTER_T).^2));
                        if testMag > bestMag
                            baseix = ix;baseiy = iy;bestMag = testMag;
                        end
                    end
                end
                omega = omega.*[baseix baseiy]';
                SliceList{iIndex}.dynamic{jIndex}.Omega = omega;
                SliceList{iIndex}.dynamic{jIndex}.TAG_Angle = atan2(omega(1),omega(2));
                SliceList{iIndex}.dynamic{jIndex}.Filter_Shift = SliceList{iIndex}.dynamic{jIndex}.Filter_Shift*0;
            end
            SliceList{iIndex}.ProcessInfo.FilterStatus.flag = 'dirty';
        end
    case 3 % Fine Fix for Center of Filter --> all
        % Added by -HA-
        for iIndex=1:nSlices,
            NXY = SliceList{iIndex}.ReconResolution;
            nDynamics = length(SliceList{iIndex}.dynamic);
            for jIndex=1:nDynamics,
                omega = SliceList{iIndex}.dynamic{jIndex}.Omega;
                testImage = abs(fftshift(fft2(squeeze(SliceList{iIndex}.dynamic{jIndex}.Data(1,:,:)))));
                [Wx, Wy] = meshgrid(-NXY(2)/2:(NXY(2)/2 -1),-NXY(1)/2:(NXY(1)/2 -1));
                netMotion = 100;
                while(netMotion > .01)
                    FILTER_T = (design_filter(omega, ...
                        SliceList{iIndex}.dynamic{jIndex}.Ratios,6, ...
                        SliceList{iIndex}.dynamic{jIndex}.Rotation,NXY(1),NXY(2),'gausswin '));
                    omegaNew(1) = mean(mean(Wy.*testImage.*FILTER_T))/mean(mean(testImage.*FILTER_T));
                    omegaNew(2) = mean(mean(Wx.*testImage.*FILTER_T))/mean(mean(testImage.*FILTER_T));
                    netMotion = norm(omega-omegaNew');
                    omega = omegaNew';
                end
                SliceList{iIndex}.dynamic{jIndex}.Filter_Shift = omega-SliceList{iIndex}.dynamic{jIndex}.Omega;
            end
        end
end

setappdata(0,'SliceList',SliceList);
UpdateLST_SLICES;
UpdatePreviewROIGUI;
set(hObject,'Value',1);
UpdatePreviewROIGUI;


function EDT_RATIOS_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_RATIOS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_RATIOS as text
%        str2double(get(hObject,'String')) returns contents of EDT_RATIOS as a double

SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end
UpdateSliceAndDynamicProperties;
UpdatePreviewROIGUI;


function EDT_DECAY_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_DECAY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_DECAY as text
%        str2double(get(hObject,'String')) returns contents of EDT_DECAY as a double

SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end
UpdateSliceAndDynamicProperties;
UpdatePreviewROIGUI;


function EDT_ROTATION_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_ROTATION (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_ROTATION as text
%        str2double(get(hObject,'String')) returns contents of EDT_ROTATION as a double

SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end
UpdateSliceAndDynamicProperties;
UpdatePreviewROIGUI;


% --- Executes on button press in BTN_FLTER_ALL.
function BTN_FLTER_ALL_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_FLTER_ALL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

hImportHARPGui = getappdata(0,'hImportHARPGui');
SlicesWithFilterResults = 0;

SliceListLength = length(SliceList);
hWaitBar1 = waitbar(0,['Filtering a total number of ' num2str(SliceListLength) ' slices...']);

for SliceIndex = 1:SliceListLength
    
    waitbar(SliceIndex/SliceListLength,hWaitBar1,['Filtering slice ' num2str(SliceIndex) ' of ' num2str(SliceListLength) ' slices. Please wait...'])
    
    FilterSlice(SliceIndex);
    
    % Storing Filter Properties:
    SliceList = getappdata(0,'SliceList');
    if strcmp(SliceList{SliceIndex}.ProcessInfo.FilterStatus.flag,'clean')
        NumDyn = length(SliceList{SliceIndex}.dynamic);
        for DynIndex = 1:NumDyn
            SliceList{SliceIndex}.ProcessInfo.FilterStatus.Properties.Omega(:,DynIndex) = ...
                SliceList{SliceIndex}.dynamic{DynIndex}.Omega;
            SliceList{SliceIndex}.ProcessInfo.FilterStatus.Properties.Filter_Shift(:,DynIndex) = ...
                SliceList{SliceIndex}.dynamic{DynIndex}.Filter_Shift;
            SliceList{SliceIndex}.ProcessInfo.FilterStatus.Properties.Ratios(:,DynIndex) = ...
                SliceList{SliceIndex}.dynamic{DynIndex}.Ratios;
            SliceList{SliceIndex}.ProcessInfo.FilterStatus.Properties.Rotation(1,DynIndex) = ...
                SliceList{SliceIndex}.dynamic{DynIndex}.Rotation;
            SliceList{SliceIndex}.ProcessInfo.FilterStatus.Properties.Decay(1,DynIndex) = ...
                SliceList{SliceIndex}.dynamic{DynIndex}.Decay;
            SliceList{SliceIndex}.ProcessInfo.FilterStatus.Properties.TLC(:,DynIndex) = ...
                SliceList{SliceIndex}.dynamic{DynIndex}.TLC;
        end
        
        % trueHARP filter property
        if isfield(SliceList{SliceIndex},'dSampFactor')
            SliceList{SliceIndex}.ProcessInfo.FilterStatus.Properties.dSampFactor = ...
                SliceList{SliceIndex}.dSampFactor;
        end
        
        SliceList{SliceIndex}.ProcessInfo.FilterStatus.Properties.NXY = ...
            SliceList{SliceIndex}.NXY;
        SliceList{SliceIndex}.ProcessInfo.FilterStatus.Properties.TAG_Spacing = ...
            SliceList{SliceIndex}.TAG_Spacing;
        SliceList{SliceIndex}.ProcessInfo.FilterStatus.Properties.ScanType = ...
            SliceList{SliceIndex}.ScanType;
        
        SliceList{SliceIndex}.ProcessInfo.FilterStatus.time = datestr(clock,0);
        
        SlicesWithFilterResults = SlicesWithFilterResults + 1;
    end
    SliceList{SliceIndex}.ProcessInfo.Strain2DStatus.flag     = 'dirty';
    SliceList{SliceIndex}.ProcessInfo.TrackStatus.flag        = 'dirty';
    SliceList{SliceIndex}.ProcessInfo.TrackStatus.Harp.flag   = 'dirty';
    SliceList{SliceIndex}.ProcessInfo.TrackStatus.Refine.flag = 'dirty';
    
    setappdata(0,'SliceList',SliceList);
    if isfield(SliceList{SliceIndex},'Strain2DInfoSet')
        SliceList{SliceIndex} = rmfield(SliceList{SliceIndex},'Strain2DInfoSet');
    end
    SliceList{SliceIndex}.AvailResults = UpdateAvailableResults(SliceIndex);
    setappdata(0,'SliceList',SliceList);
    
end

close(hWaitBar1);

hLST_PRVWTYPE = findobj(hImportHARPGui,'Tag','LST_PRVWTYPE');
hLST_SLICES = findobj(hImportHARPGui,'Tag','LST_SLICES');
CurrentSlice = get(hLST_SLICES,'Value');
set(hLST_PRVWTYPE,'String',SliceList{CurrentSlice}.AvailResults);
if ~isfield(SliceList{CurrentSlice},'HarpMagnitude')
    set(hLST_PRVWTYPE,'Value',FindStringIndex(SliceList{CurrentSlice}.AvailResults,'ROI Magnitude'));
else
    set(hLST_PRVWTYPE,'Value',FindStringIndex(SliceList{CurrentSlice}.AvailResults,'HarpMagnitude'));
end

UpdatePreviewGUI;
UpdatePreviewROIGUI;
UpdatePropertiesGUI;


function MyList = AddToAvailableResults(MyList,New2Add)

New2AddLength   = length(New2Add);
for i = 1:New2AddLength
    OldLength = length(MyList);
    found = false;
    for j = 1:OldLength,
        if strcmp(MyList(j),New2Add(i))
            found=true;
            break;
        end
    end
    if ~found
        MyList(OldLength+1)=New2Add(i);
    end
end


function NewAvailResults = UpdateAvailableResults(SlcIndx)

SliceList  = getappdata(0,'SliceList');
AllResults = SliceList{SlcIndx}.AvailResults;

if ~strcmpi(SliceList{SlcIndx}.ProcessInfo.FilterStatus.flag,'clean')
    
    TempInd = 1:length(AllResults);
    Found = 0;
    for index = 1:length(AllResults)
        if strcmpi(AllResults{index},'HarpMagnitude')
            Found = index;
        end
    end
    % to be removed:
    %    [{'HarpMagnitude'};{'HarpPhase1'};...
    %    {'HarpPhase2'};{'HarpPhase3'};{'Inhomogeneity'};
    if Found
        TempInd([Found:(Found+4)]) = [];
        AllResults = {AllResults{TempInd}}';
    end
end

if ~strcmpi(SliceList{SlcIndx}.ProcessInfo.Strain2DStatus.flag,'clean')
    TempInd = 1:length(AllResults);
    Found = 0;
    for index = 1:length(AllResults)
        if strcmpi(AllResults{index},'Ecc')
            Found = index;
        end
    end
    % to be removed:
    %              [{'Ecc'};{'Err'};{'Exx'};{'Eyy'};{'zEcc'};{'zErr'};...
    %              {'zExx'};{'zEyy'};{'CircAng1_2D'};{'CircAng2_2D'};{'Ep1_2D'};{'Ep2_2D'}]);
    if Found
        TempInd([Found:(Found+11)]) = [];
        AllResults = {AllResults{TempInd}}';
    end
    
    TempInd = 1:length(AllResults);
    Found = 0;
    for index = 1:length(AllResults)
        if strcmpi(AllResults{index},'FA_2D')
            Found = index;
        end
    end
    % to be removed:
    %            [{'FA_2D'},{'Indx_2D'},{'MD_2D'}]
    if Found
        TempInd([Found:(Found+2)]) = [];
        AllResults = {AllResults{TempInd}}';
    end
end

NewAvailResults = AllResults;


function StrIndex = FindStringIndex(MyList, StringToFind)
StrLstLength   = length(MyList);
for i = 1:StrLstLength,
    if strcmp(MyList(i),StringToFind)
        break;
    end
end
StrIndex=i;


function EDT_FLTR_SHFT_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_FLTR_SHFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_FLTR_SHFT as text
%        str2double(get(hObject,'String')) returns contents of EDT_FLTR_SHFT as a double

SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end
UpdateSliceAndDynamicProperties;
UpdatePreviewROIGUI;


% --- Executes on button press in BTN_TRACK_HARP.
function BTN_TRACK_HARP_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_TRACK_HARP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the entire data structure that stores everything.
SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

SlicesWithTrackResults = 0;
hWaitBar1 = waitbar(0,['Conventional HARP tracking for ' num2str(length(SliceList)) ' slices. Please wait...']);

for SlcIndex = 1:length(SliceList)
    
    if ~strcmp(SliceList{SlcIndex}.ProcessInfo.FilterStatus.flag,'clean')
        close(hWaitBar1);
        errordlg('Please filter all slices first before tracking.');
        return;
    else
        
        RefPhase = 1;
        NXY = SliceList{SlcIndex}.NXY;
        
        xRange = 1:NXY(2);
        yRange = 1:NXY(1);
        TotalFrames = SliceList{SlcIndex}.NO_FRAMES;
        
        maskForReferenceFrame = ones(NXY(1),NXY(2));
        
        SliceList{SlcIndex}.HarpTrackInfoSet.TrackForwardData = ...
            permute(VectorTrack(SliceList{SlcIndex},RefPhase, xRange,...
            yRange, TotalFrames, maskForReferenceFrame, NXY),[3 1 2 4]);
        
        SlicesWithTrackResults = SlicesWithTrackResults + 1;
        
        SliceList{SlcIndex}.ProcessInfo.TrackStatus.flag = 'clean';
        SliceList{SlcIndex}.ProcessInfo.TrackStatus.Harp.flag = 'clean';
        SliceList{SlcIndex}.ProcessInfo.TrackStatus.Harp.RefPhase = RefPhase;
        SliceList{SlcIndex}.ProcessInfo.TrackStatus.Harp.time = datestr(clock,0);
        
    end
    waitbar(SlcIndex/length(SliceList),hWaitBar1);
end

close(hWaitBar1);
setappdata(0,'SliceList',SliceList);

UpdatePreviewROIGUI;
UpdatePropertiesGUI;


% --- Executes on button press in BTN_TRACK_REFINE.
function BTN_TRACK_REFINE_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_TRACK_REFINE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the entire data structure that stores everything.
SliceList       = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

hImportHARPGui = getappdata(0,'hImportHARPGui');

hLST_PHS = findobj(hImportHARPGui,'Tag','LST_PHS');
CurrentPhase = get(hLST_PHS,'Value');
str = 'HARP refinement will use current yellow circle(s) as seed(s).';
str = strcat(str, ' If you want to select different seed(s), please cancel.');
answer = questdlg(str, 'Seed Confirmation');

if ~strcmp(answer, 'Yes')
    return;
end
Sliceswithtrackresults = 0;

for SlcIndex = 1:length(SliceList)
    SliceList = getappdata(0,'SliceList');
    SliceList{SlcIndex}.ProcessInfo.TrackStatus.Refine.flag = 'dirty';
    setappdata(0,'SliceList',SliceList);
    
    if ~strcmp(SliceList{SlcIndex}.ProcessInfo.FilterStatus.flag,'clean')
        errordlg('Please filter all slices first before tracking.');
        return;
    else
        if(TrackRefineSlice(SlcIndex,CurrentPhase))
            
            Sliceswithtrackresults = Sliceswithtrackresults + 1;
            
            SliceList = getappdata(0,'SliceList');
            SliceList{SlcIndex}.ProcessInfo.TrackStatus.flag = 'clean';
            SliceList{SlcIndex}.ProcessInfo.TrackStatus.Refine.flag = 'clean';
            setappdata(0,'SliceList',SliceList);
        end
    end
end

setappdata(0,'SliceList',SliceList);
UpdatePreviewROIGUI;
UpdatePropertiesGUI;


function TrackStatus = TrackRefineSlice(WhichSlice,CurrentPhase)

%-SS- 01/28/10
% HarpTrackingRefinement.cpp and HarpTrackingRefinement.h is written by Xiaofeng Liu.

SliceList       = getappdata(0,'SliceList');
CurrentSlice    = WhichSlice;
numPhase        = size(SliceList{CurrentSlice}.dynamic{1}.Data, 1);
TrackStatus     = 0;

% If you have not HARP-filtered, then return.
if ~isfield(SliceList{CurrentSlice},'HarpPhase1')
    errordlg('HARP value is not available. Please filter the image first!');
    return;
end

sHarp1 = size(SliceList{CurrentSlice}.HarpPhase1);
sHarp2 = size(SliceList{CurrentSlice}.HarpPhase2);

if length(sHarp1) ~= 3
    msgbox('The size of field ''HarpPhase1'' is not right. Please check.');
    return;
end
if sHarp1(1) ~= numPhase || sHarp1(2) ~= SliceList{CurrentSlice}.NXY(1) || ...
        sHarp1(3) ~= SliceList{CurrentSlice}.NXY(2)
    msgbox('The size of field ''HarpPhase1'' is not right. Please check.');
    return;
end

dim = 2;
% You have a choice of doing 1-D or 2-D tracking. The default is 2-D.

if length(sHarp2) ~= 3 ...
        || (sHarp2(1) ~= sHarp1(1) || sHarp2(2) ~= sHarp1(2) || sHarp2(3) ~= sHarp1(3))
    answer = questdlg('HARP value is only available in one tag direction. Do you want to do 1D tracking?');
    
    if ~strcmp(answer, 'yes')
        return;
    end
    dim = 1;
end

% Find out whether the current tag is in x or y direction, if dim = 1
seed = SliceList{CurrentSlice}.CurrentPoint(:, CurrentPhase);
t1 = wrap(0.5*wrap(SliceList{CurrentSlice}.HarpPhase1(CurrentPhase, round(seed(1))+(-1:1),round(seed(2)))));
t2 = wrap(0.5*wrap(SliceList{CurrentSlice}.HarpPhase1(CurrentPhase, round(seed(1)),round(seed(2))+(-1:1))));
diff1 = abs(wrap(t1(2)-t1(1))) + abs(wrap(t1(3)-t1(2)));
diff2 = abs(wrap(t2(2)-t2(1))) + abs(wrap(t2(3)-t2(2)));
if diff1 > diff2
    dir1 = 1;
else
    dir1 = 2;
end

tagSep      = SliceList{CurrentSlice}.TAG_Spacing;
pixSep      = SliceList{CurrentSlice}.PxlSpacing;

% Calculate the reference phase
[y,x] = meshgrid(1:SliceList{CurrentSlice}.NXY(2), 1:SliceList{CurrentSlice}.NXY(1));

% Replaced by Fangxu based on Jonghye's modification. Source algorithm is in Xiaofeng's thesis page 134 Tecelao [74].
%     if dir1 == 1
%         RefPhase1 = wrap(x*SliceList{CurrentSlice}.RefPhaseInfo1.slope + SliceList{CurrentSlice}.RefPhaseInfo1.intercept);
%         RefPhase2 = wrap(y*SliceList{CurrentSlice}.RefPhaseInfo1.slope + SliceList{CurrentSlice}.RefPhaseInfo1.intercept);
%     else
%         RefPhase1 = wrap(y*SliceList{CurrentSlice}.RefPhaseInfo2.slope + SliceList{CurrentSlice}.RefPhaseInfo2.intercept);
%         RefPhase2 = wrap(x*SliceList{CurrentSlice}.RefPhaseInfo2.slope + SliceList{CurrentSlice}.RefPhaseInfo2.intercept);
%     end
% Later we found reference phase computed this way was not consistent for
% the computation of IDEA, so we ended up using the old code again to make
% up a reference phase. And for IDEA, we use the 1st time frame as the
% reference instead. 06/22/2011.
%     %%%%%%%%%% Below is the old code to calculate RefPhase %%%%%%%%%%%%
if dir1 == 1
    RefPhase1 = wrap(x*2*pi/(tagSep/pixSep(1)));
    RefPhase2 = wrap(y*2*pi/(tagSep/pixSep(2)));
else
    RefPhase1 = wrap(y*2*pi/(tagSep/pixSep(2)));
    RefPhase2 = wrap(x*2*pi/(tagSep/pixSep(1)));
end
%     %%%%%%%%%% Above is the old code to calculate RefPhase %%%%%%%%%%%%

% The reference phase is the set of straight tag lines that can be
% analytically got from the tag spacing.
SliceList{CurrentSlice}.RefPhase1 = RefPhase1;
SliceList{CurrentSlice}.RefPhase2 = RefPhase2;
setappdata(0, 'SliceList', SliceList);

% track the seed point;
seed_track = TrackPoint(seed, WhichSlice, CurrentPhase, dim, 0);

% At this point, I think we should show the user the track of this point,
% the make them confirm whether they want to use this seed point or not.

SliceList = getappdata(0,'SliceList');
SliceList{CurrentSlice}.CurrentPoint = seed_track;
SliceList{CurrentSlice}.ProcessInfo.TrackStatus.Refine.SeedPhase = CurrentPhase;

% Fangxu -JW- 1/9/2011 Track the seed to the reference time frame was
% added.
seed_cur = SliceList{CurrentSlice}.CurrentPoint(:,1); % 1st frame to ref frame
seed_ref = TrackPointHarp(seed_cur, ...
    squeeze(SliceList{CurrentSlice}.HarpPhase1(1, :,:)), ...
    RefPhase1, ...
    squeeze(SliceList{CurrentSlice}.HarpPhase2(1, :,:)), ...
    RefPhase2);

SliceList{CurrentSlice}.SeedInReference = seed_ref; % -JW- 1/9/2011 add new field (SeedInReference, point in reference frame)
NXY = SliceList{CurrentSlice}.NXY;

SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardX      = zeros(numPhase,NXY(1),NXY(2));
SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardY      = zeros(numPhase,NXY(1),NXY(2));
SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardIndex  = zeros(numPhase,NXY(1),NXY(2));
SliceList{CurrentSlice}.RefineTrackInfoSet.MotionBackwardX     = zeros(numPhase,NXY(1),NXY(2));
SliceList{CurrentSlice}.RefineTrackInfoSet.MotionBackwardY     = zeros(numPhase,NXY(1),NXY(2));
SliceList{CurrentSlice}.RefineTrackInfoSet.MotionBackwardIndex = zeros(numPhase,NXY(1),NXY(2));

% Now do the refinement slice by slice
% first track from the reference frame to all the other time frames.
hWaitbar = waitbar(0,['Performing HARP refinement on slice ' num2str(CurrentSlice) '. Please wait...']);
waitbar(0,hWaitbar);

for phase = 1:numPhase-1
    %phase
    seed_cur = SliceList{CurrentSlice}.CurrentPoint(:, phase);
    seed_next = SliceList{CurrentSlice}.CurrentPoint(:, phase+1);
    
    hx1 = wrap(0.5*wrap(squeeze(SliceList{CurrentSlice}.HarpPhase1(phase, :,:))));
    hx2 = wrap(0.5*wrap(squeeze(SliceList{CurrentSlice}.HarpPhase1(phase+1, :,:))));
    
    % made comment by -SS- 10/04/10
    %     if strcmp(SliceList{CurrentSlice}.ScanType, 'zHARP')
    hx1 = hx1*2;
    hx2 = hx2*2;
    %     end
    
    if dim == 1
        hy1 = RefPhase2;
        hy2 = RefPhase2;
    else
        hy1 = wrap(0.5*wrap(squeeze(SliceList{CurrentSlice}.HarpPhase2(phase, :,:))));
        hy2 = wrap(0.5*wrap(squeeze(SliceList{CurrentSlice}.HarpPhase2(phase+1, :,:))));
        %         if strcmp(SliceList{CurrentSlice}.ScanType, 'zHARP')
        hy1 = hy1*2;
        hy2 = hy2*2;
        %         end
    end
    
    mag1 = squeeze(SliceList{CurrentSlice}.HarpMagnitude(phase,:,:));
    mag2 = squeeze(SliceList{CurrentSlice}.HarpMagnitude(phase+1,:,:));
    
    %-SS- 01/28/10
    % Tracking is only possible in ROI since filter results are only stored
    % in ROI.
    
    %seed_ref_roi = seed_ref;
    seed_cur_roi = seed_cur;
    seed_next_roi = seed_next;
    
    %         [mx_roi, my_roi, index_roi] = HarpTrackingRefinement_FastMarching(...
    %             seed_cur_roi, ...
    %             hx1(TLC(1)+(1:NXY), TLC(2)+(1:NXY)), ...
    %             hx2(TLC(1)+(1:NXY), TLC(2)+(1:NXY)), ...
    %             hy1(TLC(1)+(1:NXY), TLC(2)+(1:NXY)), ...
    %             hy2(TLC(1)+(1:NXY), TLC(2)+(1:NXY)), ...
    %             seed_next_roi);
    %
    %         mx = zeros(size(RefPhase1));
    %         my = zeros(size(RefPhase1));
    %         index = zeros(size(RefPhase1));
    %         mx(TLC(1)+(1:NXY), TLC(2)+(1:NXY))      = mx_roi;
    %         my(TLC(1)+(1:NXY), TLC(2)+(1:NXY))      = my_roi;
    %         index(TLC(1)+(1:NXY), TLC(2)+(1:NXY))   = index_roi;
    [mx_roi, my_roi, index_roi] = HarpTrackingRefinement_FastMarching(...
        seed_cur_roi, ...
        hx1((1:NXY(1)), (1:NXY(2))), ...
        hx2((1:NXY(1)), (1:NXY(2))), ...
        hy1((1:NXY(1)), (1:NXY(2))), ...
        hy2((1:NXY(1)), (1:NXY(2))), ...
        mag1((1:NXY(1)), (1:NXY(2))), ...
        mag2((1:NXY(1)), (1:NXY(2))), ...
        seed_next_roi);
    
    mx = zeros(size(RefPhase1));
    my = zeros(size(RefPhase1));
    index = zeros(size(RefPhase1));
    mx((1:NXY(1)), (1:NXY(2)))      = mx_roi;
    my((1:NXY(1)), (1:NXY(2)))      = my_roi;
    index((1:NXY(1)), (1:NXY(2)))   = index_roi;
    
    SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardX(phase, :, :)      = mx;
    SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardY(phase, :, :)      = my;
    SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardIndex(phase, :, :)  = index;
    
    waitbar(phase/(3*numPhase), hWaitbar);
end
SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardX(numPhase,:,:) = 0;
SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardY(numPhase,:,:) = 0;
SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardIndex(phase, :, :)  = 0;

%  track from every time frame back to the reference frame
for phase = numPhase:-1:2
    %     phase
    seed_cur = SliceList{CurrentSlice}.CurrentPoint(:, phase);%-(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TLC)';
    seed_prev = SliceList{CurrentSlice}.CurrentPoint(:, phase-1);%-(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.TLC)';
    hx1 = wrap(0.5*wrap(squeeze(SliceList{CurrentSlice}.HarpPhase1(phase, :,:))));
    hx2 = wrap(0.5*wrap(squeeze(SliceList{CurrentSlice}.HarpPhase1(phase-1, :,:))));
    % made comment by -SS- 10/04/10
    %     if strcmp(SliceList{CurrentSlice}.ScanType, 'zHARP')
    hx1 = hx1*2;
    hx2 = hx2*2;
    %     end
    if dim == 1
        hy1 = RefPhase2;
        hy2 = RefPhase2;
    else
        hy1 = wrap(0.5*wrap(squeeze(SliceList{CurrentSlice}.HarpPhase2(phase, :,:))));
        hy2 = wrap(0.5*wrap(squeeze(SliceList{CurrentSlice}.HarpPhase2(phase-1, :,:))));
        %         if strcmp(SliceList{CurrentSlice}.ScanType, 'zHARP')
        hy1 = hy1*2;
        hy2 = hy2*2;
        %         end
    end
    
    mag1 = squeeze(SliceList{CurrentSlice}.HarpMagnitude(phase,:,:));
    mag2 = squeeze(SliceList{CurrentSlice}.HarpMagnitude(phase-1,:,:));
    
    %seed_ref_roi  = seed_ref;
    seed_cur_roi  = seed_cur;
    seed_prev_roi = seed_prev;
    
    %         [mx_roi, my_roi, index_roi] = HarpTrackingRefinement(...
    %             seed_cur_roi, ...
    %             hx1(TLC(1)+(1:NXY), TLC(2)+(1:NXY)), ...
    %             hx2(TLC(1)+(1:NXY), TLC(2)+(1:NXY)), ...
    %             hy1(TLC(1)+(1:NXY), TLC(2)+(1:NXY)), ...
    %             hy2(TLC(1)+(1:NXY), TLC(2)+(1:NXY)), ...
    %             seed_prev_roi);
    %
    %         mx = zeros(size(RefPhase1));
    %         my = zeros(size(RefPhase1));
    %         index = zeros(size(RefPhase1));
    %         mx(TLC(1)+(1:NXY), TLC(2)+(1:NXY))      = mx_roi;
    %         my(TLC(1)+(1:NXY), TLC(2)+(1:NXY))      = my_roi;
    %         index(TLC(1)+(1:NXY), TLC(2)+(1:NXY))   = index_roi;
    [mx_roi, my_roi, index_roi] = HarpTrackingRefinement_FastMarching(...
        seed_cur_roi, ...
        hx1((1:NXY(1)), (1:NXY(2))), ...
        hx2((1:NXY(1)), (1:NXY(2))), ...
        hy1((1:NXY(1)), (1:NXY(2))), ...
        hy2((1:NXY(1)), (1:NXY(2))), ...
        mag1((1:NXY(1)), (1:NXY(2))), ...
        mag2((1:NXY(1)), (1:NXY(2))), ...
        seed_prev_roi);
    
    mx = zeros(size(RefPhase1));
    my = zeros(size(RefPhase1));
    index = zeros(size(RefPhase1));
    mx((1:NXY(1)), (1:NXY(2)))      = mx_roi;
    my((1:NXY(1)), (1:NXY(2)))      = my_roi;
    index((1:NXY(1)), (1:NXY(2)))   = index_roi;
    
    
    SliceList{CurrentSlice}.RefineTrackInfoSet.MotionBackwardX(phase, :, :)      = mx;
    SliceList{CurrentSlice}.RefineTrackInfoSet.MotionBackwardY(phase, :, :)      = my;
    SliceList{CurrentSlice}.RefineTrackInfoSet.MotionBackwardIndex(phase, :, :)  = index;
    
    waitbar(0.33 + (numPhase-phase+2)/(3*numPhase), hWaitbar);
end

SliceList{CurrentSlice}.RefineTrackInfoSet.MotionBackwardX(1, :, :)      = 0;
SliceList{CurrentSlice}.RefineTrackInfoSet.MotionBackwardY(1, :, :)      = 0;
SliceList{CurrentSlice}.RefineTrackInfoSet.MotionBackwardIndex(1, :, :)  = 0;

%  Fangxu 06/13/2011 modified the following section to use the 1st time frame
%  as reference, since the current reference analysis brought from Harp3D
%  is unreliable.
%  JW 1/9/2011 track from every time frame to the reference frame to
%  calculate motionToreference
for phase = 1:numPhase
    
    seed_cur = SliceList{CurrentSlice}.CurrentPoint(:,phase);
    seed_ref = SliceList{CurrentSlice}.CurrentPoint(:,1);
    hx1 = wrap(0.5*wrap(squeeze(SliceList{CurrentSlice}.HarpPhase1(phase,:,:))));
    hx2 = wrap(0.5*wrap(squeeze(SliceList{CurrentSlice}.HarpPhase1(1,:,:)))); % by Fangxu
    %     hx2 = RefPhase1; % by Jonghye
    
    hx1 = hx1*2;
    hx2 = hx2*2;    % added by Fangxu brought from -SS- comment
    
    if dim == 1
        hy1 = RefPhase2;
        hy2 = RefPhase2; % these two lines for 1D tracking I'm not sure -Fangxu
    else
        hy1 = wrap(0.5*wrap(squeeze(SliceList{CurrentSlice}.HarpPhase2(phase,:,:))));
        hy2 = wrap(0.5*wrap(squeeze(SliceList{CurrentSlice}.HarpPhase2(1,:,:))));
        
        hy1 = hy1*2;
        hy2 = hy2*2;     % added by Fangxu brought from -SS- comment
        
    end
    
    mag1 = squeeze(SliceList{CurrentSlice}.HarpMagnitude(phase,:,:));
    mag2 = squeeze(SliceList{CurrentSlice}.HarpMagnitude(1,:,:)); % ref frame now is 1st time frame
    
    seed_cur_roi  = seed_cur;
    seed_ref_roi = seed_ref;
    
    [mx_roi, my_roi, index_roi] = HarpTrackingRefinement_FastMarching(...
        seed_cur_roi, ...
        hx1((1:NXY(1)), (1:NXY(2))), ...
        hx2((1:NXY(1)), (1:NXY(2))), ...
        hy1((1:NXY(1)), (1:NXY(2))), ...
        hy2((1:NXY(1)), (1:NXY(2))), ...
        mag1((1:NXY(1)), (1:NXY(2))), ...
        mag2((1:NXY(1)), (1:NXY(2))), ...
        seed_ref_roi);
    
    mx = zeros(size(RefPhase1));
    my = zeros(size(RefPhase1));
    index = zeros(size(RefPhase1));
    mx((1:NXY(1)), (1:NXY(2)))      = mx_roi;
    my((1:NXY(1)), (1:NXY(2)))      = my_roi;
    index((1:NXY(1)), (1:NXY(2)))   = index_roi;
    
    SliceList{CurrentSlice}.RefineTrackInfoSet.MotionToReferenceX(phase, :, :)      = mx;
    SliceList{CurrentSlice}.RefineTrackInfoSet.MotionToReferenceY(phase, :, :)      = my;
    SliceList{CurrentSlice}.RefineTrackInfoSet.MotionToReferenceIndex(phase, :, :)  = index;
    
    waitbar(0.67 + phase/(3*numPhase), hWaitbar);
end

close(hWaitbar);

setappdata(0, 'SliceList', SliceList);
TrackStatus = 1;


function [new_pt, motion] = TrackPointHarp(pt, hx1, hx2, hy1, hy2)
% Track one point using traditional HARP tracking over one frame

iter = 1;
y = pt;
[ny,nx] = size(hx1);

% the phase value of pt at the first time frame
current_A(1) = getLValue(hx1,pt);
current_A(2) = getLValue(hy1,pt);
h2 = [hx2, hy2];

ephi1 = ones(2,1);
while ((norm(ephi1(:,iter))>.0001)&&(iter<15))
    iter=iter+1;
    
    m=fix(y);
    dm=y-m;
    if (y(1)>=ny)
        m(1)=ny-1;dm(1)=.99;
    end;
    
    if (y(1)<1)
        m(1)=1; dm(1)=0;
    end;
    
    if (y(2)>=nx)
        m(2)=nx-1;dm(2)=.99;
    end;
    
    if (y(2)<1)
        m(2)=1; dm(2)=0;
    end;
    
    for k=1:2,
        at=-wrap(h2(m(1)+1,m(2)+(k-1)*nx)-h2(m(1),m(2)+(k-1)*nx))-...
            wrap(h2(m(1),m(2)+1+(k-1)*nx)-h2(m(1),m(2)+(k-1)*nx))+...
            wrap(h2(m(1)+1,m(2)+1+(k-1)*nx)-h2(m(1),m(2)+(k-1)*nx));
        bt=wrap(h2(m(1),m(2)+1+(k-1)*nx)-h2(m(1),m(2)+(k-1)*nx));
        ct=wrap(h2(m(1)+1,m(2)+(k-1)*nx)-h2(m(1),m(2)+(k-1)*nx));
        va=wrap(at*dm(1)*dm(2)+bt*dm(2)+ct*dm(1)+h2(m(1),m(2)+(k-1)*nx));
        ephi1(k,iter)=wrap(va-current_A(k));
        grad(k,1)=(at*dm(1)+bt);
        grad(k,2)=(at*dm(2)+ct);
    end;
    
    vt=-inv(grad'*grad)*grad'*ephi1(1:2,iter);
    
    if (abs(vt(1))>1),
        vt(1)=sign(vt(1));
    end;
    
    if (abs(vt(2))>1),
        vt(2)=sign(vt(2));
    end;
    y=y+[vt(2);vt(1)];
    
end % while

new_pt = y;
motion = new_pt - pt;


% --- Executes on button press in BTN_PLAY.
function BTN_PLAY_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_PLAY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SliceList = getappdata(0,'SliceList');
if isempty(SliceList)
    return;
end

hImportHARPGui = getappdata(0,'hImportHARPGui');

hLST_SLICES = findobj(hImportHARPGui,'Tag','LST_SLICES');
CurrentSlice  = get(hLST_SLICES,'Value');

hLST_DYNAMICS = findobj(hImportHARPGui,'Tag','LST_DYNAMICS');
CurrentDynamic = get(hLST_DYNAMICS,'Value');

hLST_PHS = findobj(hImportHARPGui,'Tag','LST_PHS');
NoPhases = size(SliceList{CurrentSlice}.dynamic{CurrentDynamic}.Data,1);

if NoPhases>1
    INDX=[1:2:NoPhases  (NoPhases-1):-2:1];
    for my_index=INDX(1:end)
        set(hLST_PHS,'Value',my_index);
        UpdatePreviewImage('LST_PHS');
    end
end


% --- Executes on button press in BTN_ANALYSIS.
function BTN_ANALYSIS_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_ANALYSIS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SliceList = getappdata(0,'SliceList');

if isempty(SliceList)
    return;
end

for i = 1:length(SliceList)
    if ~strcmp(SliceList{i}.ProcessInfo.TrackStatus.flag,'clean')
        errordlg('Use HARP or refinement to track all slices first!');
        return;
    end
end

visualizeHARP;


% --- Executes on button press in BTN_DEL_ALL_SLICES.
function BTN_DEL_ALL_SLICES_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_DEL_ALL_SLICES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

setappdata(0,'SliceList',[]);
UpdateLST_SLICES;
axes(handles.AXIS_IMG)
cla;


% --- Executes on button press in BTN_SHOW_3D_SLICES.
function BTN_SHOW_3D_SLICES_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_SHOW_3D_SLICES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

hImportHARPGui     = getappdata(0,'hImportHARPGui');
hLST_SLICES        = findobj(hImportHARPGui,'Tag','LST_SLICES');
SliceListName      = get(hLST_SLICES,'String');

ClrMap = cool(length(SliceList));

h = figure;
for SlcIndx = 1:length(SliceList)
    Planes3D=SliceList{SlcIndx}.Planes3D;
    x = [Planes3D.lt(1) Planes3D.rt(1) Planes3D.rb(1) Planes3D.lb(1) Planes3D.lt(1)]';
    y = [Planes3D.lt(2) Planes3D.rt(2) Planes3D.rb(2) Planes3D.lb(2) Planes3D.lt(2)]';
    z = [Planes3D.lt(3) Planes3D.rt(3) Planes3D.rb(3) Planes3D.lb(3) Planes3D.lt(3)]';
    %C = [1 1 1.0000 1 1.000];
    figure(h);hold on;fill3(x,y,z,ClrMap(SlcIndx,:));axis equal; axis vis3d;alpha(0.5);
    
    % If the slices in the LST_SLICES are not named based on position
    % in SliceList.
    
    SlcName = strtok(SliceListName{SlcIndx},'Slice');
    SlcStr  = SlcName(1:regexp(SlcName,'-')+1);
    % SlcName is of format "%d-scan type name" example: 1-zHarp
    
    % If slices are named as the position in SliceList uncomment the
    % following command:
    % SlcStr  = num2str(SlcIndx);
    
    text(Planes3D.lb(1),Planes3D.lb(2),Planes3D.lb(3),sprintf('lb %s',SlcStr),'Color','k');
    text(Planes3D.rt(1),Planes3D.rt(2),Planes3D.rt(3),sprintf('rt %s',SlcStr),'Color','b');
    text(Planes3D.origin(1),Planes3D.origin(2),Planes3D.origin(3),sprintf('org %s',SlcStr),'Color','k');
    plot3(Planes3D.center(1),Planes3D.center(2),Planes3D.center(3),'r*');
end


function PNL_SHOW_TRACKS_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to PNL_EL_TRACE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateImageShow;


% --- Executes on button press in BTN_EXPRT_IMG.
function BTN_EXPRT_IMG_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_EXPRT_IMG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = figure;
UpdateImageShow(gca(h));
axis image;
set(gcf,'color','w');axis off


% --- Executes during object creation, after setting all properties.
function EDT_CLR_RANGE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_CLR_RANGE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function SLDR_IMG_MASK_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SLDR_IMG_MASK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in BTN_TRACK_PT.
function BTN_TRACK_PT_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_TRACK_PT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Added by SYM 10/12/10
SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    return;
end

hImportHARPGui  = getappdata(0,'hImportHARPGui');
hLST_SLICES = findobj(hImportHARPGui,'Tag','LST_SLICES');
CurrentSlice = get(hLST_SLICES,'Value');

hLST_PHS = findobj(hImportHARPGui,'Tag','LST_PHS');
CurrentPhase = get(hLST_PHS,'Value');

dim = 2;
seed = SliceList{CurrentSlice}.CurrentPoint(:, CurrentPhase);
trackedPoints = TrackPoint(seed, CurrentSlice, CurrentPhase, dim, 0);
disp('Point track through all time frames:'); disp(trackedPoints);

totalNumberPhases = max(1:size(SliceList{CurrentSlice}.dynamic{1}.Data,1));
hImageAxes = getappdata(0,'hAXIS_IMG');

hold on;
if trackedPoints(2,:) ~= 0
    % If you need the track in single color uncomment the following
    % lines
    %         ht=plot(trackedPoints(2*(point-1)+2,:),trackedPoints(2*(point-1)+1,:),'g');
    %         set(ht,'LineWidth',1);
    
    % To denote the tracked points using different colors.
    % Divide the FRAME into six segments.
    color_seg = floor(totalNumberPhases/6);
    
    % Split the tracked points into 6 segments and color as
    % different color on the rainbow.
    
    ht1=plot(hImageAxes,trackedPoints(2,1:color_seg+1),trackedPoints(1,1:color_seg+1));
    set(ht1,'LineWidth',2,'Color',[255 0 0]./255);
    
    ht2=plot(hImageAxes,trackedPoints(2,color_seg+1:2*color_seg+1),trackedPoints(1,color_seg+1:2*color_seg+1));
    set(ht2,'LineWidth',2,'Color',[255 155 0]./255);
    
    ht3=plot(hImageAxes,trackedPoints(2,2*color_seg+1:3*color_seg+1),trackedPoints(1,2*color_seg+1:3*color_seg+1));
    set(ht3,'LineWidth',2,'Color',[255 255 0]./255);
    
    ht4=plot(hImageAxes,trackedPoints(2,3*color_seg+1:4*color_seg+1),trackedPoints(1,3*color_seg+1:4*color_seg+1));
    set(ht4,'LineWidth',2,'Color',[0 255 0]./255);
    
    ht5=plot(hImageAxes,trackedPoints(2,4*color_seg+1:5*color_seg+1),trackedPoints(1,4*color_seg+1:5*color_seg+1));
    set(ht5,'LineWidth',2,'Color',[0 0 255]./255);
    
    ht6=plot(hImageAxes,trackedPoints(2,5*color_seg+1:totalNumberPhases),trackedPoints(1,5*color_seg+1:totalNumberPhases));
    set(ht6,'LineWidth',2,'Color',[155 0 255]./255);
    
    ht1=plot(hImageAxes,trackedPoints(2,CurrentPhase),trackedPoints(1,CurrentPhase),'ro');
    set(ht1,'LineWidth',2);
end


% --- Executes during object creation, after setting all properties.
function LST_PRVWTYPE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LST_PRVWTYPE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function EDT_DATSET_NAME_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_DATSET_NAME (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function EDT_PATH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_PATH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function AXIS_IMG_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to AXIS_IMG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in BTN_ESTRAIN_2D.
function BTN_ESTRAIN_2D_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_ESTRAIN_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SliceList = getappdata(0,'SliceList');
if(isempty(SliceList))
    set(hObject,'Value',1);
    return;
end

SlicesWith2DResults = 0;

SliceListLength = length(SliceList);
hWaitBar1 = waitbar(0,['Calculating 2D Eulerian strain for ' num2str(SliceListLength) ' Slices. Please wait...']);

for SliceIndex=1:SliceListLength
    if ~strcmp(SliceList{SliceIndex}.ProcessInfo.FilterStatus.flag,'clean')
        errordlg('Please filter all slices first before calculating 2D Eulerian strain.');
        SliceList{SliceIndex}.ProcessInfo.Strain2DStatus.flag = 'dirty';
        setappdata(0,'SliceList',SliceList);
    else
        SlicesWith2DResults = SlicesWith2DResults+1;
        waitbar(SliceIndex/SliceListLength,hWaitBar1,['Calculating 2D Eulerian strain for slice ' num2str(SliceIndex) ' of ' num2str(SliceListLength) ' slices.'])
        
        EStrain2D_4(SliceIndex);
        SliceList = getappdata(0,'SliceList');
        
        SliceList{SliceIndex}.AvailResults    = ...
            AddToAvailableResults(SliceList{SliceIndex}.AvailResults,...
            [{'Ecc'};{'Err'};{'Exx'};{'Eyy'};{'CircAng1_2D'};{'CircAng2_2D'};{'Ep1_2D'};{'Ep2_2D'}]);
        
        SliceList{SliceIndex}.ProcessInfo.Strain2DStatus.flag = 'clean';
        SliceList{SliceIndex}.ProcessInfo.Strain2DStatus.time = datestr(clock,0);
        
        setappdata(0,'SliceList',SliceList);
    end
end
close(hWaitBar1);

hImportHARPGui  = getappdata(0,'hImportHARPGui');
hLST_PRVWTYPE   = findobj(hImportHARPGui,'Tag','LST_PRVWTYPE');
hLST_SLICES     = findobj(hImportHARPGui,'Tag','LST_SLICES');
CurrentSlice    = get(hLST_SLICES,'Value');
set(hLST_PRVWTYPE,'String',SliceList{CurrentSlice}.AvailResults);

if ~strcmp(SliceList{SliceIndex}.ProcessInfo.Strain2DStatus.flag,'clean') || ~isfield(SliceList{CurrentSlice}.Strain2DInfoSet,'EC_STRAIN')
    set(hLST_PRVWTYPE,'Value',FindStringIndex(SliceList{CurrentSlice}.AvailResults,'ROI Magnitude'));
else
    set(hLST_PRVWTYPE,'Value',FindStringIndex(SliceList{CurrentSlice}.AvailResults,'Ecc'));
end

UpdatePreviewGUI;
UpdatePreviewROIGUI;
UpdatePropertiesGUI;
Update_Colorbar;


% --- Executes during object creation, after setting all properties.
function EDT_OMEGA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_OMEGA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CHK_TWODYNAM.
function CHK_TWODYNAM_Callback(hObject, eventdata, handles)
% hObject    handle to CHK_TWODYNAM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CHK_TWODYNAM
