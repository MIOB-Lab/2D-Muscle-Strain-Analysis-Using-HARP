function varargout = visualizeHARP(varargin)
% VISUALIZEHARP MATLAB code for visualizeHARP.fig
%      VISUALIZEHARP, by itself, creates a new VISUALIZEHARP or raises the existing
%      singleton*.
%
%      H = VISUALIZEHARP returns the handle to a new VISUALIZEHARP or the handle to
%      the existing singleton*.
%
%      VISUALIZEHARP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUALIZEHARP.M with the given input arguments.
%
%      VISUALIZEHARP('Property','Value',...) creates a new VISUALIZEHARP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before visualizeHARP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to visualizeHARP_OpeningFcn via varargin.
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
% Copyright (c) 1999-2012 Johns Hopkins University, Image Analysis
% and Communications Lab (IACL)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @visualizeHARP_OpeningFcn, ...
                   'gui_OutputFcn',  @visualizeHARP_OutputFcn, ...
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


% --- Executes during object creation, after setting all properties.
function Visualize_HARP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Visualize_HARP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
setappdata(0,'hVisualizeHARPGui',hObject);


% --- Executes just before visualizeHARP is made visible.
function visualizeHARP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to visualizeHARP (see VARARGIN)

% Choose default command line output for visualizeHARP
handles.output = hObject;

% Update handles structure
guidata(hObject,handles);
set(handles.Visualize_HARP,'resize','on');

% Make the initial plot.
UpdatePreviewGUI(handles.AXES_CURRENT_FIG);


function Visualize_HARP_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to Visualize_HARP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hVisualizeHARPGui = getappdata(0,'hVisualizeHARPGui');
delete(hVisualizeHARPGui);


% --- Outputs from this function are returned to the command line.
function varargout = visualizeHARP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function Current_Frame_Slider_Callback(hObject, eventdata, handles)
% hObject    handle to Current_Frame_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

hVisualizeHARPGui = getappdata(0,'hVisualizeHARPGui');
hCurrent_Frame_Slider = findobj(hVisualizeHARPGui,'Tag','Current_Frame_Slider');
hCurrent_Frame_Editbox = findobj(hVisualizeHARPGui,'Tag','Current_Frame_Editbox');

% Set the slider value to the edit box.
set(hCurrent_Frame_Editbox,'String',num2str(round(get(hCurrent_Frame_Slider,'Value'))));

UpdatePreviewGUI(handles.AXES_CURRENT_FIG);


% --- Executes during object creation, after setting all properties.
function Current_Frame_Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Current_Frame_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% Set the initial parameters of the slider.
SliceList = getappdata(0,'SliceList');

% Set the minimum, maximum, slider step and the initial value.
if (SliceList{1}.NO_FRAMES > 1)
    NO_FRAMES = SliceList{1}.NO_FRAMES;
    set(hObject,'Visible','on');
    set(hObject,'Min',1);
    set(hObject,'Max',NO_FRAMES);
    set(hObject,'SliderStep',([1/(NO_FRAMES-1) 1/(NO_FRAMES-1)]));
    set(hObject,'Value',1);
else
    set(hObject,'Visible','off');
end


function Current_Frame_Editbox_Callback(hObject, eventdata, handles)
% hObject    handle to Current_Frame_Editbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Current_Frame_Editbox as text
%        str2double(get(hObject,'String')) returns contents of Current_Frame_Editbox as a double
hVisualizeHARPGui = getappdata(0,'hVisualizeHARPGui');

% Set the edit box value to the slider.
hCurrent_Frame_Editbox = findobj(hVisualizeHARPGui,'Tag','Current_Frame_Editbox');
hCurrent_Frame_Slider  = findobj(hVisualizeHARPGui,'Tag','Current_Frame_Slider');

if str2double(get(hCurrent_Frame_Editbox,'String')) > get(hCurrent_Frame_Slider,'Max') || ...
        str2double(get(hCurrent_Frame_Editbox,'String')) < get(hCurrent_Frame_Slider,'Min')
    errordlg('Invalid entry of time frame. It should be within 1 and the total number of frames.');
    set(hCurrent_Frame_Editbox,'String',num2str(get(hCurrent_Frame_Slider,'Value')));
    return;
else
    set(hCurrent_Frame_Slider,'Value',round(str2double(get(hCurrent_Frame_Editbox,'String'))));
    set(hCurrent_Frame_Editbox,'String',round(str2double(get(hCurrent_Frame_Editbox,'String'))));
end

UpdatePreviewGUI(handles.AXES_CURRENT_FIG);


% --- Executes during object creation, after setting all properties.
function Current_Frame_Editbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Current_Frame_Editbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in LST_SLICES.
function LST_SLICES_Callback(hObject, eventdata, handles)
% hObject    handle to LST_SLICES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LST_SLICES contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LST_SLICES
SliceList  = getappdata(0, 'SliceList');
NewSlice = get(hObject,'Value'); % current slice

if (SliceList{NewSlice}.NO_FRAMES == 1)
    set(hObject,'Value',NewSlice);
    set(handles.Current_Frame_Editbox,'String','1');
    set(handles.Current_Frame_Slider,'Value',1);
end

if ~isfield(SliceList{NewSlice},'HarpMagnitude')
    errordlg('This slice was not processed by HARP!');
end

UpdatePreviewGUI(handles.AXES_CURRENT_FIG);


% --- Executes during object creation, after setting all properties.
function LST_SLICES_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LST_SLICES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% The list of slices should be in string of the 'Slice' List box in
% importHARPGui. Go to hImportHARPGui, search for the Listbox containing 
% the Slices, and then get the string and fill this popupmenu.

hImportHARPGui = getappdata(0,'hImportHARPGui');
hLST_SLICES = findobj(hImportHARPGui,'Tag','LST_SLICES');
set(hObject,'String',get(hLST_SLICES,'String'));
set(hObject,'Value',1);


% --- Executes on button press in BTN_EXPRT_IMG.
function BTN_EXPRT_IMG_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_EXPRT_IMG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = figure;
UpdatePreviewGUI(gca(h));

% Make a title.
hVisualizeHARPGui = getappdata(0,'hVisualizeHARPGui');
% Get current slice and current phase.
hLST_SLICES = findobj(hVisualizeHARPGui,'Tag','LST_SLICES');
CurrentSlice = get(hLST_SLICES,'Value');
hCurrent_Frame_Editbox = findobj(hVisualizeHARPGui, 'Tag', 'Current_Frame_Editbox');
CurrentFrame = str2double(get(hCurrent_Frame_Editbox,'String'));

titleHandle = title(gca,['Slice ' num2str(CurrentSlice) ' Time Frame ' num2str(CurrentFrame)]);
set(titleHandle ,'FontWeight','Bold','FontSize',11);

% This part for future purpose to fix save sequence. FX
% for CurrentFrame = 1:26
% set(handles.Current_Frame_Slider,'Value',CurrentFrame);
% set(handles.Current_Frame_Editbox,'String',num2str(CurrentFrame));
% saveas(gcf,[DataInfo ' Image - Time Frame ' num2str(CurrentFrame)],'png');
% close(11);

% --- Executes on button press in BTN_EXPRT_IMG_All.
function BTN_EXPRT_IMG_All_Callback(hObject, eventdata, handles)
% hObject    handle to BTN_EXPRT_IMG_All (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hImportHARPGui = getappdata(0,'hImportHARPGui');
SliceList = getappdata(0,'SliceList');

CurrentSlice = get(handles.LST_SLICES,'Value');
NoPhases = SliceList{CurrentSlice}.NO_FRAMES;

% Generate default filename and pathname.
hEDT_PATH = findobj(hImportHARPGui,'Tag','EDT_PATH');
pathname = get(hEDT_PATH,'string');
if isfield(SliceList{CurrentSlice},'DatasetInfo')
    % Check if there is a date extension.
    % date format: 'yyyy-mm-dd'
    filename = SliceList{1}.DatasetInfo.FileName;
    [match start] = regexpi(filename,'-(\d{4})-(\d{2})-(\d{2})','match','start');
    if ~isempty(match)
        filename = filename(1:start-1);
    end
end
if isempty(pathname)
    pathname = pwd; % current directory
end
if isempty(filename)
    filename = 'New';
end
filename = [filename '-Slc' num2str(CurrentSlice) '.png'];

% Use define filename and pathname.
[filename, pathname, dummy]    = uiputfile({'*.png','PNG files (*.png)';'*.jpg','JPEG images (*.jpg)';...
    '*.gif','GIF images (*.gif)'},'Save Images',[pathname filesep filename]);
if isequal(filename,0) || isequal(pathname,0) % cancel pressed
    return;
end

WaitBar = waitbar(0,'Saving Images of All Time Frames...');
   
for my_index = 1:NoPhases
    set(handles.Current_Frame_Slider,'Value',my_index);
    set(handles.Current_Frame_Editbox,'String',num2str(my_index));
    
    % Plot in a new figure.
    h = figure;
    UpdatePreviewGUI(gca(h));
  
    % Capture the figure frame.
    ClrMap = frame2im(getframe(gca(h)));
    [indImage,Map] = rgb2ind(ClrMap,1024);
    
    % Write and save the image.
    fullName = [pathname filename(1:end-4) '-Tf' num2str(my_index) filename(end-3:end)];
    imwrite(indImage,Map,fullName);
    
    % Update waitbar and close figure.
    close(h);    
    waitbar(my_index/NoPhases,WaitBar);
end

close(WaitBar);


function EDT_VEL_VEC_VSTRETCH_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_VEL_VEC_VSTRETCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_VEL_VEC_VSTRETCH as text
%        str2double(get(hObject,'String')) returns contents of EDT_VEL_VEC_VSTRETCH as a double
UpdatePreviewGUI(handles.AXES_CURRENT_FIG);


% --- Executes during object creation, after setting all properties.
function EDT_VEL_VEC_VSTRETCH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_VEL_VEC_VSTRETCH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function EDT_VEL_VEC_VDS_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_VEL_VEC_VDS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_VEL_VEC_VDS as text
%        str2double(get(hObject,'String')) returns contents of EDT_VEL_VEC_VDS as a double
UpdatePreviewGUI(handles.AXES_CURRENT_FIG);


% --- Executes during object creation, after setting all properties.
function EDT_VEL_VEC_VDS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_VEL_VEC_VDS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function EDT_VEL_VEC_VTHRESH_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_VEL_VEC_VTHRESH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_VEL_VEC_VTHRESH as text
%        str2double(get(hObject,'String')) returns contents of EDT_VEL_VEC_VTHRESH as a double
UpdatePreviewGUI(handles.AXES_CURRENT_FIG);


% --- Executes during object creation, after setting all properties.
function EDT_VEL_VEC_VTHRESH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_VEL_VEC_VTHRESH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function SLD_MSK_THRE1_Callback(hObject, eventdata, handles)
% hObject    handle to SLD_MSK_THRE1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

hVisualizeHARPGui = getappdata(0,'hVisualizeHARPGui');

maskThre = get(hObject,'Value') - mod(get(hObject,'Value'),0.01);

% Update window display
hMask_Threshold_Editbox = findobj(hVisualizeHARPGui,'Tag','EDT_MSK_THRE1');
set(hMask_Threshold_Editbox,'String',num2str(maskThre));
% Do the 'mod' operation to restrict the threshold to two decimal places

UpdatePreviewGUI(handles.AXES_CURRENT_FIG);


% --- Executes during object creation, after setting all properties.
function SLD_MSK_THRE1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SLD_MSK_THRE1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function EDT_MSK_THRE1_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_MSK_THRE1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_MSK_THRE1 as text
%        str2double(get(hObject,'String')) returns contents of EDT_MSK_THRE1 as a double

hVisualizeHARPGui = getappdata(0,'hVisualizeHARPGui');

maskThre = round(str2double(get(hObject,'string'))*100)/100;

% update window display
hMask_Threshold_Slider = findobj(hVisualizeHARPGui,'Tag','SLD_MSK_THRE1');
set(hMask_Threshold_Slider,'value',maskThre);
% Do this 'mod' thing to restrict the threshold to two decimal places

UpdatePreviewGUI(handles.AXES_CURRENT_FIG);


% --- Executes during object creation, after setting all properties.
function EDT_MSK_THRE1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_MSK_THRE1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function SLD_MSK_THRE2_Callback(hObject, eventdata, handles)
% hObject    handle to SLD_MSK_THRE2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

hVisualizeHARPGui = getappdata(0,'hVisualizeHARPGui');

maskThre = get(hObject,'Value') - mod(get(hObject,'Value'),0.01);

% Update window display
hMask_Threshold_Editbox = findobj(hVisualizeHARPGui,'Tag','EDT_MSK_THRE2');
set(hMask_Threshold_Editbox,'String',num2str(maskThre));
% Do the 'mod' operation to restrict the threshold to two decimal places

UpdatePreviewGUI(handles.AXES_CURRENT_FIG);


% --- Executes during object creation, after setting all properties.
function SLD_MSK_THRE2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SLD_MSK_THRE2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function EDT_MSK_THRE2_Callback(hObject, eventdata, handles)
% hObject    handle to EDT_MSK_THRE2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EDT_MSK_THRE2 as text
%        str2double(get(hObject,'String')) returns contents of EDT_MSK_THRE2 as a double

hVisualizeHARPGui = getappdata(0,'hVisualizeHARPGui');

maskThre = round(str2double(get(hObject,'string'))*100)/100;

% update window display
hMask_Threshold_Slider = findobj(hVisualizeHARPGui,'Tag','SLD_MSK_THRE2');
set(hMask_Threshold_Slider,'value',maskThre);
% Do this 'mod' thing to restrict the threshold to two decimal places

UpdatePreviewGUI(handles.AXES_CURRENT_FIG);


% --- Executes during object creation, after setting all properties.
function EDT_MSK_THRE2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EDT_MSK_THRE2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function UpdatePreviewGUI(Haxis)
% This function is reserved for future expansion.
VelocityVectorsPlot(Haxis);


function VelocityVectorsPlot(Haxis)

% Clear current figure to avoid figure overlay.
axes(Haxis);
cla;

SliceList = getappdata(0,'SliceList');
hVisualizeHARPGui = getappdata(0,'hVisualizeHARPGui');

% Get current phase, slice, totalframe and NXY.
hCurrent_Frame_Editbox = findobj(hVisualizeHARPGui,'Tag','Current_Frame_Editbox');
CurrentPhase = str2double(get(hCurrent_Frame_Editbox,'String'));
hLST_SLICES = findobj(hVisualizeHARPGui,'Tag','LST_SLICES');
CurrentSlice = get(hLST_SLICES,'Value');
TotalFrames = SliceList{CurrentSlice}.NO_FRAMES;

% Get stretch, downsample ratio and length threshold.
Stretch = str2double(get(findobj(hVisualizeHARPGui,'Tag','EDT_VEL_VEC_VSTRETCH'),'String'));
DSratio = str2double(get(findobj(hVisualizeHARPGui,'Tag','EDT_VEL_VEC_VDS'),'String'));
Thresh = str2double(get(findobj(hVisualizeHARPGui,'Tag','EDT_VEL_VEC_VTHRESH'),'String'));

% Create mask.
hSLD_MSK_THRE1 = findobj(hVisualizeHARPGui,'Tag','SLD_MSK_THRE1');
hSLD_MSK_THRE2 = findobj(hVisualizeHARPGui,'Tag','SLD_MSK_THRE2');
Thresh1 = get(hSLD_MSK_THRE1,'Value');
Thresh2 = get(hSLD_MSK_THRE2,'Value');
MaskThreshold = Thresh1+(CurrentPhase-1)/(TotalFrames-1)*(Thresh2-Thresh1);
HarpMagnitude = squeeze(SliceList{CurrentSlice}.HarpMagnitude(CurrentPhase,:,:));
HarpMagnitudeNormalized = HarpMagnitude./max(max(HarpMagnitude));
Mask = (HarpMagnitudeNormalized > MaskThreshold);

% Detect track method.
hRAD_SHOW_TRK_HARP = findobj(hVisualizeHARPGui,'Tag','RAD_SHOW_TRK_HARP');
hRAD_SHOW_TRK_REFINE = findobj(hVisualizeHARPGui,'Tag','RAD_SHOW_TRK_REFINE');
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

% Plot the HARP magnitue image.
UpdateImageShow(Haxis);
hold on;

% Get the velocity fields. NOTE: These are actually the displacement fields
% between two consecutive time frames. They do not have a /time factor comparing to
% true physical "velocity" or "speed". This is only for visualization as it
% is going to be stretched. The physical velocity value should be
% considered more carefully if needed.
[X,Y] = meshgrid(1:SliceList{1}.NXY(2),1:SliceList{1}.NXY(1));

% Select vector field.
hPNL_SHOW_TRACKS = findobj(hVisualizeHARPGui,'Tag','PNL_SHOW_TRACKS');
hObject = get(hPNL_SHOW_TRACKS,'SelectedObject');
Selection = get(hObject,'String');
if strcmp(Selection,'HARP')
    SliceList{CurrentSlice}.HarpTrackInfoSet.TrackForwardData(end+1,:,:,:) = ...
        SliceList{CurrentSlice}.HarpTrackInfoSet.TrackForwardData(end,:,:,:);
    DY = squeeze(SliceList{CurrentSlice}.HarpTrackInfoSet.TrackForwardData(CurrentPhase+1,:,:,1) - ...
        SliceList{CurrentSlice}.HarpTrackInfoSet.TrackForwardData(CurrentPhase,:,:,1));
    DX = squeeze(SliceList{CurrentSlice}.HarpTrackInfoSet.TrackForwardData(CurrentPhase+1,:,:,2) - ...
        SliceList{CurrentSlice}.HarpTrackInfoSet.TrackForwardData(CurrentPhase,:,:,2));
elseif strcmp(Selection,'HARP Refinement')
    DY = squeeze(SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardX(CurrentPhase,:,:));
    DX = squeeze(SliceList{CurrentSlice}.RefineTrackInfoSet.MotionForwardY(CurrentPhase,:,:));
else
    return;
end

% Downsample the velocity field.
if DSratio ~= 1
    DownSampMask = ones(size(DX));
    for i = 1:size(DX,1)
        if rem(i,DSratio) ~= 1
            DownSampMask(i,:) = zeros(size(DownSampMask(i,:)));
        end
    end
    for i = 1:size(DX,2)
        if rem(i,DSratio) ~= 1
            DownSampMask(:,i) = zeros(size(DownSampMask(:,i)));
        end
    end
    DX = DX.*DownSampMask;
    DY = DY.*DownSampMask;
end

% Preserve only velocities whose magnitudes are less than the specified threshold.
DX = (sqrt(DX.^2+DY.^2) < Thresh) .* DX;
DY = (sqrt(DX.^2+DY.^2) < Thresh) .* DY;

% Apply the mask.
X(~Mask) = NaN;
Y(~Mask) = NaN;

% Plot the velocity vectors.
axes(Haxis);
s_pos = quiver_harp(X,Y,Stretch*DX,Stretch*DY,0);
set(s_pos,'Color','red','linewidth',1.8);
set(gca,'YDir','reverse');
axis image; axis off;


% --- Executes on slider movement.
function GrayScale_Intensity_Slider_Callback(hObject, eventdata, handles)
% hObject    handle to GrayScale_Intensity_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

hVisualizeHARPGui = getappdata(0,'hVisualizeHARPGui');
hSlider = findobj(hVisualizeHARPGui,'Tag','GrayScale_Intensity_Slider');
hEditbox = findobj(hVisualizeHARPGui,'Tag','GrayScale_Intensity_Editbox');

% Set the slider value to the edit box.
set(hEditbox,'String',num2str(round(get(hSlider,'Value'))));

UpdatePreviewGUI(handles.AXES_CURRENT_FIG);


% --- Executes during object creation, after setting all properties.
function GrayScale_Intensity_Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GrayScale_Intensity_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function GrayScale_Intensity_Editbox_Callback(hObject, eventdata, handles)
% hObject    handle to GrayScale_Intensity_Editbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GrayScale_Intensity_Editbox as text
%        str2double(get(hObject,'String')) returns contents of GrayScale_Intensity_Editbox as a double

hVisualizeHARPGui = getappdata(0,'hVisualizeHARPGui');

% Set the edit box value to the slider.
hEditbox = findobj(hVisualizeHARPGui,'Tag','GrayScale_Intensity_Editbox');
hSlider  = findobj(hVisualizeHARPGui,'Tag','GrayScale_Intensity_Slider');

if str2double(get(hEditbox,'String')) > 100 || str2double(get(hEditbox,'String')) < 0
    errordlg('Invalid entry. The percentage should be within 0 and 100.');
    set(hEditbox,'String',num2str(get(hSlider,'Value')));
    return;
else
    set(hSlider,'Value',round(str2double(get(hEditbox,'String'))));
    set(hEditbox,'String',round(str2double(get(hEditbox,'String'))));
end

UpdatePreviewGUI(handles.AXES_CURRENT_FIG);


% --- Executes during object creation, after setting all properties.
function GrayScale_Intensity_Editbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GrayScale_Intensity_Editbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function UpdateImageShow(Haxis)

hVisualizeHARPGui = getappdata(0,'hVisualizeHARPGui');

% Get current slice and current phase.
hLST_SLICES = findobj(hVisualizeHARPGui,'Tag','LST_SLICES');
CurrentSlice = get(hLST_SLICES,'Value');
hCurrent_Frame_Editbox = findobj(hVisualizeHARPGui, 'Tag', 'Current_Frame_Editbox');
CurrentPhase = str2double(get(hCurrent_Frame_Editbox,'String'));

SliceList = getappdata(0,'SliceList');
totalNumberPhases = SliceList{CurrentSlice}.NO_FRAMES;

% Get mask threshholds.
hSLD_MSK_THRE1 = findobj(hVisualizeHARPGui,'Tag','SLD_MSK_THRE1');
hSLD_MSK_THRE2 = findobj(hVisualizeHARPGui,'Tag','SLD_MSK_THRE2');
Thresh1 = get(hSLD_MSK_THRE1,'Value');
Thresh2 = get(hSLD_MSK_THRE2,'Value');

% Get gray value intensity.
hGrayScale_Intensity_Editbox = findobj(hVisualizeHARPGui,'Tag','GrayScale_Intensity_Editbox');
IntensityThreshold = str2double(get(hGrayScale_Intensity_Editbox,'String'));

% Generate the mask from current frame threshold.
MaskThresholdValue = Thresh1+(CurrentPhase-1)/(totalNumberPhases-1)*(Thresh2-Thresh1);
HarpMagnitude = squeeze(SliceList{CurrentSlice}.HarpMagnitude(CurrentPhase,:,:));
HarpMagnitudeNormalized = HarpMagnitude./max(max(HarpMagnitude));
Mask = HarpMagnitudeNormalized > MaskThresholdValue;

% Get gray scale image.
[Image,ImClrMap,dummy] = ChooseDataType(CurrentSlice);
MyImage = squeeze(Image(CurrentPhase,:,:));

% Modify intensity threshold.
ImRange = [min(MyImage(:)) max(MyImage(:))];
if IntensityThreshold == 0
    IntensityThreshold = 0.001; % Prevent 0 intensity range errors.
end
ImRange = ImRange*(IntensityThreshold/100);

% Generate final gray scale image with proper mask.
MyImage = MyImage.*Mask;

% Display image.
axes(Haxis);
imagesc(MyImage);
colormap gray;
set(gca,'Clim',ImRange);
axis image; axis off;


function [MyImage,ImClrMap,DataInfo] = ChooseDataType(CurrentSlice)
% This function is reserved for future expansion use.
SliceList = getappdata(0, 'SliceList');

MyImage = SliceList{CurrentSlice}.HarpMagnitude; % image
ImClrMap= gray(1024); % colormap
DataInfo = 'HarpMagnitude'; % data type


% --- Executes when selected object is changed in PNL_SHOW_TRACKS.
function PNL_SHOW_TRACKS_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in PNL_SHOW_TRACKS 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

UpdatePreviewGUI(handles.AXES_CURRENT_FIG);
