function varargout = ThreshSeg(varargin)
% THRESHSEG MATLAB code for ThreshSeg.fig
%      THRESHSEG, by itself, creates a new THRESHSEG or raises the existing
%      singleton*.
%
%      H = THRESHSEG returns the handle to a new THRESHSEG or the handle to
%      the existing singleton*.
%
%      THRESHSEG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THRESHSEG.M with the given input arguments.
%
%      THRESHSEG('Property','Value',...) creates a new THRESHSEG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ThreshSeg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ThreshSeg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ThreshSeg

% Last Modified by GUIDE v2.5 05-Jul-2016 08:33:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ThreshSeg_OpeningFcn, ...
                   'gui_OutputFcn',  @ThreshSeg_OutputFcn, ...
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


% --- Executes just before ThreshSeg is made visible.
function ThreshSeg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ThreshSeg (see VARARGIN)

% Add lib
addpath('./lib/');

% Show welcome image
handles.I_original = imread('./resources/startup.jpg');
handles.I = handles.I_original;
axes(handles.axes1);
imshow(handles.I);

% Choose default command line output for ThreshSeg
handles.output = hObject;

% Set default values for parameters
handles.normalization_method = Normalize.None;
handles.noise_level = 0;

handles.dt_original = 0.03;
handles.dt       = handles.dt_original;
handles.lambda   = 0.02;
handles.n_phases = 2;
handles.resize1  = 0;
handles.resize2  = 0;
handles.max_iter = 10000;
handles.damping  = 0;

handles.g = im2double(handles.I);
[handles.M, handles.N, handles.n_channels] = size(handles.g);

handles.u = ones(handles.M, handles.N, handles.n_phases);
handles.u_init = handles.u;
handles.energy = [];

handles.initialized = false;
handles.transient_fill = false;
handles.time_step_number = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ThreshSeg wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ThreshSeg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DrawInit.
function DrawInit_Callback(hObject, eventdata, handles)
% hObject    handle to DrawInit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

choice = questdlg(['Please draw ' int2str(handles.n_phases-1) ' region(s) on the image.'], ...
	'Choose Initials', ...
	'OK','Cancel','OK');
switch choice
    case 'OK'
        handles.initialized = true;
        axes(handles.axes1);
        imshow(handles.g);
        hold on;
        
        handles.u = zeros(handles.M,handles.N,handles.n_phases);
        
        for phase = 1:handles.n_phases-1
            rect = getrect; % rect = [xmin, ymin, width, height]
            % lmt = [ymin, ymax, xmin, xmax] due to the transpose of
            % image show.
            rect = floor(rect);
            lmt = [rect(2), rect(2) + rect(4), rect(1), rect(1) + rect(3)];
            lmt = max(lmt, 1);
            % Handle out-of-image selections
            lmt(1:2) = min(lmt(1:2), handles.M);
            lmt(3:4) = min(lmt(3:4), handles.N);
            handles.u(lmt(1):lmt(2), lmt(3):lmt(4), phase) = 1;
            contour(handles.u(:,:,phase), [0.5 0.5], 'r', 'LineWidth',2);
            drawnow;
        end
              
        % Handle intersection
        intersection_mask = ones(handles.M,handles.N) - (sum(handles.u,3) > 1);
        for phase = 1:handles.n_phases-1
            handles.u(:,:,phase) = min(handles.u(:,:,phase), intersection_mask);
        end;
        
        % Background phase
        handles.u(:,:,handles.n_phases) = 1 - sum(handles.u,3);
        
        handles.u_init = handles.u;
        
    case 'Cancel'
        return
end

guidata(handles.figure1, handles);



% --- Executes on button press in Help.
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

message = {['Load an image to get started. ' ...
           'To specify inital conditions, either draw rectangles ' ...
           'or load from files saved in "Polygon" format.']; ''; ...
           'Hint: set max_iter to be 1 to run step by step.'; ''; ...
           'Click on "More Info.." for more information.'; ''; ...
           };
header = 'ThreshSeg v0.2 (2016)';
    
react = questdlg(message, header, 'More Info..','Back','Back');

if strcmp(react,'More Info..')
    web('http://www.math.ust.hk/~mawang/');
end


% --- Executes on button press in LoadImage.
function LoadImage_Callback(hObject, eventdata, handles)
% hObject    handle to LoadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {'*.jpg', 'All jpg-Files (*.jpg)'; '*.bmp','All bmp-Files (*.bmp)' ; '*.png' , 'All png-Files (*.png)' ; '*.*' , 'All Files'}, 'Select Address Book');
% If "Cancel" is selected then return
if isequal([filename,pathname],[0,0])
    return
% Otherwise construct the fullfilename and Check and load the file
else
    File = fullfile(pathname,filename);
    handles.I_original = imread(File) ;
    handles.I = handles.I_original;
    axes(handles.axes1);
    cla;
    imshow(handles.I);
end

handles.g = im2double(handles.I);
[handles.M, handles.N, handles.n_channels] = size(handles.g);

guidata(handles.figure1, handles);




% --- Executes on button press in LoadInit.
function LoadInit_Callback(hObject, eventdata, handles)
% hObject    handle to LoadInit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile( ...
    {'*.txt', 'All txt-Files (*.txt)'; '*.*' , 'All Files'}, 'Select Address Book');
% If "Cancel" is selected then return
if isequal([filename,pathname],[0,0])
    return
% Otherwise construct the fullfilename and Check and load the file
else
    File = fullfile(pathname,filename); 
    A = dlmread(File);
                
    % Handle abnormal inputs
    if handles.n_phases ~= size(A,1)/2 + 1
        errordlg ('Invalid input file! (Number of shapes must equal to n_phases-1!)', 'Initials Error');
        return
    end
    A = max(A, 0);
    A = min(A, 1);
    
    handles.u = zeros(handles.M, handles.N, handles.n_phases);
    x = (1:handles.N) ./ handles.N;
    y = (1:handles.M) ./ handles.M;
    [xx,yy] = meshgrid(x, y);
    
    for phase = 1:handles.n_phases-1
        n_vertices = nnz( A(2*phase-1,:) + A(2*phase,:) );
        handles.u(:,:,phase) = inpolygon(xx,yy, ...
            A(2*phase-1,1:n_vertices),A(2*phase,1:n_vertices));
    end
    
    % Handle intersection
    intersection_mask = ones(handles.M,handles.N) - (sum(handles.u,3) > 1);
    for phase = 1:handles.n_phases-1
        handles.u(:,:,phase) = min(handles.u(:,:,phase), intersection_mask);
    end;
    
    % Background phase
    handles.u(:,:,handles.n_phases) = 1 - sum(handles.u,3);
    
    handles.u_init = handles.u;
    
    handles.initialzed = true;
end

axes(handles.axes1);
ColOrd   = get(gca,'ColorOrder');
n_colors = size(ColOrd, 1);
imshow(handles.g);
hold on;
for phase = 1:handles.n_phases-1
    Col = ColOrd(mod(phase,n_colors)+1,:);
    contour(handles.u(:,:,phase), [0.5 0.5], 'Color', Col, 'LineWidth',2); 
end

guidata(handles.figure1, handles);


% --- Executes on selection change in Normalize.
function Normalize_Callback(hObject, eventdata, handles)
% hObject    handle to Normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Normalize contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Normalize

contents = cellstr(get(hObject,'String'));
switch contents{get(hObject,'Value')}
    case 'None'
        handles.normalization_method = Normalize.None;
    case 'Max'
        handles.normalization_method = Normalize.Max;
    case 'MinMax'
        handles.normalization_method = Normalize.MinMax;
    case 'MeanVar'
        handles.normalization_method = Normalize.MeanVar;
end

handles.g = Normalize.apply (im2double(handles.I) + handles.noise_level .* randn(size(handles.g)), ...
                             handles.normalization_method);
axes(handles.axes1);
imshow(handles.g);

guidata(handles.figure1, handles);



% --- Executes during object creation, after setting all properties.
function Normalize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TimeStepSize_Callback(hObject, eventdata, handles)
% hObject    handle to TimeStepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeStepSize as text
%        str2double(get(hObject,'String')) returns contents of TimeStepSize as a double

handles.dt_original = str2double(get(hObject,'String'));
handles.dt = handles.dt_original;
guidata(handles.figure1, handles);

% --- Executes during object creation, after setting all properties.
function TimeStepSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeStepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lambda_Callback(hObject, eventdata, handles)
% hObject    handle to TimeStepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeStepSize as text
%        str2double(get(hObject,'String')) returns contents of TimeStepSize as a double

handles.lambda = str2double(get(hObject,'String'));
guidata(handles.figure1, handles);

% --- Executes during object creation, after setting all properties.
function Lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeStepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowPreprocessedImage.
function ShowPreprocessedImage_Callback(hObject, eventdata, handles)
% hObject    handle to ShowPreprocessedImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);
cla;
hold off;
imshow(handles.g);


% --- Executes on button press in ShowContour.
function ShowContour_Callback(hObject, eventdata, handles)
% hObject    handle to ShowContour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);
cla;
hold off;
ColOrd   = get(gca,'ColorOrder');
n_colors = size(ColOrd, 1);
imshow(handles.g);
hold on;
for phase = 1:handles.n_phases-1
    Col = ColOrd(mod(phase,n_colors)+1,:);
    contour(handles.u(:,:,phase), [0.5 0.5], 'Color', Col, 'LineWidth',2); 
end
hold off;

handles.transient_fill = false;
guidata(handles.figure1, handles);


% --- Executes on button press in ShowCartoonization.
function ShowCartoonization_Callback(hObject, eventdata, handles)
% hObject    handle to ShowCartoonization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);
cla;
hold off;
cartoon = zeros(handles.M,handles.N);
for phase = 1:handles.n_phases
    cartoon = cartoon + phase * handles.u(:,:,phase);
end
contourf(cartoon);
axis equal;
set(gca,'ydir','reverse');

% No box
axis off;
hold off;

handles.transient_fill = true;
guidata(handles.figure1, handles);


% --- Executes on button press in ShowEnergy.
function ShowEnergy_Callback(hObject, eventdata, handles)
% hObject    handle to ShowEnergy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);
cla;
hold off;

plot(handles.energy);
grid on;

% No box
axis off;
hold off;



% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.initialized
    choice = questdlg(['Please Specify initial conditions!'], ...
        'Missing Initials?', ...
        'Still Run','Cancel','Cancel');
    switch choice
        case 'Still Run'
            true;
        case 'Cancel'
            return
    end
end

iter = handles.time_step_number;
is_moving = true;
min_index = 1;
axes(handles.axes1);
ColOrd   = get(gca,'ColorOrder');
n_colors = size(ColOrd, 1);
while get(hObject, 'Value')
    iter = iter + 1;
    set(hObject, 'String', ['Running ' int2str(iter)]);
    
 %% Main Iteration
 
    % Update data term
    handles.f = compute_data_term (handles.g, handles.u);
    
    % Apply heat kernel convolution
    handles.uh = apply_heat_convolution (handles.dt, handles.u);
    
    % Thresholding
    is_moving = false;
    for j=1:handles.N
        for i=1:handles.M
            % Get scores for each phase
            min_score = Inf;
            for phase = 1:handles.n_phases
                score = handles.f(i,j,phase) - 2 * handles.lambda / sqrt(handles.dt) * sqrt(pi) * handles.uh(i,j,phase);
                if score < min_score
                    min_score = score;
                    min_index = phase;
                end;
            end
            if handles.u(i,j,min_index)==0
                is_moving = true;
            end
            handles.u(i,j,:) = 0;
            handles.u(i,j,min_index) = 1;
        end
    end
    
    % Damping
    handles.dt = handles.dt * (1-handles.damping);
    set(handles.TimeStepSize,   'String', num2str(handles.dt));
    
    % Update energy    
    handles.energy(iter) = compute_energy (handles.u, handles.uh,...
                                handles.f, handles.dt, handles.lambda);
    
    if iter >= handles.max_iter
        set(hObject, 'Value', false);
        % handles.time_step_number = iter - 1;
        set(hObject, 'String', 'Run More');
    elseif ~is_moving
        set(hObject, 'Value', false);
        set(hObject, 'String', ['Done: ' int2str(iter)]);
    end
    
    cla;
    hold off;
    if handles.transient_fill
        % Update cartoonization
        cartoon = zeros(handles.M,handles.N);
        for phase = 1:handles.n_phases
            cartoon = cartoon + phase * handles.u(:,:,phase);
        end
        contourf(cartoon);
        axis equal;
        set(gca,'ydir','reverse');
        axis off;
    else
        % Update contour
        imshow(handles.g);
        hold on;
        for phase = 1:handles.n_phases-1
            Col = ColOrd(mod(phase,n_colors)+1,:);
            contour(handles.u(:,:,phase), [0.5 0.5], 'Color', Col, 'LineWidth',2);
        end
    end
    hold off;
    
    drawnow;
end

if is_moving && ~get(hObject, 'Value')
    set(hObject, 'String', ['Paused ' int2str(iter)]);
    handles.time_step_number = max(iter - 1,0);
end

guidata(handles.figure1, handles);


% --- Executes on button press in Restart.
function Restart_Callback(hObject, eventdata, handles)
% hObject    handle to Restart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

choice = questdlg(['Really restart? Contours will be reset to initials.'], ...
	'Confirm Restart', ...
	'Restart','Cancel','Restart');
switch choice
    case 'Restart'
        handles.u = handles.u_init;
        handles.time_step_number = 0;
        handles.dt = handles.dt_original;
        handles.energy = [];
        set(handles.Run, 'String', 'Run');
    case 'Cancel'
        return
end

axes(handles.axes1);
cla;
hold off;
ColOrd   = get(gca,'ColorOrder');
n_colors = size(ColOrd, 1);
imshow(handles.g);
hold on;
for phase = 1:handles.n_phases-1
    Col = ColOrd(mod(phase,n_colors)+1,:);
    contour(handles.u(:,:,phase), [0.5 0.5], 'Color', Col, 'LineWidth',2); 
end
hold off;

guidata(handles.figure1, handles);


function NPhases_Callback(hObject, eventdata, handles)
% hObject    handle to NPhases (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NPhases as text
%        str2double(get(hObject,'String')) returns contents of NPhases as a double

handles.n_phases = str2double(get(hObject,'String'));
guidata(handles.figure1, handles);


% --- Executes during object creation, after setting all properties.
function NPhases_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NPhases (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowImage.
% Display image
function ShowImage_Callback(hObject, eventdata, handles)
% hObject    handle to ShowImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
imshow(handles.I);



function Resize1_Callback(hObject, eventdata, handles)
% hObject    handle to Resize1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Resize1 as text
%        str2double(get(hObject,'String')) returns contents of Resize1 as a double

handles.resize1 = floor(str2double(get(hObject,'String')));

if handles.resize1~=0 && handles.resize2~=0
    handles.I = imresize(handles.I_original, [handles.resize1, handles.resize2]);
else
    handles.I = handles.I_original;
end
[handles.M, handles.N, ~] = size(handles.I);
handles.g = Normalize.apply (im2double(handles.I) + handles.noise_level .* randn(size(handles.I)), ...
                             handles.normalization_method);
axes(handles.axes1);
imshow(handles.g);

guidata(handles.figure1, handles);


% --- Executes during object creation, after setting all properties.
function Resize1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Resize1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Resize2_Callback(hObject, eventdata, handles)
% hObject    handle to Resize2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Resize2 as text
%        str2double(get(hObject,'String')) returns contents of Resize2 as a double
handles.resize2 = floor(str2double(get(hObject,'String')));

if handles.resize1~=0 && handles.resize2~=0
    handles.I = imresize(handles.I_original, [handles.resize1, handles.resize2]);
else
    handles.I = handles.I_original;
end
[handles.M, handles.N, ~] = size(handles.I);
handles.g = Normalize.apply (im2double(handles.I) + handles.noise_level .* randn(size(handles.I)), ...
                             handles.normalization_method);
axes(handles.axes1);
imshow(handles.g);

guidata(handles.figure1, handles);


% --- Executes during object creation, after setting all properties.
function Resize2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Resize2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxIter_Callback(hObject, eventdata, handles)
% hObject    handle to MaxIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxIter as text
%        str2double(get(hObject,'String')) returns contents of MaxIter as a double

handles.max_iter = floor(str2double(get(hObject,'String')));
guidata(handles.figure1, handles);



% --- Executes during object creation, after setting all properties.
function MaxIter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DampingRate_Callback(hObject, eventdata, handles)
% hObject    handle to DampingRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DampingRate as text
%        str2double(get(hObject,'String')) returns contents of DampingRate as a double

handles.damping = str2double(get(hObject,'String'));
guidata(handles.figure1, handles);


% --- Executes during object creation, after setting all properties.
function DampingRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DampingRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NoiseLevel_Callback(hObject, eventdata, handles)
% hObject    handle to NoiseLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NoiseLevel as text
%        str2double(get(hObject,'String')) returns contents of NoiseLevel as a double

handles.noise_level = str2double(get(hObject,'String'));

handles.g = Normalize.apply (im2double(handles.I) + handles.noise_level .* randn(size(handles.g)), ...
                             handles.normalization_method);
axes(handles.axes1);
imshow(handles.g);

guidata(handles.figure1, handles);


% --- Executes during object creation, after setting all properties.
function NoiseLevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NoiseLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
