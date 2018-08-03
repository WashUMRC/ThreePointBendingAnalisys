function varargout = ThreePointBendingAnalysisGUI(varargin)
% THREEPOINTBENDINGANALYSISGUI MATLAB code for ThreePointBendingAnalysisGUI.fig
%      THREEPOINTBENDINGANALYSISGUI, by itself, creates a new THREEPOINTBENDINGANALYSISGUI or raises the existing
%      singleton*.
%
%      H = THREEPOINTBENDINGANALYSISGUI returns the handle to a new THREEPOINTBENDINGANALYSISGUI or the handle to
%      the existing singleton*.
%
%      THREEPOINTBENDINGANALYSISGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THREEPOINTBENDINGANALYSISGUI.M with the given input arguments.
%
%      THREEPOINTBENDINGANALYSISGUI('Property','Value',...) creates a new THREEPOINTBENDINGANALYSISGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ThreePointBendingAnalysisGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ThreePointBendingAnalysisGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ThreePointBendingAnalysisGUI

% Last Modified by GUIDE v2.5 23-May-2017 16:04:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ThreePointBendingAnalysisGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ThreePointBendingAnalysisGUI_OutputFcn, ...
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


% --- Executes just before ThreePointBendingAnalysisGUI is made visible.
function ThreePointBendingAnalysisGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ThreePointBendingAnalysisGUI (see VARARGIN)

% Choose default command line output for ThreePointBendingAnalysisGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

handles.beamWidth = 7;

handles.rig = 'Dynamight';
% handles.pathstrResults = pwd;

guidata(hObject, handles);
% UIWAIT makes ThreePointBendingAnalysisGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ThreePointBendingAnalysisGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editBeamWidth_Callback(hObject, eventdata, handles)
% hObject    handle to editBeamWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBeamWidth as text
%        str2double(get(hObject,'String')) returns contents of editBeamWidth as a double
handles.beamWidth = str2num(get(handles.editBeamWidth,'String'));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editBeamWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBeamWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonLoadRigFiles.
function pushbuttonLoadRigFiles_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadRigFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get directory where Rig files are located (should all be in one
%level)
handles.pathstrRig = uigetdir(pwd,'Please select the directory containing your Rig force/displacement files');

%get load order option and generate list according to that order
str = get(handles.popupmenuRigLoadOptions,'String');
val = get(handles.popupmenuRigLoadOptions,'Value');
switch str{val}
    case 'Filename'
        handles.loadOptionsRig = 'Filename';
    case 'CreationDate'
        handles.loadOptionsRig = 'CreationDate';
    case 'Size'
        handles.loadOptionsRig = 'Size';
end

if strcmp(handles.rig,'Dynamight') == 1
    handles.filesRig = dir([handles.pathstrRig '\*.lvm']);
    if isempty(handles.filesRig)
        handles.filesRig = dir([handles.pathstrRig '\*.txt']);
    end
elseif strcmp(handles.rig,'Electropulse') == 1
    handles.filesRig = dir([handles.pathstrRig '\*.csv']);
elseif strcmp(handles.rig,'5866') == 1
    handles.filesRig = dir([handles.pathstrRig '\*.txt']);
end
    
if strcmpi(handles.loadOptionsRig,'CreationDate') == 1
    for i = 1:length(handles.filesRig)
        dNum(i) = handles.filesRig(i).datenum;
    end
    [~,idx] = sort(dNum,'ascend');
    handles.filesRig = handles.filesRig(idx);
elseif strcmpi(handles.loadOptionsRig,'Size') == 1
    for i = 1:length(handles.filesRig)
        dNum(i) = handles.filesRig(i).fileSize;
    end
    [~,idx] = sort(dNum,'ascend');
    handles.filesRig = handles.filesRig(idx);
elseif strcmpi(handles.loadOptionsRig,'Filename') == 1
    for i = 1:length(handles.filesRig)
        dNum{i} = handles.filesRig(i).name;
    end
    [~,idx] = sort(dNum);
    handles.filesRig = handles.filesRig(idx);
end

%update text field with list of Rig files in order

for i = 1:length(handles.filesRig)
    dFiles{i} = handles.filesRig(i).name;
end
set(handles.textRigFileList,'String',dFiles);

guidata(hObject, handles);

% --- Executes on selection change in popupmenuRigLoadOptions.
function popupmenuRigLoadOptions_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuRigLoadOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuRigLoadOptions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuRigLoadOptions


% --- Executes during object creation, after setting all properties.
function popupmenuRigLoadOptions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuRigLoadOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonLoadDICOMFiles.
function pushbuttonLoadDICOMFiles_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadDICOMFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.pathstrDICOM = uigetdir(pwd,'Please select the directory containing your DICOM force/displacement files');

%get load order option and generate list according to that order
str = get(handles.popupmenuDICOMLoadOptions,'String');
val = get(handles.popupmenuDICOMLoadOptions,'Value');
switch str{val}
    case 'Filename'
        handles.loadOptionsDICOM = 'Filename';
    case 'CreationDate'
        handles.loadOptionsDICOM = 'CreationDate';
    case 'Size'
        handles.loadOptionsDICOM = 'Size';
end
handles.filesDICOM = dir([handles.pathstrDICOM '\*.']);
handles.filesDICOM = handles.filesDICOM(3:end);
if strcmpi(handles.loadOptionsDICOM,'CreationDate') == 1
    for i = 1:length(handles.filesDICOM)
        dNum(i) = handles.filesDICOM(i).datenum;
    end
    [~,idx] = sort(dNum,'ascend');
    handles.filesDICOM = handles.filesDICOM(idx);
elseif strcmpi(handles.loadOptionsDICOM,'Size') == 1
    for i = 1:length(handles.filesDICOM)
        dNum(i) = handles.filesDICOM(i).fileSize;
    end
    [~,idx] = sort(dNum,'ascend');
    handles.filesDICOM = handles.filesDICOM(idx);
elseif strcmpi(handles.loadOptionsDICOM,'Filename') == 1
    for i = 1:length(handles.filesDICOM)
        dNum{i} = handles.filesDICOM(i).name;
    end
    [~,idx] = sort(dNum);
    handles.filesDICOM = handles.filesDICOM(idx);
end

%update text field with list of Rig files in order

for i = 1:length(handles.filesDICOM)
    dFiles{i} = handles.filesDICOM(i).name;
end
set(handles.textDICOMFileList,'String',dFiles);

guidata(hObject, handles);

% --- Executes on selection change in popupmenuDICOMLoadOptions.
function popupmenuDICOMLoadOptions_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuDICOMLoadOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuDICOMLoadOptions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuDICOMLoadOptions


% --- Executes during object creation, after setting all properties.
function popupmenuDICOMLoadOptions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuDICOMLoadOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonStartAnalysis.
function pushbuttonStartAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonStartAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'pathstrDICOM')
    %generate 2D data from CT scans
    files = handles.filesRig;
    for i = 1:length(handles.filesRig)
        
        
        data = importdata([handles.pathstrRig '\' files(i).name]);%time,position,force
        if strcmp(handles.rig,'Electropulse') == 1
            data = data.data;
            data(:,2) = data(:,7);
            data(:,3) = data(:,8);
            data(:,3) = data(:,3) .* -1;
        end
        data(:,2) = windowAverage(data(:,2),21);
        data(:,3) = data(:,3) .* -1;%correct force direction for later calculations
        moment = data(:,3) .* (handles.beamWidth/4);

        plot(data(:,3),'Parent',handles.axes1);%plot displacement to choose data start
        set(handles.textCurrentObjective,'String',['Select the first point to use for ' files(i).name ' (where the real data begins)']);
%         scrsz = get(0,'ScreenSize');
%         set(gcf,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
        [pointslist,xselect,yselect] = selectdata('selectionmode','closest');
        % close all;

        d0 = data(1,2);
        startFrame = xselect;

        data2 = data(startFrame:end,:);
        moment = moment(startFrame:end);



        plot(data2(:,3),'Parent',handles.axes1);%plot displacement to choose data start
        set(handles.textCurrentObjective,'String',['Select the last point to use for ' files(i).name ' (where the real data ends)']);
%         scrsz = get(0,'ScreenSize');
%         set(gcf,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
        [pointslist,xselect,yselect] = selectdata('selectionmode','closest');
        % close all;

        endFrame = xselect;

        data2 = data2(1:endFrame,:);
        moment = moment(1:endFrame);


        data2(:,2) = data2(:,2) - d0;%subtract offset from displacement
        if mean(data2(:,2)) < 0
            data2(:,2) = data2(:,2) .* -1;
        end

        plot(data2(:,3),'Parent',handles.axes1);
        set(handles.textCurrentObjective,'String','Select the fracture location for this bone');
%         set(gcf,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
        [pointslist,xselect,yselect] = selectdata('selectionmode','closest');
        % close all;

        fractureForce(i) = yselect;
        fractureDisplacement(i) = data2(xselect,2);
        fractureMoment(i) = moment(xselect);
        [ultimateForce(i) I] = max(data2(1:xselect,3));
        ultimateMoment(i) = max(moment(1:xselect));
        ultimateDisplacement(i) = data2(I,2);
        energyAtUltimateForce(i) = sum(data2(1:I,3) .* mean(diff(data2(1:I,2))));
        energyAtFractureForce(i) = sum(data2(1:xselect,3) .* mean(diff(data2(1:xselect,2))));

        plot(data2(:,2),data2(:,3),'Parent',handles.axes1);
        set(handles.textCurrentObjective,'String','Select the points to use in stiffness calculations for this bone (circle them)');
%         set(gcf,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
        [pointslist,xselect,yselect] = selectdata('selectionmode','lasso');
        % close all;

        coeffs = polyfit(xselect,yselect,1);
        stiffness(i) = coeffs(1);

        x1 = linspace(0,max(data2(:,2)));
        y1 = polyval(coeffs,x1);
        aa = length(data2(:,2));
        plot(x1,y1,'red','Parent',handles.axes1);
        hold on;
        plot(data2(1:aa/2,2),data2(1:aa/2,3),'blue','Parent',handles.axes1);
        xlim(handles.axes1,[min(data2(1:aa/2,2)) max(data2(1:aa/2,2))]);
        hold off;
%         set(gcf,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
        set(handles.textCurrentObjective,'String','Select the yield point; ideally, this is where the force v displacement curve deflects, and the stiffness line intersects the curve');
        [pointslist,xselect,yselect] = selectdata('selectionmode','closest');

        yieldForce(i) = yselect{1};
        yieldDisplacement(i) = xselect{1};

        postYieldDisplacement(i) = fractureDisplacement(i) - yieldDisplacement(i);
        yieldMoment(i) = yieldForce(i) * handles.beamWidth / 4;
    
    end
    twoDResults = twoDAnalysisAutomationSub(handles.pathstrDICOM);
    answer = inputdlg('Do you need to use a vector to reorder your 2D data? y/n');
    if strcmpi(answer{1},'y') == 1
        answer = inputdlg('Please enter a vector representing the new order for your 2D data (e.g. 1,3,5,2,4)');
        order = str2num(cell2mat(answer));
        twoDResults = twoDResults(order,:);
    else
    end
    a = length(twoDResults);
    for i = 1:a
        ultimateStress(i) = ultimateMoment(i) * twoDResults(i).data(5) / twoDResults(i).data(3);
        fractureStress(i) = fractureMoment(i) * twoDResults(i).data(5) / twoDResults(i).data(3);
        yieldStress(i) = yieldMoment(i) * twoDResults(i).data(5) / twoDResults(i).data(3);
        rigidity(i) = stiffness(i) * handles.beamWidth^3 / 48;
        youngsModulus(i) = rigidity(i) / twoDResults(i).data(3);
    end
    norms = length(whos('*Stress*'));

    header = {'Rig File Name',...
        'Reorder Vector',...
        'Stiffness (N/mm)',...
        'Yield Load (or Yield Force) (N)',...
        'Maximum Load (or Ultimate Load) (N)',...   
        'Post-Yield Displacement (mm)',...
        'Work to Fracture (N*mm)',...
        'Total Area (Tt.Ar) (mm^2)',...
        'Bone Area (Ct.Ar) (mm^2)',...
        'Medullary Area (Ma.Ar) (mm^2)',...
        'Cortical Thickness (Ct.Th) (mm)',...
        'Polar Moment of Inertia (pMOI)'...
        'Minimum Moment of Inertia (Imin) (mm^4)',...
        'Maximum Moment of Inertia (Imax) (mm^4)',...
        'c or Ymax (mm)',...
        'Young''s Modulus'...
        'Yield Stress (N/mm^2)',...
        'Ultimate Stress (N/mm^2)',...
        };

    if ~isfield(handles,'pathstrResults')
        fid = fopen([handles.pathstrRig '\BendingResults.txt'],'w');
    else
        fid = fopen([handles.pathstrResults '\BendingResults.txt'],'w');
    end

    for i = 1:length(header)
        fprintf(fid,'%s\t',header{i});
    end
    fprintf(fid,'%s\n','');

    for i = 1:length(files)
        fprintf(fid,'%s\t',files(i).name);
        if exist('order') == 1
            fprintf(fid,'%s\t',num2str(order(i)));
        else
            fprintf(fid,'%s\t',num2str(i));
        end
        fprintf(fid,'%s\t',num2str(stiffness(i)));
        fprintf(fid,'%s\t',num2str(yieldForce(i)));
        fprintf(fid,'%s\t',num2str(ultimateForce(i)));
        fprintf(fid,'%s\t',num2str(postYieldDisplacement(i)));
        fprintf(fid,'%s\t',num2str(energyAtFractureForce(i)));
        fprintf(fid,'%s\t',num2str(twoDResults(i).data(6)));
        fprintf(fid,'%s\t',num2str(twoDResults(i).data(1)));
        fprintf(fid,'%s\t',num2str(twoDResults(i).data(7)));
        fprintf(fid,'%s\t',num2str(twoDResults(i).data(14)));
        fprintf(fid,'%s\t',num2str(twoDResults(i).data(2)+twoDResults(i).data(3)));
        fprintf(fid,'%s\t',num2str(twoDResults(i).data(3)));
        fprintf(fid,'%s\t',num2str(twoDResults(i).data(2)));
        fprintf(fid,'%s\t',num2str(twoDResults(i).data(5)));
        fprintf(fid,'%s\t',num2str(youngsModulus(i)));
        fprintf(fid,'%s\t',num2str(yieldStress(i)));
        fprintf(fid,'%s\t',num2str(ultimateStress(i)));

        fprintf(fid,'%s\n','');
    end

    fclose(fid);
    % close all;
    % close all;
else
    %generate 2D data from CT scans
    files = handles.filesRig;
    for i = 1:length(handles.filesRig)
        
        
        data = importdata([handles.pathstrRig '\' files(i).name]);%time,position,force
        if strcmp(handles.rig,'Electropulse') == 1
            data = data.data;
            data(:,2) = data(:,7);
            data(:,3) = data(:,8);
            data(:,3) = data(:,3) .* -1;
        end
        data(:,2) = windowAverage(data(:,2),21);
        data(:,3) = data(:,3) .* -1;%correct force direction for later calculations
        moment = data(:,3) .* (handles.beamWidth/4);

        plot(data(:,3),'Parent',handles.axes1);%plot displacement to choose data start
        set(handles.textCurrentObjective,'String',['Select the first point to use for ' files(i).name ' (where the real data begins)']);
%         scrsz = get(0,'ScreenSize');
%         set(gcf,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
        [pointslist,xselect,yselect] = selectdata('selectionmode','closest');
        % close all;

        d0 = data(1,2);
        startFrame = xselect;

        data2 = data(startFrame:end,:);
        moment = moment(startFrame:end);



        plot(data2(:,3),'Parent',handles.axes1);%plot displacement to choose data start
        set(handles.textCurrentObjective,'String',['Select the last point to use for ' files(i).name ' (where the real data ends)']);
%         scrsz = get(0,'ScreenSize');
%         set(gcf,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
        [pointslist,xselect,yselect] = selectdata('selectionmode','closest');
        % close all;

        endFrame = xselect;

        data2 = data2(1:endFrame,:);
        moment = moment(1:endFrame);


        data2(:,2) = data2(:,2) - d0;%subtract offset from displacement
        if mean(data2(:,2)) < 0
            data2(:,2) = data2(:,2) .* -1;
        end

        plot(data2(:,3),'Parent',handles.axes1);
        set(handles.textCurrentObjective,'String','Select the fracture location for this bone');
%         set(gcf,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
        [pointslist,xselect,yselect] = selectdata('selectionmode','closest');
        % close all;

        fractureForce(i) = yselect;
        fractureDisplacement(i) = data2(xselect,2);
        fractureMoment(i) = moment(xselect);
        [ultimateForce(i) I] = max(data2(1:xselect,3));
        ultimateMoment(i) = max(moment(1:xselect));
        ultimateDisplacement(i) = data2(I,2);
        energyAtUltimateForce(i) = sum(data2(1:I,3) .* mean(diff(data2(1:I,2))));
        energyAtFractureForce(i) = sum(data2(1:xselect,3) .* mean(diff(data2(1:xselect,2))));

        plot(data2(:,2),data2(:,3),'Parent',handles.axes1);
        set(handles.textCurrentObjective,'String','Select the points to use in stiffness calculations for this bone (circle them)');
%         set(gcf,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
        [pointslist,xselect,yselect] = selectdata('selectionmode','lasso');
        % close all;

        coeffs = polyfit(xselect,yselect,1);
        stiffness(i) = coeffs(1);

        aa = length(data2(:,2));
        x1 = linspace(0,max(data2(:,2)));
        y1 = polyval(coeffs,x1);
        plot(data2(1:aa/2,2),data2(1:aa/2,3),'blue');
        hold on;
        plot(x1,y1,'red','Parent',handles.axes1);
        xlim(handles.axes1,[min(data2(1:aa/2,2)) max(data2(1:aa/2,2))]);
        hold off;
%         set(gcf,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
        set(handles.textCurrentObjective,'String','Select the yield point; ideally, this is where the force v displacement curve deflects, and the stiffness line intersects the curve');
        [pointslist,xselect,yselect] = selectdata('selectionmode','closest');

        yieldForce(i) = yselect{2};
        yieldDisplacement(i) = xselect{2};

        postYieldDisplacement(i) = fractureDisplacement(i) - yieldDisplacement(i);
        yieldMoment(i) = yieldForce(i) * handles.beamWidth / 4;
    
    end
%     twoDResults = twoDAnalysisAutomationSub(handles.pathstrDICOM);
%     answer = inputdlg('Do you need to use a vector to reorder your 2D data? y/n');
%     if strcmpi(answer{1},'y') == 1
%         answer = inputdlg('Please enter a vector representing the new order for your 2D data (e.g. 1,3,5,2,4)');
%         order = str2num(cell2mat(answer));
%         twoDResults = twoDResults(order,:);
%     else
%     end
%     a = length(twoDResults);
%     for i = 1:a
%         ultimateStress(i) = ultimateMoment(i) * twoDResults(i).data(5) / twoDResults(i).data(3);
%         fractureStress(i) = fractureMoment(i) * twoDResults(i).data(5) / twoDResults(i).data(3);
%         yieldStress(i) = yieldMoment(i) * twoDResults(i).data(5) / twoDResults(i).data(3);
%         rigidity(i) = stiffness(i) * handles.beamWidth^3 / 48;
%         youngsModulus(i) = rigidity(i) / twoDResults(i).data(3);
%     end
%     norms = length(whos('*Stress*'));

    header = {'LVM File Name',...
        'Sample',...
        'Measurement',...
        'Slices',...
        'Reorder Vector',...
        'Stiffness (N/mm)',...
        'Yield Load (or Yield Force) (N)',...
        'Maximum Load (or Ultimate Load) (N)',...   
        'Post-Yield Displacement (mm)',...
        'Work to Fracture (N*mm)',...
%         'Total Area (Tt.Ar) (mm^2)',...
%         'Bone Area (Ct.Ar) (mm^2)',...
%         'Medullary Area (Ma.Ar) (mm^2)',...
%         'Cortical Thickness (Ct.Th) (mm)',...
%         'Polar Moment of Inertia (pMOI)'...
%         'Minimum Moment of Inertia (Imin) (mm^4)',...
%         'Maximum Moment of Inertia (Imax) (mm^4)',...
%         'c or Ymax (mm)',...
%         'Young''s Modulus'...
%         'Yield Stress (N/mm^2)',...
%         'Ultimate Stress (N/mm^2)',...
        };

    if ~isfield(handles,'pathstrResults')
        fid = fopen([handles.pathstrRig '\BendingResults.txt'],'w');
    else
        fid = fopen([handles.pathstrResults '\BendingResults.txt'],'w');
    end

    for i = 1:length(header)
        fprintf(fid,'%s\t',header{i});
    end
    fprintf(fid,'%s\n','');

    for i = 1:length(files)
        fprintf(fid,'%s\t',files(i).name);
        fprintf(fid,'%s\t','');
        fprintf(fid,'%s\t','');
        fprintf(fid,'%s\t','');
        if exist('order') == 1
            fprintf(fid,'%s\t',num2str(order(i)));
        else
            fprintf(fid,'%s\t',num2str(i));
        end
        fprintf(fid,'%s\t',num2str(stiffness(i)));
        fprintf(fid,'%s\t',num2str(yieldForce(i)));
        fprintf(fid,'%s\t',num2str(ultimateForce(i)));
        fprintf(fid,'%s\t',num2str(postYieldDisplacement(i)));
        fprintf(fid,'%s\t',num2str(energyAtFractureForce(i)));
%         fprintf(fid,'%s\t',num2str(twoDResults(i).data(6)));
%         fprintf(fid,'%s\t',num2str(twoDResults(i).data(1)));
%         fprintf(fid,'%s\t',num2str(twoDResults(i).data(7)));
%         fprintf(fid,'%s\t',num2str(twoDResults(i).data(14)));
%         fprintf(fid,'%s\t',num2str(twoDResults(i).data(2)+twoDResults(i).data(3)));
%         fprintf(fid,'%s\t',num2str(twoDResults(i).data(3)));
%         fprintf(fid,'%s\t',num2str(twoDResults(i).data(2)));
%         fprintf(fid,'%s\t',num2str(twoDResults(i).data(5)));
%         fprintf(fid,'%s\t',num2str(youngsModulus(i)));
%         fprintf(fid,'%s\t',num2str(yieldStress(i)));
%         fprintf(fid,'%s\t',num2str(ultimateStress(i)));

        fprintf(fid,'%s\n','');
    end

    fclose(fid);
    % close all;
    % close all;
end




% --- Executes on button press in pushbuttonSetResultsPath.
function pushbuttonSetResultsPath_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSetResultsPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.pathstrResults = uigetdir(pwd,'Select a folder to output your results');
set(handles.textResultsPath,'String',handles.pathstrResults);

guidata(hObject, handles);


function [twoDResults] = twoDAnalysisAutomationSub(pathstr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%This script automates the 2D analysis Silva lab has always done for 3
%%point bending. You'll need a stack of DICOM files for each midshaft you
%%want to analyze.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all');

answer = inputdlg('Please enter Gaussian sigma (0 for none)');
sigma = str2num(answer{1});
answer = inputdlg('Please enter Gaussian support (0 for none)');
support = str2num(answer{1});

dirs = dir(pathstr);

for i = 3:length(dirs)
    if dirs(i).isdir == 1
        files = dir([pathstr '\' dirs(i).name '\*.dcm*']);
        clear img mask* bw;
        info = dicominfo([pathstr '\' dirs(i).name '\' files(1).name]);
        a = info.Height;
        b = info.Width;
        c = length(files);
      

        img = zeros(a,b,length(files),'int16');
        if sigma == 0 && support == 0
            for l = 1:length(files)
                img(:,:,l) = int16(dicomread([pathstr '\' dirs(i).name '\' files(l).name]));
            end
        elseif sigma ~= 0 && support ~= 0
            for l = 1:length(files)
                img(:,:,l) = int16(imgaussian(double(dicomread([pathstr '\' dirs(i).name '\' files(l).name])),sigma,support));
            end
        else
            error('Sigma and support must make sense! You cannot make only one of those values zero.');
        end
        startSlice = 1;
        endSlice = 50;
        twoDHeader = Analysis_DICOM_Input_Sub_Auto_Thresh(startSlice,endSlice,1,[pathstr '\' dirs(i).name],'y',[pathstr '\' dirs(i).name '2DResults.txt'],1,'n');
    end
end

sysLine = ['del "' pathstr '\*cortical*.txt"'];
system(sysLine);
txtFiles = dir([pathstr '\*.txt']);
for i = 1:length(txtFiles)
    twoDResults(i) = importdata([pathstr '\' txtFiles(i).name]);
    header = ['Sample' twoDHeader];
end

%print out final results
fid = fopen([pathstr '\FinalCorticalResults.txt'],'w');
for i = 1:length(header)
    if i ~= length(header)
        fprintf(fid,'%s\t',header{i});
    else
        fprintf(fid,'%s\n',header{i});
    end
end


for i = 1:length(twoDResults)
    junkText = strfind(txtFiles(i).name,'2DResults.txt');
    sampleText = txtFiles(i).name(1:junkText-1);
    fprintf(fid,'%s\t',sampleText);
    fprintf(fid,'%s\t','');
    for k = 1:length(twoDResults(i).data)-1
        fprintf(fid,'%s\t',num2str(twoDResults(i).data(k)));
    end
    for k = length(twoDResults(i).data)
        fprintf(fid,'%s\n',num2str(twoDResults(i).data(k)));
    end
end

fclose(fid);



function editFileExtension_Callback(hObject, eventdata, handles)
% hObject    handle to editFileExtension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFileExtension as text
%        str2double(get(hObject,'String')) returns contents of editFileExtension as a double


% --- Executes during object creation, after setting all properties.
function editFileExtension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFileExtension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuRig.
function popupmenuRig_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuRig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuRig contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuRig
str = get(handles.popupmenuRig,'String');
val = get(handles.popupmenuRig,'Value');
switch str{val}
    case 'Rig'
        handles.rig = 'Rig';
    case 'Electropulse'
        handles.rig = 'Electropulse';
    case '5866'
        handles.rig = '5866';
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenuRig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuRig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
