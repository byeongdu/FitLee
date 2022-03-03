function varargout = FormFactor_calculator(varargin)
% FORMFACTOR_CALCULATOR MATLAB code for FormFactor_calculator.fig
%      FORMFACTOR_CALCULATOR, by itself, creates a new FORMFACTOR_CALCULATOR or raises the existing
%      singleton*.
%
%      H = FORMFACTOR_CALCULATOR returns the handle to a new FORMFACTOR_CALCULATOR or the handle to
%      the existing singleton*.
%
%      FORMFACTOR_CALCULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FORMFACTOR_CALCULATOR.M with the given input arguments.
%
%      FORMFACTOR_CALCULATOR('Property','Value',...) creates a new FORMFACTOR_CALCULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FormFactor_calculator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FormFactor_calculator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FormFactor_calculator

% Last Modified by GUIDE v2.5 12-Sep-2020 16:22:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FormFactor_calculator_OpeningFcn, ...
                   'gui_OutputFcn',  @FormFactor_calculator_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargin && ischar(varargin{2})
    f = str2func(varargin{2});
    if contains_replace(varargin{2}, 'Callback')
        f([],[], varargin{3});
    else
        f(varargin{3:end});
    end
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FormFactor_calculator is made visible.
function FormFactor_calculator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FormFactor_calculator (see VARARGIN)

% Choose default command line output for FormFactor_calculator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FormFactor_calculator wait for user response (see UIRESUME)
% uiwait(handles.formfactor_calculator);


% --- Outputs from this function are returned to the command line.
function varargout = FormFactor_calculator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in pmshape.
function pmshape_Callback(hObject, eventdata, handles)
% hObject    handle to pmshape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pmshape contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pmshape
%model = get(gcbf, 'userdata');
%N = get(handles.pmnumber, 'value');
shape = get(handles.pmshape, 'string');
Nshape = get(handles.pmshape, 'value');
set(handles.edsize, 'enable', 'on');
set(handles.edsizewidth, 'enable', 'on');
set(handles.edrho, 'enable', 'on');
set(handles.edLshell, 'enable', 'on');
set(handles.edLshell, 'enable', 'on');
switch lower(shape{Nshape})
    case 'sphere'
        set(handles.txtsize, 'string', 'radius');
    case {'pdb', 'boxmodel'}
        pb_loadPDB_Callback(lower(shape{Nshape}))
        set(handles.edsize, 'enable', 'off');
        set(handles.edsizewidth, 'enable', 'off');
        set(handles.edrho, 'enable', 'off');
        set(handles.edLshell, 'enable', 'off');
        set(handles.edLshell, 'enable', 'off');
        %set(handles.edsize, 'string', 1);
        %set(handles.txtsize, 'string', 'edge length');
%        cprintf('*blue', 'Recommended to compute P(q) and Load it\n');
    otherwise 
        set(handles.txtsize, 'string', 'edge length');
end

setvalues(handles)

% --- Executes during object creation, after setting all properties.
function pmshape_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmshape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
polyheds = particle2polyhedron;
% In order to add a polyhedron to FormFactor_calculator, add it to
% particle2polyhedron.m and also F(q) function to AnisoFormFactor.m

shape = {'none', 'sphere', polyheds{:}, 'rho-spherical', 'PDB', 'BOXmodel'};
set(hObject, 'string', shape);

%initialize(handles)



function edposition_Callback(hObject, eventdata, handles)
% hObject    handle to edposition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edposition as text
%        str2double(get(hObject,'String')) returns contents of edposition as a double
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function edposition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edposition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edsize_Callback(hObject, eventdata, handles)
% hObject    handle to edsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edsize as text
%        str2double(get(hObject,'String')) returns contents of edsize as a double
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function edsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'enable', 'off')


function edrho_Callback(hObject, eventdata, handles)
% hObject    handle to edrho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edrho as text
%        str2double(get(hObject,'String')) returns contents of edrho as a double
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function edrho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edrho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'enable', 'off')


function edLshell_Callback(hObject, eventdata, handles)
% hObject    handle to edLshell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edLshell as text
%        str2double(get(hObject,'String')) returns contents of edLshell as a double
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function edLshell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edLshell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'enable', 'off')


function edsizewidth_Callback(hObject, eventdata, handles)
% hObject    handle to edsizewidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edsizewidth as text
%        str2double(get(hObject,'String')) returns contents of edsizewidth as a double
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function edsizewidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edsizewidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'enable', 'off')

% --- Executes on button press in pbdraw.
function pbdraw_Callback(hObject, eventdata, handles)
% hObject    handle to pbdraw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setvalues(handles)
%model = get(gcbf, 'userdata');
model = evalin('base', 'SgAtoms');
%t = numel(model.particle);
N = 1;
if isempty(model{N})
    return;
end
pt{1} = model{N};
pt{1}.position = [0, 0, 0];
%figTag = getappdata(gcbf, 'figTag');
figTag = 'SingleParticle_Fig';
drawparticle(pt, figTag);

% --- Executes on button press in pbadd.
function pbadd_Callback(hObject, eventdata, handles)
% hObject    handle to pbadd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%model = get(gcbf, 'userdata');
model = evalin('base', 'SgAtoms');
%t = numel(model.particle);
t = numel(get(handles.pmnumber, 'string'));
n = numel(model);
model{n+1} = model{n};
%num = numel(md.edensity);
num = numlist2cellstr((1:t+1));
%set(gcbf, 'userdata', model);
assignin('base', 'SgAtoms', model);
set(handles.pmnumber, 'string', num);
%pmnumber_Callback([],[],handles)
%defaultparticle(handles);



function edphi_Callback(hObject, eventdata, handles)
% hObject    handle to edphi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edphi as text
%        str2double(get(hObject,'String')) returns contents of edphi as a double
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function edphi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edphi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtheta_Callback(hObject, eventdata, handles)
% hObject    handle to edtheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtheta as text
%        str2double(get(hObject,'String')) returns contents of edtheta as a double
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function edtheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edpsi_Callback(hObject, eventdata, handles)
% hObject    handle to edpsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edpsi as text
%        str2double(get(hObject,'String')) returns contents of edpsi as a double
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function edpsi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edpsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbdrawinlattice.
function pbdrawinlattice_Callback(hObject, eventdata, handles)
% hObject    handle to pbdrawinlattice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setvalues(handles)
%model = get(gcbf, 'userdata');
model = evalin('base', 'SgAtoms');
%t = numel(model.particle);
currentModelN = get(handles.pmnumber, 'value');
if isempty(model{currentModelN})
    return;
end

pt{1} = model{currentModelN};
figTag = getappdata(gcbf, 'figTag');
sginfo = evalin('base', 'sginfo');
cellinfo = evalin('base', 'cellinfo');
if contains_replace(sginfo.LatticeSystem, 'trigonal')
    if (cellinfo.A==cellinfo.B) && (cellinfo.A==cellinfo.B) && (cellinfo.alpha==cellinfo.beta) && (cellinfo.alpha==cellinfo.gamma)
        sginfo.SymMatrices = sginfo.SymMatricesR;
        sginfo.LatticeSystem = 'trigonal(R)';
        sginfo.NoLatticeCenteringVector = 1;
    end
end
figh = findobj('tag', figTag);
if isempty(figh)
    figh = figure;
    set(figh, 'tag', figTag);
end
azel = get(findobj(figh, 'type', 'axes'), 'view');
housekeepfigure(handles, currentModelN)
setappdata(figh, 'myviewangle', azel);

[p, b, ph] = drawparticle(pt, figTag, sginfo, cellinfo);
set(ph, 'tag', strcat('Model', num2str(currentModelN)));

check_duplicatedparticles(figTag, b, p{1}.color)


function check_duplicatedparticles(figTag, newpo, newparticlecolor)
figh = findobj('tag', figTag);
allparticle = findobj(figh, 'type', 'surface');
newparticle = findobj(figh, 'type', 'surface', 'Facecolor', newparticlecolor);
otherparticle = setdiff(allparticle, newparticle);
if isempty(otherparticle)
    return
end
m = get(otherparticle, 'userdata');

for count1 = 1:numel(newpo)
    for count2 = 1:numel(newpo{count1}.position)/3
        po = newpo{count1}.position(count2, :);
        for i=1:numel(m)
            if iscell(m)
                mn = m{i};
            else
                mn = m;
            end
               
            if isempty(mn)
                break;
            end
            if (sum((mn.positionincell - po).^2) == 0)
                cprintf('Errors', sprintf('Position of particles (%s and %s) overlap at [%i, %i, %i]\n',...
                mn.particle.color, newparticlecolor, po(1), po(2), po(3)));
            end
        end
    end
end


% --- Executes on selection change in pmnumber.
function pmnumber_Callback(hObject, eventdata, handles)
% hObject    handle to pmnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pmnumber contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pmnumber
%model = get(gcbf, 'userdata');
model = evalin('base', 'SgAtoms');
N = get(handles.pmnumber, 'value');
Nparticleall = numel(get(handles.pmnumber, 'string'));
if Nparticleall < numel(model)
    model(Nparticleall+1:end) = [];
    assignin('base', 'SgAtoms', model);
end
C = cellstr(get(handles.pmshape, 'string'));
vshape = findcellstr(C, model{N}.shape);
set(handles.pmshape, 'value', vshape);

% str = [num2str(model{N}.positioninput(1)),',',num2str(model{N}.positioninput(2)),',',...
%     num2str(model{N}.positioninput(3))];
% set(handles.edposition, 'string', ['[', str, ']']);
set(handles.edposition, 'string', model{N}.positioninput);

if isfield(model{N}, 'centering')
    set(handles.ed_centering, 'string', model{N}.centering);
else
    set(handles.ed_centering, 'string', '');
end

if isfield(model{N}, 'position2')
    set(handles.ed_position_LabCoordinate, 'string', ['[',num2str(model{N}.position2), ']']);
else
    set(handles.ed_position_LabCoordinate, 'string', '[]');
end

if isfield(model{N}, 'orientationfactor')
    set(handles.ed_degofrandomorientation, 'string', model{N}.orientationfactor);
else
    set(handles.ed_degofrandomorientation, 'string', 0);
end


set(handles.edsize, 'enable', 'on');
set(handles.edsizewidth, 'enable', 'on');
set(handles.edLshell, 'enable', 'on');
set(handles.ed_PDBfilename, 'enable', 'off');
%set(handles.pb_loadPDB, 'enable', 'on');
% text of txtsize
switch model{N}.shape
    case {'sphere'}
        shapefield = 'radius';
        sigfield = 'radius_sig';
        set(handles.txtsize, 'string', 'radius');
    otherwise
        shapefield = 'edgelength';
        sigfield = 'edgelength_sig';
        set(handles.txtsize, 'string', 'edge length');
        set(handles.ed_degofrandomorientation, 'enable', 'on')
end
switch model{N}.shape
    case {'PDB', 'BOXmodel', 'rho-spherical'}
        set(handles.edsize, 'string', 1)
        set(handles.edsize, 'enable', 'off');
        set(handles.edsizewidth, 'enable', 'off');
        set(handles.edLshell, 'enable', 'off');
        set(handles.ed_PDBfilename, 'enable', 'on');
        cprintf('*red', 'Should load a model using the button "Load".\n')
    otherwise
%        set(handles.ed_PDBfilename, 'enable', 'off');
%        set(handles.pb_loadPDB, 'enable', 'off');
        if isfield(model{N}, 'pdb')
            model{N} = rmfield(model{N}, 'pdb');
        end
        if isfield(model{N}, 'BOXmodel')
            model{N} = rmfield(model{N}, 'BOXmodel');
        end
end

if numel(model{N}.(shapefield)) == 2
    str = [num2str(model{N}.(shapefield)(1)),',',num2str(model{N}.(shapefield)(2))];
    set(handles.edsize, 'string', ['[', str, ']']);
else
    set(handles.edsize, 'string', num2str(model{N}.(shapefield)));
end

set(handles.edsizewidth, 'string', num2str(model{N}.(sigfield)/model{N}.(shapefield)));
if numel(model{N}.rho) == 1
    set(handles.edrho, 'string', num2str(model{N}.rho));
else
    set(handles.edrho, 'string', '[', num2str(model{N}.rho), ']');
end
if isempty(model{N}.Lshell)
    set(handles.edLshell, 'string', '[]');
else
    set(handles.edLshell, 'string', num2str(model{N}.Lshell));
end

check_Pqoption(handles)

if isfield(model{N}, 'pdbfilename')
    set(handles.ed_PDBfilename, 'string', model{N}.pdbfilename);
else
    set(handles.ed_PDBfilename, 'string', '');
end
if isfield(model{N}, 'rhoname')
    set(handles.ed_PDBfilename, 'string', model{N}.rhoname);
else
    set(handles.ed_PDBfilename, 'string', '');
end

% setvalues(handles);
h_radiobutton = findobj(handles.uibg_ParticleOrientation, 'style', 'radiobutton');
orienttype = get(h_radiobutton, 'string');
if isfield(model{N}, 'OrientationType')
    sel = findcellstr(orienttype,model{N}.OrientationType);
else
     sel = findcellstr(orienttype, 'Defined');
end
set(h_radiobutton(sel), 'value', 1);

set(handles.edSOF, 'string', num2str(model{N}.SOF));
set(handles.edcolor, 'string', model{N}.color);
set(handles.edphi, 'string', model{N}.Rotstring1);
set(handles.edtheta, 'string', model{N}.Rotstring2);
set(handles.edpsi, 'string', model{N}.Rotstring3);
% switch lower(model{N}.OrientationType)
%     case 'euler'
%         set(handles.edphi, 'string', num2str(model{N}.rotangle(1)));
%         set(handles.edtheta, 'string', num2str(model{N}.rotangle(2)));
%         set(handles.edpsi, 'string', num2str(model{N}.rotangle(3)));
%     case 'around a vector'
%         set(handles.edphi, 'string', ['[', num2str(model{N}.rotvector), ']']);
%         set(handles.edtheta, 'string', num2str(model{N}.rotangle));
%     case 'xyz->hkl'
%         set(handles.edphi, 'string', ['[', num2str(model{N}.rotvector), ']']);
%         set(handles.edtheta, 'string', num2str(model{N}.rotangle));
% end


% --- Executes during object creation, after setting all properties.
function pmnumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbremove.
function pbremove_Callback(hObject, eventdata, handles)
% hObject    handle to pbremove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
t = numel(get(handles.pmnumber, 'string'));
if t<2
    return;
end
sellayer = get(handles.pmnumber, 'value');
%model = get(gcbf, 'userdata');
model = evalin('base', 'SgAtoms');
%model = gf_model(model, 'remove', 'particle', sellayer);
model(sellayer) = [];
%particle_update(model, sellayer, handles);
%t = numel(model.particle);
%num = numel(md.edensity);
num = numlist2cellstr((1:t-1));
set(handles.pmnumber, 'value', numel(num))
set(handles.pmnumber, 'string', num);

if sellayer == numel(model)
    set(handles.pmnumber, 'value', sellayer-1);
end
%set(gcbf, 'userdata', model);
assignin('base', 'SgAtoms', model);

% display the parameters corresponding to the selected particle number.
pmnumber_Callback([], [], handles)

% --- Executes on button press in pberase.
function pberase_Callback(hObject, eventdata, handles)
% hObject    handle to pberase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function defaultparticle(handles)
maxN = numel(get(handles.pmnumber, 'string'));
set(handles.pmnumber, 'value', maxN);
set(handles.pmshape, 'value', 1);
set(handles.edposition, 'string', '[0, 0, 0]');
set(handles.edsize, 'string', '10');
set(handles.edsizewidth, 'string', '0');
set(handles.edrho, 'string', '1');
set(handles.edLshell, 'string', '[]');
set(handles.edSOF, 'string', '1');
set(handles.edcolor, 'string', 'r');
set(handles.ed_PDBfilename, 'string', '');

%set(handles.ed
function housekeepfigure(handles, currentModelN)
    figTag = 'particlemaker_figure';
    figh = findobj('Tag', figTag);
    delete(findobj(figh, 'type', 'line'));
    
    ax = findobj(figh, 'type', 'axes');
    newparticle = findobj(ax, 'tag', strcat('Model', num2str(currentModelN)));
    delete(findobj(ax, 'type', 'light'));
    allparticle = get(ax, 'children');
    if ~isempty(allparticle)
        allparticle(~isprop(allparticle, 'facealpha')) = [];
    end
    
    if numel(get(handles.pmnumber, 'string')) == 1
        if ~isempty(figh)
            figure(figh)
            clf
        else
            figh = figure;
            set(figh, 'tag', figTag);
        end
    else
        %allparticle = findobj(figh, 'type', 'surface');
        %colr = get(allparticle{1}, 'Facecolor');
        %newparticle = findobj(figh, 'type', 'surface', 'Facecolor', colr);
        if ~isempty(allparticle)
            otherparticle = setdiff(allparticle, newparticle);
            set(otherparticle, 'Facealpha', 0.2);
            delete(newparticle);
        end
    end

function setvalues(handles)
%model = get(gcbf, 'userdata');
% try
%     model = evalin('base', 'SgAtoms');
% catch
%     model = [];
% end
model = [];
N = 1;
shape = get(handles.pmshape, 'string');
Nshape = get(handles.pmshape, 'value');
model{N}.shape = shape{Nshape};
%model{N}.positioninput = eval(get(handles.edposition, 'string'));

model{N}.orientationfactor = str2double(get(handles.ed_degofrandomorientation, 'string'));

switch model{N}.shape
    case 'sphere'
        shapefield = 'radius';
        sigfield = 'radius_sig';
    case 'PDB'
        shapefield = 'edgelength';
        sigfield = 'edgelength_sig';
        set(handles.edsize, 'string', num2str(1));
    otherwise
        shapefield = 'edgelength';
        sigfield = 'edgelength_sig';
end

model{N}.(shapefield) = eval(get(handles.edsize, 'string'));
sig = eval(get(handles.edsizewidth, 'string'));
if sig > 1
    error('sig should be less than 1')
end

if numel(model{N}.(shapefield)) == numel(sig)
    model{N}.(sigfield) = model{N}.(shapefield).*sig;
elseif numel(sig) == 1
    model{N}.(sigfield) = model{N}.(shapefield)*sig;
elseif numel(model{N}.(shapefield)) == 1
    error('The number of sig is larger than that of edgelength');
end
model{N}.rho = eval(get(handles.edrho, 'string'));
model{N}.Lshell = eval(get(handles.edLshell, 'string'));
if (model{N}.(shapefield) == 0)
    cprintf('red', '%s of the model %d is 0\n', shapefield, N);
    model{N} = [];
    %set(gcbf, 'userdata', model);
    assignin('base', 'SgAtoms', model);
    return
end
if ~isempty(model{N}.Lshell)
    if model{N}.Lshell == 0
        model{N}.Lshell = [];
    end
end

model{N}.color = get(handles.edcolor, 'string');
model{N}.TfactorB = 0;

h = get(handles.uibg_ParticleOrientation, 'SelectedObject');
model{N}.OrientationType = get(h, 'string');
model{N}.Rotmat = eye(3);
model{N}.Rotstring1 = get(handles.edphi, 'string');
model{N}.Rotstring2 = get(handles.edtheta, 'string');
model{N}.Rotstring3 = get(handles.edpsi, 'string');
switch lower(model{N}.OrientationType)
    case 'euler'
        phi = str2double(get(handles.edphi, 'string'));
        theta = str2double(get(handles.edtheta, 'string'));
        psi = str2double(get(handles.edpsi, 'string'));

        model{N}.rotangle = [phi, theta, psi];
        if (phi~=0)|(theta~=0)|(psi~=0)
            Rmat = eulerAnglesToRotation3d(phi, theta, psi);
            model{N}.Rotmat = Rmat(1:3, 1:3);
        end
    case 'around a vector'
        vector = eval(get(handles.edphi, 'string'));
        angle = str2double(get(handles.edtheta, 'string'));
        model{N}.Rotmat = rotate_around_vector(vector, angle);
        model{N}.rotvector = vector;
        model{N}.rotangle = angle;
    case 'xyz->hkl'
        try
            hklx = eval(get(handles.edphi, 'string'));
        catch
            hklx = [];
        end
        try
            hkly = eval(get(handles.edtheta, 'string'));
        catch
            hkly = [];
        end
        try
            hklz = eval(get(handles.edpsi, 'string'));
        catch
            hklz = [];
        end

        XYZ = [1, 0, 0; 0, 1, 0; 0, 0, 1];
        if isempty(hklx) || numel(hklx) ~= 3
            hklx = [0, 0, 0];
        end
        if isempty(hkly) || numel(hkly) ~= 3
            hkly = [0, 0, 0];
        end
        if isempty(hklz) || numel(hklz) ~= 3
            hklz = [0, 0, 0];
        end
        isempty_rvector = [(isempty(hklx)|| sum(hklx.^2) == 0),...
            (isempty(hkly)|| sum(hkly.^2) == 0),...
            isempty(hklz)|| sum(hklz.^2) == 0];
        
        hkl = [hklx; hkly; hklz]; % a*; b*; c*
        hkl = hkl(~isempty_rvector, :);
        XYZ = XYZ(~isempty_rvector, :);
        cellinfo = evalin('base', 'cellinfo');
        HKL = cellinfo.recimat*hkl';
        HKL = HKL';
        R = eye(3);
        for i=1:size(HKL, 1)
            Ri = rotate_between_vectors(HKL(i, :), XYZ(i, :));
            R = R*Ri;
        end
        model{N}.Rotmat = R;
end
%model{N}.Rotmat = R2*model{N}.Rotmat;

%set(gcbf, 'userdata', model);
assignin('base', 'SgAtoms', model);
% Since the following function also load and save SgAtoms, it should be
% after assignin... line.
uibg_ParticleOrientation_SelectionChangedFcn(h, [], handles);

% --- Executes during object creation, after setting all properties.
function particlemaker_pdb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to formfactor_calculator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
figTag = 'particlemaker_figure';
setappdata(gcbf, 'figTag', figTag);



function edcolor_Callback(hObject, eventdata, handles)
% hObject    handle to edcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edcolor as text
%        str2double(get(hObject,'String')) returns contents of edcolor as a double
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function edcolor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'enable', 'off')

% --- Executes on button press in pbloadaset.
function pbloadaset_Callback(hObject, eventdata, handles)
% hObject    handle to pbloadaset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [filename, pathname] = uigetfile( ...
       {'*.mat','mat file (*.dat)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a file', ...
        'MultiSelect', 'off');
 
    % This code checks if the user pressed cancel on the dialog.
 
    if isequal(filename,0) || isequal(pathname,0)
       disp('User pressed cancel')
       return
    else
       disp(['User selected ', fullfile(pathname, filename)])
    end
model = load(fullfile(pathname, filename));
%set(gcbf, 'userdata', model);
assignin('base', 'SgAtoms', model);
%t = numel(model.particle);
n = numel(model);
%num = numel(md.edensity);
num = numlist2cellstr((1:n));
assignin('base', 'SgAtoms', model);
set(handles.pmnumber, 'string', num);
%defaultparticle(handles);


% --- Executes on button press in pbsave.
function pbsave_Callback(hObject, eventdata, handles)
% hObject    handle to pbsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% [filename, pathname] = uiputfile( ...
%    '*.mat', ...
%     'Save as', ...
%     'MultiSelect', 'off');
% 
% % This code checks if the user pressed cancel on the dialog.
% 
% if isequal(filename,0) || isequal(pathname,0)
%    disp('User pressed cancel')
%    return
% else
%    disp(['User selected ', fullfile(pathname, filename)])
% end

%fn = fullfile(pathname, filename);
%model = get(gcbf, 'userdata');
model = evalin('base', 'SgAtoms');
%N = get(handles.pmnumber, 'value');
N = 1;
%save(fn, model{N});
model = model{N};
uisave('model', [model.shape, '.mat']);


% --- Executes on button press in pbload.
function pbload_Callback(hObject, eventdata, handles)
% hObject    handle to pbload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [filename, pathname] = uigetfile( ...
       {'*.mat','mat file (*.dat)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a file', ...
        'MultiSelect', 'off');
 
    % This code checks if the user pressed cancel on the dialog.
 
    if isequal(filename,0) || isequal(pathname,0)
       disp('User pressed cancel')
       return
    else
       disp(['User selected ', fullfile(pathname, filename)])
    end
    
    p = load(fullfile(pathname, filename));
%model = get(gcbf, 'userdata');
try
    model = evalin('base', 'SgAtoms');
catch
    model = [];
end

N = get(handles.pmnumber, 'value');
model{N} = p.model;
%set(gcbf, 'userdata', model);
assignin('base', 'SgAtoms', model);

pmnumber_Callback(hObject, eventdata, handles)



function ed_Pqmononame_Callback(hObject, eventdata, handles)
% hObject    handle to ed_Pqmononame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_Pqmononame as text
%        str2double(get(hObject,'String')) returns contents of ed_Pqmononame as a double
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function ed_Pqmononame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_Pqmononame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_loadPq.
function pb_loadPq_Callback(hObject, eventdata, handles)
% hObject    handle to pb_loadPq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model = evalin('base', 'SgAtoms');
N = get(handles.pmnumber, 'value');

    [filename, pathname] = uigetfile( ...
       {'*.mat','mat file (*.dat)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a file', ...
        'MultiSelect', 'off');
 
    % This code checks if the user pressed cancel on the dialog.
 
    if isequal(filename,0) || isequal(pathname,0)
       disp('User pressed cancel')
       return
    else
       disp(['User selected ', fullfile(pathname, filename)])
    end
    
p = load(fullfile(pathname, filename));

model{N}.Pq = p;

assignin('base', 'SgAtoms', model);


% --- Executes on button press in pb_PickPqname.
function pb_PickPqname_Callback(hObject, eventdata, handles)
% hObject    handle to pb_PickPqname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [filename, pathname] = uigetfile( ...
       {'*.dat','dat file (*.dat)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a file', ...
        'MultiSelect', 'off');
 
    % This code checks if the user pressed cancel on the dialog.
 
    if isequal(filename,0) || isequal(pathname,0)
       disp('User pressed cancel')
       return
    else
       disp(['User selected ', fullfile(pathname, filename)])
    end
    
pb_clearPqfields_Callback(hObject, eventdata, handles)
model = evalin('base', 'SgAtoms');
N = get(handles.pmnumber, 'value');
[a, b] = hdrload(fullfile(pathname, filename));
[keyw, remw] = strtok(a);
Pqmono.size = str2double(remw(2:end));
Pqmono.Pq = b(:, 1:2);
Pqmono.Fq_o = [b(:,1), b(:,3)+sqrt(-1)*b(:,4)];
model{N}.Pqmono = Pqmono;
assignin('base', 'SgAtoms', model);

% fnstr = get(handles.ed_Pqmononame, 'string');
% 
% if isempty(fnstr)
%     try
%     model{N} = rmfield(model{N}, 'Pqname');
%     end
% else
%    model{N}.Pqname = fnstr;
% end
    
% set(handles.ed_Pqmononame, 'string', fullfile(pathname, filename));
% setvalues(handles)


% --- Executes on button press in pb_loadPDB.
function pb_loadPDB_Callback(vshape)
% hObject    handle to pb_loadPDB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%N = get(handles.pmnumber, 'value');
N = 1;
%Nsh = get(handles.pmshape, 'value');
%C = cellstr(get(handles.pmshape, 'string'));
%vshape = C{Nsh};
switch lower(vshape)
    case 'pdb'
        fileselectstr = {'*.pdb','PDB file (*.pdb)'; ...
   '*.mat','MAT file (*.mat)'; ...
    '*.*',  'All Files (*.*)'};
    case 'boxmodel'
        fileselectstr = {'*.mat','MAT file (*.mat)'; ...
    '*.*',  'All Files (*.*)'};
    case 'rho-spherical'
        fileselectstr = {'*.dat','DAT file (*.dat)'; ...
    '*.*',  'All Files (*.*)'};
    otherwise
        cprintf('*red', 'Use this only for PDB, BOXmodel and rho-spherical.\n');
        return
end
%set(handles.pmshape, 'value', vshape);

[filename, pathname] = uigetfile( ...
    fileselectstr, ...
    'Pick a file', ...
    'MultiSelect', 'off');

% This code checks if the user pressed cancel on the dialog.

if isequal(filename,0) || isequal(pathname,0)
   disp('User pressed cancel')
   return
else
   disp(['User selected ', fullfile(pathname, filename)])
end
FN = fullfile(pathname, filename);
[~, PDBname, ext] = fileparts(FN);

model = evalin('base', 'SgAtoms');

switch vshape
    case 'rho-spherical'
        dt = load(FN);
        model{N}.rhomodel = dt;
        model{N}.rhoname = PDBname;
    case 'PDB'
        switch ext
            case '.pdb'
                cprintf('*blue', 'Wait for loading.....\n');
                pdb = PDB(FN);
                cprintf('*blue', 'File "%s" are loaded.\n', FN);
            case '.mat'
                dt = load(FN);
                pdb = dt.pdb;
        end
        model{N}.pdb = pdb;
        model{N}.pdbfilename = PDBname;
    case 'BOXmodel'
        dt = load(FN);
        dtfieldName = fieldnames(dt);
        model{N}.BOXmodel = dt.(dtfieldName{1});
        model{N}.BOXname = PDBname;
end

set(handles.ed_PDBfilename, 'string', PDBname);
assignin('base', 'SgAtoms', model);



function ed_PDBfilename_Callback(hObject, eventdata, handles)
% hObject    handle to ed_PDBfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_PDBfilename as text
%        str2double(get(hObject,'String')) returns contents of ed_PDBfilename as a double


% --- Executes during object creation, after setting all properties.
function ed_PDBfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_PDBfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_calPq.
function pb_calPq_Callback(hObject, eventdata, handles)
% hObject    handle to pb_calPq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model = evalin('base', 'SgAtoms');
N = get(handles.pmnumber, 'value');
%[Sq, r,pr] = debyepr(q, x, y, z, 500, atm, Fq)
q = evalin('base', 'q');
cprintf('*blue', 'Computation started...Wait.\n');
%Pq = pq_pdb(q, pdb.atom);
p{1} = model{N};

[Pqall, ~, Fqo, ~, Fq_os] = Pq_general(p, q);
model{N}.Pq = [q(:), Pqall];
model{N}.Fq_o = [q(:), Fqo];
model{N}.Fq_os = [q(:), Fq_os];

%[~, Pq] = AnisoFormFactor(model{N}, q, [], 3);
cprintf('*blue', 'Computation done.\n');
%model{N}.Pq = [q(:), Pq(:)];
figure;
loglog(q, Pqall);
k = [q(:), Pqall];
assignin('base', 'SgAtoms', model);
% save
% This code checks if the user pressed cancel on the dialog.
[filename, pathname] = uiputfile( ...
   {'*.dat','Dat file (*.dat)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Save data');

% This code checks if the user pressed cancel on the dialog.

if isequal(filename,0) || isequal(pathname,0)
   disp('User pressed cancel')
   return
end
FN = fullfile(pathname, filename);
save(FN, 'k', '-ascii');


% --- Executes when selected object is changed in uibg_ParticleOrientation.
function uibg_ParticleOrientation_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibg_ParticleOrientation 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    model = evalin('base', 'SgAtoms');
    N = get(handles.pmnumber, 'value');
    model{N}.OrientationType = get(hObject, 'string');
catch
    model = [];
end

switch lower(get(hObject, 'string'))
    case 'random'
        if ~isempty(model)
            if ~isfield(model{N}, 'Fq')
                if isfield(model{N}, 'Fqmean')
                    Fqmean = model{N}.Fqmean;
                else
                    q = evalin('base', 'q');
                    cprintf('*blue', 'Wait for computing <F(q)>.\n');
                    [~, ~, Fqmean] = AnisoFormFactor(model{N}, q, [], 3);
                    cprintf('*blue', 'Done.\n');
                    Fqmean = [sqrt(q(:,1).^2+q(:,2).^2+q(:,3).^2), Fqmean(:)];
                end
                model{N}.Fq = Fqmean;
            end
        end
    case 'around a vector'
        set(handles.txtphi, 'string', 'vector')
        set(handles.txttheta, 'string', 'Angle(deg)')
        set(handles.txtpsi, 'enable', 'off')
    case 'xyz->hkl'
        set(handles.txtphi, 'string', 'X')
        set(handles.txttheta, 'string', 'Y')
        set(handles.txtpsi, 'string', 'Z')
        set(handles.txtpsi, 'enable', 'on')
    case 'euler'
        set(handles.txtphi, 'string', 'Rotangle: phi')
        set(handles.txttheta, 'string', 'Rotangle2: theta')
        set(handles.txtpsi, 'enable', 'on')
        set(handles.txtpsi, 'string', 'Rotangle3: psi')
    otherwise
        if ~isempty(model)
            if isfield(model{N}, 'Fq')
                model{N}.Fqmean = model{N}.Fq;
                model{N} = rmfield(model{N}, 'Fq');
            end
        end
end
if ~isempty(model)
    assignin('base', 'SgAtoms', model);
end


% --- Executes on button press in pb_Fqcal.
function pb_Fqcal_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Fqcal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model = evalin('base', 'SgAtoms');
N = get(handles.pmnumber, 'value');
q = evalin('base', 'q');
cprintf('*blue', 'Wait for computing <F(q)>.\n');
[~, ~, Fqmean] = AnisoFormFactor(model{N}, q, [], 3);
cprintf('*blue', 'Done.\n');
Fqmean = [q(:), Fqmean(:)];
model{N}.Fqmean = Fqmean;
model{N}.Fq = Fqmean;
assignin('base', 'SgAtoms', model);



function ed_degofrandomorientation_Callback(hObject, eventdata, handles)
% hObject    handle to ed_degofrandomorientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_degofrandomorientation as text
%        str2double(get(hObject,'String')) returns contents of ed_degofrandomorientation as a double
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function ed_degofrandomorientation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_degofrandomorientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'enable', 'off')

% --- Executes on button press in pb_clearPqfields.
function pb_clearPqfields_Callback(hObject, eventdata, handles)
% hObject    handle to pb_clearPqfields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model = evalin('base', 'SgAtoms');
N = get(handles.pmnumber, 'value');
%[Sq, r,pr] = debyepr(q, x, y, z, 500, atm, Fq)
%Pq = pq_pdb(q, pdb.atom);
if isfield(model{N}, 'Pq')
    model{N} = rmfield(model{N}, 'Pq');
    model{N} = rmfield(model{N}, 'Fq_o');
    model{N} = rmfield(model{N}, 'Fq_os');
end
if isfield(model{N}, 'Pqpdb')
    model{N} = rmfield(model{N}, 'Pqpdb');
end
if isfield(model{N}, 'Pqmono')
    model{N} = rmfield(model{N}, 'Pqmono');
end
assignin('base', 'SgAtoms', model);


% --- Executes on button press in rb_Pqoption3.
function rb_Pqoption3_Callback(hObject, eventdata, handles)
% hObject    handle to rb_Pqoption3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_Pqoption3


% --- Executes when selected object is changed in ui_SimulationMode.
function uibg_isPq_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in ui_SimulationMode 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in rb_Pqoption2.
function rb_Pqoption2_Callback(hObject, eventdata, handles)
% hObject    handle to rb_Pqoption2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_Pqoption2

function check_Pqoption(varargin)
handles = varargin{1};
try
    model = evalin('base', 'SgAtoms');
catch
    return;
end
N = get(handles.pmnumber, 'value');
option = 1;
if isfield(model{N}, 'Pqmono')
    option = 2;
    if isfield(model{N}, 'Pq')
        model{N} = rmfield(model{N}, 'Pq');
    end
    if isfield(model{N}, 'Fq_o')
        model{N} = rmfield(model{N}, 'Fq_o');
    end
    if isfield(model{N}, 'Fq_os')
        model{N} = rmfield(model{N}, 'Fq_os');
    end
elseif isfield(model{N}, 'pdb')
    option = 4;
else
    if isfield(model{N}, 'Pq')
        option = 3;
    end
    if isfield(model{N}, 'Fq_o')
        option = 3;
    end
    if isfield(model{N}, 'Fq_o')
        option = 3;
    end
end
switch option
    case 1
        set(handles.rb_Pqoption1, 'value', 1);
    case 2
        set(handles.rb_Pqoption2, 'value', 1);
    case 3
        set(handles.rb_Pqoption3, 'value', 1);
    case 4
        set(handles.rb_Pqoption4, 'value', 1);
end
assignin('base', 'SgAtoms', model);


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model = evalin('base', 'SgAtoms');
N = get(handles.pmnumber, 'value');
%[Sq, r,pr] = debyepr(q, x, y, z, 500, atm, Fq)
q = evalin('base', 'q');
cprintf('*blue', 'Computation started...Wait.\n');
%Pq = pq_pdb(q, pdb.atom);

p{1} = model{N};
switch model{N}.shape
    case 'sphere'
        shapefield = 'radius';
        sigfield = 'radius_sig';
    case 'PDB'
        shapefield = 'edgelength';
        sigfield = 'edgelength_sig';
        set(handles.edsize, 'string', num2str(1));
    otherwise
        shapefield = 'edgelength';
        sigfield = 'edgelength_sig';
end

p{1}.(sigfield) = 0;

[Pqmono, ~, Fqo] = Pq_general(p, q);
Pq = [];
Pq.Pq = [q(:), Pqmono];
Pq.Fq_o = [q(:), Fqo];
Pq.size = p{1}.(shapefield)(1);
if isfield(model{N}, 'Pq')
    model{N} = rmfield(model{N}, 'Pq');
end
model{N}.Pqmono = Pq;
%model{N}.Fq_os = [q(:), Fq_os];

%[~, Pq] = AnisoFormFactor(model{N}, q, [], 3);
cprintf('*blue', 'Computation done.\n');
assignin('base', 'SgAtoms', model);


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile( ...
   {'*.dat','Dat file (*.dat)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Save data');

% This code checks if the user pressed cancel on the dialog.

if isequal(filename,0) || isequal(pathname,0)
   disp('User pressed cancel')
   return
end

model = evalin('base', 'SgAtoms');
q = evalin('base', 'q');

N = get(handles.pmnumber, 'value');
if isfield(model{N}, 'Pqmono')
    Pqmono = model{N}.Pqmono;
    k = [q(:), Pqmono.Pq(:, 2), ...
        real(Pqmono.Fq_o(:,2)), imag(Pqmono.Fq_o(:,2))];
else
    errordlg('Dude! Calculate before trying to save.', 'My Error Dialog');
    return
end

% switch model{N}.shape
%     case 'sphere'
%         shapefield = 'radius';
%         sigfield = 'radius_sig';
%     case 'PDB'
%         shapefield = 'edgelength';
%         sigfield = 'edgelength_sig';
%         set(handles.edsize, 'string', num2str(1));
%     otherwise
%         shapefield = 'edgelength';
%         sigfield = 'edgelength_sig';
% end

FN = fullfile(pathname, filename);
fid = fopen(FN, 'w');
%fprintf(fid, 'Size : %0.5e\n', model{N}.(shapefield));
fprintf(fid, 'Size : %0.5e\n', model{N}.Pqmono.size);
fprintf(fid, '%0.8e %0.8e %0.8e %0.8e\n', k');
fclose(fid);
%save(FN, 'k', '-ascii');

% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model = evalin('base', 'SgAtoms');
N = get(handles.pmnumber, 'value');
%[Sq, r,pr] = debyepr(q, x, y, z, 500, atm, Fq)
    [filename, pathname] = uigetfile( ...
       {'*.dat','dat file (*.dat)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a file', ...
        'MultiSelect', 'off');
 
    % This code checks if the user pressed cancel on the dialog.
 
    if isequal(filename,0) || isequal(pathname,0)
       disp('User pressed cancel')
       return
    else
       disp(['User selected ', fullfile(pathname, filename)])
    end
    
k = load(fullfile(pathname, filename));

model{N}.Pq = [k(:,1), k(:,2)];
model{N}.Fq_o = [k(:,1), k(:,3)+sqrt(-1)*k(:,4)];
model{N}.Fq_os = [k(:,1), k(:,5)+sqrt(-1)*k(:,6)];

%[~, Pq] = AnisoFormFactor(model{N}, q, [], 3);
assignin('base', 'SgAtoms', model);

% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model = evalin('base', 'SgAtoms');
N = get(handles.pmnumber, 'value');
%[Sq, r,pr] = debyepr(q, x, y, z, 500, atm, Fq)
q = evalin('base', 'q');
cprintf('*blue', 'Computation started...Wait.\n');
%Pq = pq_pdb(q, pdb.atom);

p{1} = model{N};
switch model{N}.shape
    case 'sphere'
        shapefield = 'radius';
        sigfield = 'radius_sig';
    case 'PDB'
        shapefield = 'edgelength';
        sigfield = 'edgelength_sig';
        set(handles.edsize, 'string', num2str(1));
    otherwise
        shapefield = 'edgelength';
        sigfield = 'edgelength_sig';
end

%p{1}.(shapefield) = 0;

[Pqall, ~, Fqo, ~, Fq_os] = Pq_general(p, q);
model{N}.Pq = [q(:), Pqall];
model{N}.Fq_o = [q(:), Fqo];
model{N}.Fq_os = [q(:), Fq_os];

%[~, Pq] = AnisoFormFactor(model{N}, q, [], 3);
cprintf('*blue', 'Computation done.\n');
assignin('base', 'SgAtoms', model);


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile( ...
   {'*.dat','Dat file (*.dat)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Save data');

% This code checks if the user pressed cancel on the dialog.

if isequal(filename,0) || isequal(pathname,0)
   disp('User pressed cancel')
   return
end

model = evalin('base', 'SgAtoms');
N = get(handles.pmnumber, 'value');
if isfield(model{N}, 'Pq')
    k = [model{N}.Pq(:, 1), model{N}.Pq(:, 2), ...
        real(model{N}.Fq_o(:, 2)), imag(model{N}.Fq_o(:, 2)), ...
        real(model{N}.Fq_os(:, 2)), imag(model{N}.Fq_os(:, 2))];
else
    errordlg('Dude! Calculate before trying to save.', 'My Error Dialog');
    return
end

FN = fullfile(pathname, filename);
save(FN, 'k', '-ascii');


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model = evalin('base', 'SgAtoms');
N = get(handles.pmnumber, 'value');
%[Sq, r,pr] = debyepr(q, x, y, z, 500, atm, Fq)
    [filename, pathname] = uigetfile( ...
       {'*.dat','dat file (*.dat)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a file', ...
        'MultiSelect', 'off');
 
    % This code checks if the user pressed cancel on the dialog.
 
    if isequal(filename,0) || isequal(pathname,0)
       disp('User pressed cancel')
       return
    else
       disp(['User selected ', fullfile(pathname, filename)])
    end
pb_clearPqfields_Callback(hObject, eventdata, handles)
model = evalin('base', 'SgAtoms');
N = get(handles.pmnumber, 'value');
[a, b] = hdrload(fullfile(pathname, filename));
[keyw, remw] = strtok(a);
Pqmono.size = str2double(remw(2:end));
Pqmono.Pq = b(:, 1:2);
Pqmono.Fq_o = [b(:,1), b(:,3)+sqrt(-1)*b(:,4)];
model{N}.Pqpdb = Pqmono;
assignin('base', 'SgAtoms', model);
    
% k = load(fullfile(pathname, filename));
% 
% model{N}.Pq = [k(:,1), k(:,2)];
% model{N}.Fq_o = [k(:,1), k(:,3)+sqrt(-1)*k(:,4)];
% 
% %[~, Pq] = AnisoFormFactor(model{N}, q, [], 3);
% assignin('base', 'SgAtoms', model);

% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
model = evalin('base', 'SgAtoms');
N = get(handles.pmnumber, 'value');
%[Sq, r,pr] = debyepr(q, x, y, z, 500, atm, Fq)
q = evalin('base', 'q');
cprintf('*blue', 'Computation started...Wait.\n');
%Pq = pq_pdb(q, pdb.atom);
% model{N}.pdb

p{1} = model{N};
if ~isfield(p{1}, 'pdb')
    errordlg('Dude! Load PDB file.', 'My Error Dialog');
    return
end

%Pqall = Pq_general(p, q);
[Pq, Fqm] = pq_pdb(q, p{1}.pdb, 3);
Pqpdb.Pq = [q(:), Pq];
Pqpdb.Fq_o = [q(:), Fqm];
model{N}.Pqpdb = Pqpdb;
if isfield(model{N}, 'Pq')
    model{N} = rmfield(model{N}, 'Pq');
end

%model{N}.Fq_os = [q(:), Fq_os];

%[~, Pq] = AnisoFormFactor(model{N}, q, [], 3);
cprintf('*blue', 'Computation done.\n');
assignin('base', 'SgAtoms', model);

% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile( ...
   {'*.dat','Dat file (*.dat)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Save data');

% This code checks if the user pressed cancel on the dialog.

if isequal(filename,0) || isequal(pathname,0)
   disp('User pressed cancel')
   return
end

model = evalin('base', 'SgAtoms');
N = get(handles.pmnumber, 'value');
if isfield(model{N}, 'Pqpdb')
    Pqpdb = model{N}.Pqpdb;
    k = [Pqpdb.Pq(:, 1), Pqpdb.Pq(:, 2), ...
        real(Pqpdb.Fq_o(:, 2)), imag(Pqpdb.Fq_o(:, 2))];
else
    errordlg('Dude! Calculate before trying to save.', 'My Error Dialog');
    return
end

FN = fullfile(pathname, filename);
save(FN, 'k', '-ascii');


% --- Executes on button press in rb_orientationtype1.
function rb_orientationtype1_Callback(hObject, eventdata, handles)
% hObject    handle to rb_orientationtype1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_orientationtype1



function ed_centering_Callback(hObject, eventdata, handles)
% hObject    handle to ed_centering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_centering as text
%        str2double(get(hObject,'String')) returns contents of ed_centering as a double
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function ed_centering_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_centering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_position_LabCoordinate_Callback(hObject, eventdata, handles)
% hObject    handle to ed_position_LabCoordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_position_LabCoordinate as text
%        str2double(get(hObject,'String')) returns contents of ed_position_LabCoordinate as a double
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function ed_position_LabCoordinate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_position_LabCoordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rb_orientationtype4.
function rb_orientationtype4_Callback(hObject, eventdata, handles)
% hObject    handle to rb_orientationtype4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_orientationtype4


% --- Executes on selection change in pm_wyckoff.
function pm_wyckoff_Callback(hObject, eventdata, handles)
% hObject    handle to pm_wyckoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_wyckoff contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_wyckoff
% sgh = findobj('tag', 'spacegroup');
% ud = get(sgh, 'userdata');
% sgv = evalin('base', 'sginfo');
% wy = ud.wy;
% wyck = wy{sgv.Number};
wy = get(hObject, 'userdata');
set(handles.edposition, 'string', strcat('[',wy{get(hObject,'Value'), 4},']'));
try
    setvalues(handles)
catch
    
end

% --- Executes during object creation, after setting all properties.
function pm_wyckoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_wyckoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pm_rho.
function pm_rho_Callback(hObject, eventdata, handles)
% hObject    handle to pm_rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_rho contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_rho
ud = get(hObject, 'userdata');
set(handles.edrho, 'string', ud{get(hObject,'Value'),2})
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function pm_rho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
str = [{'Arbitrary', 1};
    {'Au', edensity('Au')};
    {'Ag',   edensity('Ag')};
    {'Pt', edensity('Pt')};
    {'Empty', 0}];
set(hObject, 'string', str(:,1));
set(hObject, 'userdata', str)

% --- Executes on selection change in pm_color.
function pm_color_Callback(hObject, eventdata, handles)
% hObject    handle to pm_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_color contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_color
ud = get(hObject, 'userdata');
set(handles.edcolor, 'string', ud{get(hObject,'Value'),2})
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function pm_color_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
str = [{'red', 'r'};
    {'green', 'g'};
    {'blue', 'b'};
    {'yellow', 'y'};
    {'magenta', 'm'};
    {'black', 'k'}];
set(hObject, 'string', str(:,1));
set(hObject, 'userdata', str)

% --- Executes on selection change in pm_orientationfactor.
function pm_orientationfactor_Callback(hObject, eventdata, handles)
% hObject    handle to pm_orientationfactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_orientationfactor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_orientationfactor
ud = get(hObject, 'userdata');
set(handles.ed_degofrandomorientation, 'string', ud{get(hObject,'Value'),2})
setvalues(handles)

% --- Executes during object creation, after setting all properties.
function pm_orientationfactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_orientationfactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
str = [{'Ordered', 0};
    {'Random', 1};
    {'otherwise', 'x'}];
set(hObject, 'string', str(:,1));
set(hObject, 'userdata', str);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnAddParticle_Callback(hObject, eventdata, handles)
% hObject    handle to mnAddParticle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pbadd_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function mnRemoveParticle_Callback(hObject, eventdata, handles)
% hObject    handle to mnRemoveParticle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
t = numel(get(handles.pmnumber, 'string'));
if t<2
    return;
end
sellayer = get(handles.pmnumber, 'value');
%model = get(gcbf, 'userdata');
model = evalin('base', 'SgAtoms');
%model = gf_model(model, 'remove', 'particle', sellayer);
model(sellayer) = [];
%particle_update(model, sellayer, handles);
%t = numel(model.particle);
%num = numel(md.edensity);
num = numlist2cellstr((1:t-1));
set(handles.pmnumber, 'value', numel(num))
set(handles.pmnumber, 'string', num);

if sellayer == numel(model)
    set(handles.pmnumber, 'value', sellayer-1);
end
%set(gcbf, 'userdata', model);
assignin('base', 'SgAtoms', model);

% display the parameters corresponding to the selected particle number.
pmnumber_Callback([], [], handles)

% --------------------------------------------------------------------
function mnLoadmat_Callback(hObject, eventdata, handles)
% hObject    handle to mnLoadmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pbload_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to mnSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pbsave_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function mnDrawParticle_Callback(hObject, eventdata, handles)
% hObject    handle to mnDrawParticle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pbdraw_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnDrawParticleinUnitCell_Callback(hObject, eventdata, handles)
% hObject    handle to mnDrawParticleinUnitCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pbdrawinlattice_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function mnLoadASet_Callback(hObject, eventdata, handles)
% hObject    handle to mnLoadASet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PM_load_particles(handles)


% --------------------------------------------------------------------
function mnSaveAllParticles_Callback(hObject, eventdata, handles)
% hObject    handle to mnSaveAllParticles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
particles = evalin('base', 'SgAtoms');
PM_save_particles(particles)


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gisaxsleenew


function ed_dX_Callback(hObject, eventdata, handles)
% hObject    handle to ed_dX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_dX as text
%        str2double(get(hObject,'String')) returns contents of ed_dX as a double


% --- Executes during object creation, after setting all properties.
function ed_dX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_dX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_dY_Callback(hObject, eventdata, handles)
% hObject    handle to ed_dY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_dY as text
%        str2double(get(hObject,'String')) returns contents of ed_dY as a double


% --- Executes during object creation, after setting all properties.
function ed_dY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_dY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_qx_Callback(hObject, eventdata, handles)
% hObject    handle to ed_qx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_qx as text
%        str2double(get(hObject,'String')) returns contents of ed_qx as a double


% --- Executes during object creation, after setting all properties.
function ed_qx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_qx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_qy_Callback(hObject, eventdata, handles)
% hObject    handle to ed_qy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_qy as text
%        str2double(get(hObject,'String')) returns contents of ed_qy as a double


% --- Executes during object creation, after setting all properties.
function ed_qy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_qy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_qz_Callback(hObject, eventdata, handles)
% hObject    handle to ed_qz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_qz as text
%        str2double(get(hObject,'String')) returns contents of ed_qz as a double


% --- Executes during object creation, after setting all properties.
function ed_qz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_qz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function generate_2Dpixelcoordinates(handles)
    saxs.eng = str2double(get(handles.edit23, 'String'));
    saxs.psize = str2double(get(handles.edit25, 'String'));
    saxs.SDD = str2double(get(handles.edit26, 'String'));
    saxs.tiltangle = [0,0,0];
    saxs.waveln = 12.398/saxs.eng;
    saxs.center = [0, 0];
    saxs.ai = 0;
    saxs.tthi = 0;
    x = eval(get(handles.ed_dX, 'string'));
    y = eval(get(handles.ed_dY, 'string'));
    [X, Y] = meshgrid(x, y);
    [Q, SF] = pixel2qv([X(:)+saxs.center(1), Y(:)+saxs.center(2)], saxs);
    Q.SF = SF;
    Q.size = size(X);
    Q.coordinates = 'Pixel';
    Q.axeslabels = {'dX', 'dY'};
    Q.labelx = x;
    Q.labely = y;
    assignin('base', 'Q', Q);

function generate_2Dqcoordinates(handles)
    q_x = eval(get(handles.ed_qx, 'string'));
    q_y = eval(get(handles.ed_qy, 'string'));
    q_z = eval(get(handles.ed_qz, 'string'));
    arr = {q_x, q_y, q_z};
    dim = cellfun(@(x) numel(x), arr);
    %dim = [numel(qx), numel(qy), numel(qz)];
    if numel(dim == 1) < 1
        cprintf('*red', 'There should be at least one scalar among qx, qy and qz\n');
        return
    end
    ind = find(dim ~=1);
    Q.q_x = q_x;
    Q.q_y = q_y;
    Q.q_z = q_z;
    fn = fieldnames(Q);
    Q.coordinates = 'q';
    Q.axeslabels = fn(ind);
    Q.labelx = arr{ind(1)};
    Q.size = size(arr{ind(1)});
    if numel(ind)>1
        Q.labely = arr{ind(2)};
        x = Q.labelx;
        y = Q.labely;
        [X, Y] = meshgrid(x, y);
        Q.size = size(X);
        Q.(Q.axeslabels{1}) = X(:);
        Q.(Q.axeslabels{2}) = Y(:);
    end

    Qcell = {Q.q_x, Q.q_y, Q.q_z};
    numq = cellfun(@(x) numel(x), Qcell);
    scalarq = find(numq==1);
    arrayq = find(numq>1);
    if numel(arrayq)>=1
        for i=1:numel(scalarq)
            Qcell{scalarq(i)} = Qcell{scalarq(i)}*ones(size(Qcell{arrayq(1)}));
        end
    end
    Q.q_x = Qcell{1}(:);
    Q.q_y = Qcell{2}(:);
    Q.q_z = Qcell{3}(:);
    
    assignin('base', 'Q', Q);


function generate_3Dqcoordinates(handles)
    q_x = eval(get(handles.ed_qx, 'string'));
    q_y = eval(get(handles.ed_qy, 'string'));
    q_z = eval(get(handles.ed_qz, 'string'));
    arr = {q_x, q_y, q_z};
    sz = [numel(q_x), numel(q_y), numel(q_z)];
    [q_x, q_y, q_z] = meshgrid(q_x, q_y, q_z);
    Q.q_x = q_x;
    Q.q_y = q_y;
    Q.q_z = q_z;
    Q.coordinates = 'q';
    Q.axeslabels = {'q_x', 'q_y', 'q_z'};
    Q.labelx = arr{1};
    Q.labely = arr{2};
    Q.labelz = arr{3};
    Q.size = sz;
    
    assignin('base', 'Q', Q);

function ui_SimulationMode_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in ui_SimulationMode 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch lower(get(hObject, 'tag'))
    case 'rb_pixelcoordinates'
        generate_2Dpixelcoordinates(handles)
    case 'rb_qcoordinates'
        q_x = eval(get(handles.ed_qx, 'string'));
        q_y = eval(get(handles.ed_qy, 'string'));
        q_z = eval(get(handles.ed_qz, 'string'));
        arr = {q_x, q_y, q_z};
        dim = cellfun(@(x) numel(x), arr);
        %dim = [numel(qx), numel(qy), numel(qz)];
        if numel(dim == 1) < 1
            cprintf('*red', 'There should be at least one scalar among qx, qy and qz\n');
            return
        end
        ind = find(dim ~=1);
        if numel(ind) == 2
            generate_2Dqcoordinates(handles)
        elseif numel(ind) == 3
            generate_3Dqcoordinates(handles)
        end
    otherwise
        Q = [];
        assignin('base', 'Q', Q);
end


function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mncalFF_Callback(hObject, eventdata, handles)
ui_SimulationMode_SelectionChangedFcn(handles.ui_SimulationMode.SelectedObject, ...
    eventdata, handles)

Q = evalin('base', 'Q');
%if isempty(Q)
%    generate_2Dpixelcoordinates(handles)
%    Q = evalin('base', 'Q');
%end
model = evalin('base', 'SgAtoms');
if ~isfield(model{1}, 'position')
    model{1}.position = [0,0,0];
    model{1}.SOF = 1;
end
model(2:end) = [];
%obj = polyhedra(model{1}.shape, model{1}.edgelength);
%F = saxsPolyhedronAmp(Q.q_x(:), Q.q_y(:), Q.q_z(:), obj);
F = AnisoFormFactor(model, [Q.q_x(:), Q.q_y(:), Q.q_z(:)]);
F0 = AnisoFormFactor(model, [0, 0, 0]);
F = reshape(F, Q.size);
assignin('base', 'Fq', F/abs(F0));

function mnDrawFF_Callback(hObject,eventdata,handles)
Fq = evalin('base', 'Fq');
Q = evalin('base', 'Q');
if numel(Q.size) == 2
    figure;
    imagesc(Q.labelx, Q.labely, log10(abs(Fq).^2), [-5, 0]);
    axis image;
elseif numel(Q.size) == 3
    Iq = abs(Fq).^2;
    draw_3dmap(Iq, Q.labelx, Q.labely, Q.labelz)
end

function mnSaveFF_Callback(hObject,eventdata,handles)
% save
% This code checks if the user pressed cancel on the dialog.

Q = evalin('base', 'Q');
Fq = evalin('base', 'Fq');

[filename, pathname] = uiputfile( ...
   {'*.mat','MAT file (*.mat)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Save data');

% This code checks if the user pressed cancel on the dialog.

if isequal(filename,0) || isequal(pathname,0)
   disp('User pressed cancel')
   return
end
FN = fullfile(pathname, filename);
save(FN, 'Q', 'Fq');


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mncal_Callback(hObject, eventdata, handles)
% hObject    handle to mncal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ed_pdbfilename_Callback(hObject, eventdata, handles)
% hObject    handle to ed_pdbfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_pdbfilename as text
%        str2double(get(hObject,'String')) returns contents of ed_pdbfilename as a double


% --- Executes during object creation, after setting all properties.
function ed_pdbfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_pdbfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'enable', 'off')
