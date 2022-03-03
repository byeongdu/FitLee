function [out, report] = FitLee_superposition_two(varargin)
    % Superposition of two curves, fit the ratio of each.
FitLee_helpstr = {'Fit coefficient of two states to the mixed state. ' ,...
'I(q) = a*I_a(q) + b*I_b(q) + background',...
'    I_a(q) is the curve of the state a',...
'    I_b(q) is the curve of the state b',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'Byeongdu Lee (blee@anl.gov)',...
};

if numel(varargin) > 1
    p = varargin{1};
    q = varargin{2};
    isini = 0;
elseif numel(varargin) == 1
    p = varargin{1};
    isini = 1;
    if ishandle(p)
        guiadd(p);
        return
    end
    if ischar(p)
        out = FitLee_helpstr;
        return
    end
elseif numel(varargin) == 0
    p = 1;
    isini = 1;
end

%% initialize fit parameter bestP
if isini
    bestP = [];
    bestP.a = 1;
    bestP.b = 1;
    % Need 4 parameters for background.
    bestP.poly1 = 0;
    bestP.poly2 = -2;
    bestP.poly3 = 0;
    bestP.poly4 = 0;
    bestP.string = 'FitLee_superposition_two(figH)';
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
if iscell(q)
    if numel(q) > 1
        error('FitLee_schultzsphere.m is for fitting a set of data, for now')
    end
    q = q{1};
end


try
    Ia = evalin('base', 'I_a');
    I_a = interp1(Ia(:,1), Ia(:,2), q);
    Ib = evalin('base', 'I_b');
    I_b = interp1(Ib(:,1), Ib(:,2), q);
catch
    I_a = zeros(size(q));
    I_b = zeros(size(q));
end


q = q(:);
Iq = p.a*I_a + p.b*I_b;
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
%out = Poq + pnumberfraction*r_e^2*Angstrom2Centimeter^2*(p.delta_rho0^2*Pq1.*Sq+ p.Nratio*p.delta_rho1^2*Pq2)+back;
out = [Iq + back, Iq, back];
if isnan(out)
    out = ones(size(out));
end

if nargout == 2
    x = p.a/(p.a+p.b);
    y = p.b/(p.a+p.b);
    fprintf('Fractional contribution of a is %0.3f.\n', x);
    fprintf('Fractional contribution of b is %0.3f.\n', y);
    report = '';
end
function guiadd(figH)
    uph = findobj(figH, 'type', 'uipanel');
    pos = get(figH, 'position');
    hFigHeight = pos(end);
        uicontrol(...
          'style'                     , 'text', ...
          'backgroundcolor'           , [0,0.5,0.5], ...
          'foregroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'Left', ...
          'parent'                    , uph, ...
          'string'                  ,   'Load Iqs',...
          'position'                  ,[1, hFigHeight-210,100,20], ...
          'tag'                       , 'text_Pshape');
    popup = uicontrol(...
          'parent', uph,...
          'horizontalalignment'       , 'Left', ...
          'Style', 'pushbutton',...
          'String', 'Load I_a(q)',...
          'tag', 'Iqa',...
          'Position', [70, hFigHeight-215, 130, 25],...
          'Callback', @load_Iqa);
    popup2 = uicontrol(...
          'parent', uph,...
          'horizontalalignment'       , 'Left', ...
          'Style', 'pushbutton',...
          'String', 'Load I_b(q)',...
          'tag', 'Iqb',...
          'Position', [70, hFigHeight-240, 130, 25],...
          'Callback', @load_Iqa);
function load_Iqa(varargin)
    obj = varargin{1};
    [filename, pathname] = uigetfile('*.txt; *.dat; *.*', 'Pick your file');
    if isequal(filename,0) || isequal(pathname,0)
        disp('User pressed cancel')
        return
    else
        disp(['User selected ', fullfile(pathname, filename)])
    end
    filename = fullfile(pathname, filesep, filename);
    [~, data] = hdrload(filename);
    switch lower(get(obj, 'tag'))
        case 'iqa'
            assignin('base', 'I_a', data)            
        case 'iqb'
            assignin('base', 'I_b', data)            
    end
