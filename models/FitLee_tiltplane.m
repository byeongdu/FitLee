function [out, report] = FitLee_tiltplane(varargin)
% This is to determine the normal vector of a plane rotating around the Z axis 
% from the interferometer sensor located at R and above the mirror (+Z).
global FitLee_helpstr
FitLee_helpstr = {'Position of zp ' ,...
'zp = (a*R*sin(th+th_sensor)+b*R*cos(th+th_sensor)) + off_set',...
'    zp: the distance to the mirror plane measured by the interferometer',... % 
'    Normal vector: n = (a, b, 1) of the rotating mirror plane',... % 
'    off_set: the Z position of the sensor (mm)',...
'    R: the radial distance of the sensor (mm)',...
'    The sensor is assumed to be located above the mirror (+Z)',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
};

if numel(varargin) > 1
    p = varargin{1};
    th = varargin{2};
    isini = 0;
elseif numel(varargin) == 1
    p = varargin{1};
    isini = 1;
elseif numel(varargin) == 0
    p = 1;
    isini = 1;
end

%% initialize fit parameter bestP
if isini
    Nf = p;
    bestP = [];
    bestP.a = 4E-5;
    bestP.b = 4E-5;
    bestP.th_sensor = 172; % degree
    bestP.R = 45; % mm unit
    bestP.off_set = 1E-4; % mm unit
    bestP.string = 'FitLee_tiltplane(figH)';
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
if iscell(th)
    if numel(th) > 1
        error('FitLee_tiltplane.m is for fitting a set of data, for now')
    end
    th = th{1};
end

out = p.off_set - (p.a*p.R*sind(th+p.th_sensor) + p.b*p.R*cosd(th+p.th_sensor));

if isnan(out)
    out = ones(size(out));
end

if nargout == 2
    fprintf('V angle = %0.3f urad.\n', atan(p.a)*1E6);
    fprintf('U angle = %0.3f urad.\n', atan(p.b)*1E6);
%   screw positions = [-8, 82, 172, 262] 
    screw_angle = [0, 90, 180, 270];
%   R of screws = 50*sqrt(2)/2;
    R_screw = 50*sqrt(2)/2;
    R_screw = 45;
    for i=1:numel(screw_angle)
        th = screw_angle(i);
        zpos = p.a*R_screw*sind(th+p.th_sensor) + p.b*R_screw*cosd(th+p.th_sensor);
        if zpos > 0
            direc = "Press down";
        else
            direc = "Loosen up";
        end
        fprintf("%s the %d screw by %0.3f um.\n", direc, th, abs(zpos)*1000);
    end
    report = '';
end