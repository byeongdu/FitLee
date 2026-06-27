function [out, report] = FitLee_tiltplane(varargin)
% This is to determine the normal vector of a plane rotating around the Z axis 
% from the interferometer sensor located at R and above the mirror (+Z).
global FitLee_helpstr
FitLee_helpstr = {'Interferometer reading from a tilted mirror rotating about Z.' ,...
'',...
'$z_p(\theta) = offset - \left( a \cdot R \cdot \sin(\theta + \theta_{sensor}) + b \cdot R \cdot \cos(\theta + \theta_{sensor}) \right)$',...
'',...
'Angles in degrees. Mirror normal $\vec{n} = (a, b, 1)$ with $a, b \ll 1$.',...
'Sensor is above the mirror (+Z), measuring the distance to its plane.',...
'',...
'$Parameters$',...
'$\quad  a, b$ : tilt slopes of the mirror normal (dimensionless)',...
'$\quad  \theta_{sensor}$ : azimuth of the sensor relative to the rotation reference (deg)',...
'$\quad  R$ : sensor radial distance from rotation axis (mm)',...
'$\quad  offset$ : sensor Z-position (mm) when mirror is flat',...
'',...
'Diagnostic mode (nargout==2) reports tilt angles in $\mu rad$ and the',...
'screw adjustment (press/loosen) needed at each of 4 azimuth positions.',...
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
    Nf = p;
    bestP = [];
    bestP.a = 4E-5;
    bestP.b = 4E-5;
    bestP.th_sensor = 172; % degree
    bestP.R = 45; % mm unit
    bestP.offset = 1E-4; % mm unit
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

out = p.offset - (p.a*p.R*sind(th+p.th_sensor) + p.b*p.R*cosd(th+p.th_sensor));

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