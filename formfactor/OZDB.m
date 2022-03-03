function out = OZDB(varargin)
global FitLee_helpstr
FitLee_helpstr = {'OZ and DB fit. ' ,...
'OZ = I0_oz./(1+xi_oz^2*q.^2)',...
'DB = I0_DB./(1+xi_DB^2*q.^2).^2',...
'Iq = OZ+DB+background',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'Byeongdu Lee (blee@anl.gov) ',...
'Ref: ',...
};

if numel(varargin) > 1
    p = varargin{1};
    q = varargin{2};
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
    bestP = [];
    bestP.I0_oz = 1;
    bestP.xi_oz = 20;
    bestP.I0_DB = 1;
    bestP.xi_DB = 20;
    % Need 4 parameters for background.
    bestP.poly1 = 0;
    bestP.poly2 = -2;
    bestP.poly3 = 0;
    bestP.poly4 = 0;
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
if iscell(q);
    if numel(q) > 1
        error('FitLee_schultzsphere.m is for fitting a set of data, for now')
    end
    q = q{1};
end
q = q(:);

% OZ(vari, q)
% I0 = vari(1);
% xi = vari(2);
% C = vari(3) : compressibility effect.


OZ = p.I0_oz./(1+p.xi_oz^2*q.^2);
out = OZ + p.I0_DB./(1+p.xi_DB^2*q.^2).^2;

back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
out = out + back;

if isnan(out)
    out = ones(size(out))*1E20;
end
