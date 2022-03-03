function [out, report] = FitLee_multivoigt(varargin)
% Usage:
% FitLee_multivoigt(N of peaks)
% FitLee('FitLee_multivoigt', bestP)
%
% multivoight function.
% if numel(varargin) == 1 and p is not a struct but a number, then
% it generate set of default parameters for using this function.
%
% Byeongdu Lee
if numel(varargin) > 1
    p = varargin{1};
    q = varargin{2};
    isini = 0;
elseif numel(varargin) == 1
    p = varargin{1};
    isini = 1;
    if ischar(p)
        out = FitLee_helpstr;
        return
    end
elseif numel(varargin) == 0
    p = 3;
    isini = 1;
end

%% initialize fit parameter bestP
if isini
    Nf = p;
    bestP = [];

    for i=1:Nf
        bestP.(['amp', num2str(i)]) = 1;
        bestP.(['center', num2str(i)]) = i;
        bestP.(['sigg', num2str(i)]) = 0.002;
        bestP.(['sigl', num2str(i)]) = 0.002;
    end
    % Need 4 parameters for background.
    bestP.poly1 = 1E-5;
    bestP.poly2 = -2;
    bestP.poly3 = -1;
    bestP.poly4 = 0.001;
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
%        back = zeros(size(x));
%        if isfield(p, 'poly')
%            back = p.poly(1)*q.^p.poly(2);
%            for i=numel(p.poly):-1:3
%                back = back + p.poly(i)*q.^(numel(p.poly)-i);
%            end
%        end
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;

predYs = {};
predY = zeros(size(q));
pnames = fieldnames(p);
np = 0;
for i=1:numel(pnames)
    if strfind(pnames{i}, 'amp')
        np = np + 1;
    end
end

for k=1:np
    pa = [p.(['amp', num2str(k)]), ...
        p.(['center', num2str(k)]), ...
        p.(['sigg', num2str(k)]), ...
        p.(['sigl', num2str(k)]), 0];
    pa = abs(pa);
    predYs{k} = abs(voigt(q, pa));
    predY = predY + predYs{k};
end
predY = predY + back;
out = predY;
if isnan(out)
    out = ones(size(out))*1E20;
end

if nargout == 2
    %x = 0:1:(max(p.r0, p.r1)+max(p.sig0, p.sig1)*10);
    for k=1:np
        pa = [p.(['amp', num2str(k)]), ...
            p.(['center', num2str(k)]), ...
            p.(['sigg', num2str(k)]), ...
            p.(['sigl', num2str(k)]), 0];
        xmin = pa(2)-max(pa(3:4))*10;
        xmax = pa(2)+max(pa(3:4))*10;
        x = linspace(xmin, xmax, 1000);
        parea = voigtarea(x, pa);
        W = voigtfwhm(pa(3), pa(4));
        fprintf('Area of the %ith peak is %0.3e\n', k, parea);
        fprintf('FWHM of the %ith peak is %0.3e\n', k, W);
    end    
    report = '';
end

%varargout{2} = fit;
%err = sum( (y-predY).^2);

