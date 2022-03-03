function out = FitLee_trilamellaparacrystal(varargin)
    % Paracrystal theory equation for triblock copolymer.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'ABC triblock copolymer diffraction analysis. ' ,...
'Based on Paracrystal theory',...
'',...
'Lamella structure',...
'     ^                 D            ',...
'     |       |<----------------->|  ',...
'     |                  d3          ',...
'k2-|        d1      -----         ',...
'k1-|       ---        |    |         ',...
'     |       |  |  d2  |    |         ',...
'     |       |  |        |    |         ',...
'-------------------------------> z',...
'         T    H    T    F    T    H    ',...
'k1 : electron density difference between F and T phases',...
'k2 : electron density difference between H and T phases',...
't : interfacial thickness(A)',...
'd1 : thickness(A) of H lamella',...
'd2 : thickness(A) of T lamella',...
'D : domain size(A) = d1+2*d2+d3',...
'g : std of D / D',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    1. Y. Tanaka. Polymer J., 1999, 31, 11-2, 989–994. ',...
'    There is an error in the reference.',...
'    Correction made by Byeongdu Lee to the lattice factor'};

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
    p = 1;
    isini = 1;
end

%% initialize fit parameter bestP
if isini
    Nf = p;
    bestP = [];
    bestP.I0 = 1;
    bestP.k1 = 2.0;
    bestP.k2 = 0.55;
    bestP.t = 15;
    bestP.d1 = 110;
    bestP.sigx = 40;
    bestP.d2 = 95;
    bestP.sigy = 20;
    bestP.D = 400;
    bestP.g = 0.07;
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
k1 = p.k1;
k2 = p.k2;
sigx = p.sigx;
sigy = p.sigy;
D = p.D;
g = p.g;
d1 = p.d1;
d2 = p.d2;
sig = p.t/sqrt(2*pi);

y = linspace(0, 2*d1, 20);
x = linspace(0, 2*(d1+2*d2), 20);
[X, Y, Q] = ndgrid(x, y, q);
INT = (k1*sinc(Q.*X/2).*exp(-sig^2*Q.^2/2)-k2*sinc(Q.*Y/2).*exp(-sig^2*Q.^2/2)).^2.*...
    exp(-(X-(d1+2*d2)).^2/(2*sigx^2)).*exp(-(Y-d1).^2/(2*sigy^2));
Pq = trapz(trapz(INT, 1), 2)*(x(2)-x(1))*(y(2)-y(1));
F = abs(exp(-g^2*D^2*q.^2));
F2 = abs(exp(-g^2*D^2*q.^2/2));
%Z0 = (1-F2.^2)./(1-2*F2.^2.*cos(q*D)+F2.^2);
% correction by B. Lee
% In the original paper, the equation for Z is wrong.
Z = (1-F)./(1-2*F2.*cos(q*D)+F);
%Z = (Z0-(1-0.5*exp(-(g*D)^2*q.^2))+0.5)/0.5;
% compare Z0 and Z, where Z0 is the original Z in the paper..
assignin('base', 'Z', Z)
Iq = Pq(:).*Z./q.^2;
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
out = p.I0*Iq + back;

if isnan(out)
    out = ones(size(out))*1E20;
end

function y = sinc(x)
% sinc function
% sinc = sin(x)/x
% Byeongdu Lee
t = find(abs(x) > eps);
y = ones(size(x));
y(t) = sin(x(t))./x(t);
end
end
