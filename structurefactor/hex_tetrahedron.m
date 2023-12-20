function Sq = hex_tetrahedron(varargin)
% y = hex_twophase(q, [lattice_const, L, th0], delta_rho, domainsize, microstrain, DW)
q = varargin{1};
a = varargin{2}(1);
L = varargin{2}(2);
th0 = varargin{2}(3)*pi/180;
domainsize = 5000;
microstrain = 0.02;
DW = 0.1;
gaussianwidth = 0.00025;
if numel(varargin) > 2
    domainsize = varargin{3};
end
if numel(varargin) > 3
    gaussianwidth = varargin{4};
end
if numel(varargin) > 4
    microstrain = varargin{5};
end
if numel(varargin) > 5
    DW = varargin{6};
end

wavelength = 1; % for peakwidth calculation only.
% thick = zeros(size(relcomp));
% maxh = 5;

DW = DW*a;

as = [2*pi,    2*pi/sqrt(3),   0];
bs = [0,    4*pi/sqrt(3),   0];
cs = [0,    0,   1];

hklmat = load('hklhex.mat');
hkl = hklmat.hkl;
hklhex = hklmat.hklhex;
m = hklmat.m;

Q0 = hkl'*[as;bs;cs]/(a);
%Q0 = Q0';
qx0 = Q0(:,1);
qy0 = Q0(:,2);
%qz0 = Q0(:,3);
qr = sqrt(qx0.^2+qy0.^2);
Fq2 = hex_tetrahedron_Fq_squared(L/a, th0);

y = zeros(size(q));
y = y(:);
Zq = zeros(size(y));


omega = 2*pi;
dim_cell = 2;
A = (1*sqrt(3)/2)*a^2;
cVol = A;
factor = 1/omega*(2*pi)^3/cVol;
%factor = 1/omega*(2*pi)^3/cellinfo.Vol;
Lorentzfactor = qr.^(dim_cell-1);

% final 
%Iqa = factor*Iq./Pq./Lorentzfactor;



for h = 1:numel(qr)
    qh = qr(h);
    w = peakwidth(qh, wavelength, domainsize, microstrain);
    t = pseudovoigt(q, [1, qh, gaussianwidth, w]);%trapz(q, t)
    Zq = Zq + t(:)*factor*Fq2(h)/Lorentzfactor(h);
end
Gq = exp(-(DW)^2*q.^2);
Sq = 1+(Zq - 1).*Gq;


function Fq_squared = hex_tetrahedron_Fq_squared(La_ratio, th0)
%load hklhex.mat;
m = m';
% hklhex and hkl

% % for the lattice constant 1:
% as = [6.2832,    3.6276,   0];
% bs = [0,    7.2552,   0];
% cs = [0,    0,   1];
Q = hklhex'*[as;bs;cs];
%qx = Q(:,1);
%qy = Q(:,2);
%qz = Q(:,3);

% Inversion rotation matrix R; 
% when a object rotates by R in the real space,
% the same effect can be computed by rotating inv(R) of Q coordinates..
R = [cos(th0), sin(th0), 0;-sin(th0), cos(th0),0;0,0,1];
qtmp = Q*R';
F = saxstetrahedron(qtmp(:,1),qtmp(:,2),qtmp(:,3), La_ratio);
n0 = 1;
Fq = zeros(11, 1);
mcumsum = cumsum(m);
for i=1:numel(m)
    n = n0:mcumsum(i);
    Fq(i) = sum(F(n));
    n0=mcumsum(i)+1;
end
Iq = abs(Fq).^2./m;
% form factor
Qreduced = hkl'*[as;bs;cs];
%Qreduced = Qreduced';
qx_reduced = Qreduced(:,1);
qy_reduced = Qreduced(:,2);
qz_reduced = Qreduced(:,3);

qr_reduced = sqrt(qx_reduced.^2+qy_reduced.^2);
q0 = zeros(size(qr_reduced));
Q = [qr_reduced, qz_reduced];
Pq = q0;
k=0;
for ang0=-pi/3+pi/180:pi/180:pi/3
    R = [cos(ang0), sin(ang0);-sin(ang0), cos(ang0)];
    q_reduced = Q*R';
    F = saxstetrahedron(q_reduced(:,1),q_reduced(:,2),q0, La_ratio);
    Pq = Pq+abs(F).^2;
    k = k + 1;
end
Pq = Pq/k;
Fq_squared = Iq./Pq;
end
end
