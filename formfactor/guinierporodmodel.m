function y = guinierporodmodel(varargin)
% Guinier Porod model
% B. Hammouda, J. Appl. Cryst.2010, 43, 716.
% Usage;
% Parameters: 
%   G : Guinier I0;
%   Rg : Rg of a particle;
%   d : Porod exponent;
%   s : dimension of the particle; 
% y = guinierporodmodel(q, G, d, Rg); without s means isometric particle.
% y = guinierporodmodel(q, G, d, Rg, s); % s is dimension of a particle
% For cylinder or lamella
% y = guinierporodmodel(q, G, d, 'L', 100, 'R', 5); % for a cylinder
% y = guinierporodmodel(q, G, d, 'L', 100, 'D', 10); % for a cylinder
% y = guinierporodmodel(q, G, d, 'W', 100, 'T', 5); % for a lamella
% Cylinder 
%  Rg = R/sqrt(2), Rg2 = (L^2/12+R^2/2)^(1/2);
% Lamella
%  Rg = T/sqrt(12), Rg2 = (W^2/12+T^2/12)^(1/2);

q = varargin{1};
if numel(varargin) <= 5
    G = varargin{2};
    d = varargin{3};
    Rg = varargin{4};
    s = 0;
end
if numel(varargin) == 5
    s = varargin{5};
end
if numel(varargin) == 7
    R = [];L = [];T=[];W=[];
    G = varargin{2};
    d = varargin{3};

    for i=4:2:6
        switch varargin{i}
            case 'D'
                R = varargin{i+1}/2; % cylinder diameter
            case 'R'
                R = varargin{i+1}; % cylinder radius
            case 'L'
                L = varargin{i+1}; % cylinder length
            case 'T'
                T = varargin{i+1}; % lamella thickness
            case 'W'
                W = varargin{i+1}; % lamella width
        end
    end
    if ~isempty(R) && ~isempty(L) % if cylinder
        Rg = R/sqrt(2);
        Rg2 = (L^2/12+R^2/2)^(1/2);
        s2 = 0;
        s = 1;
    end
    if ~isempty(T) && ~isempty(W) % if lamella
        Rg = T/sqrt(12);
        Rg2 = (W^2/12+T^2/12)^(1/2);
        s2 = 0;
        s = 2;
    end
end

if nargin < 5
    q1 = 1/Rg*(3*d/2)^(1/2);
    D = G*exp(-d/2)*(3*d/2)^(d/2)/Rg^d;
    qi0 = q>=q1;
    y=q;
    y(~qi0) = G*exp(-Rg^2/3*q(~qi0).^2);
    y(qi0) = D./q(qi0).^d;
    return
end
if d<=s
    error('d should be larger than s');
end
if s>=3
    error('s should be smaller than 3')
end
q1 = 1/Rg*((d-s)*(3-s)/2)^(1/2);
D = G/Rg^(d-s)*exp(-(d-s)/2)*((d-s)*(3-s)/2)^((d-s)/2);
qi0 = q>=q1;
y = q;
if nargin == 5
    y(~qi0) = G./q(~qi0).^s.*exp(-Rg^2/(3-s)*q(~qi0).^2);
    y(qi0) = D./q(qi0).^d;
    return
end
y(qi0) = D./q(qi0).^d;

q2 = ((s-s2)/(2/(3-s2)*Rg2^2-2/(3-s)*Rg^2))^(1/2);
G2 = G*exp(-q2^2*(Rg^2/(3-s)-Rg2^2/(3-s2)))*q2^(s2-s);
qi2 = q<q2;
qi1 = (q>=q2)&(q<q1);
y(qi1) = G./q(qi1).^s.*exp(-q(qi1).^2*Rg^2/(3-s));
y(qi2) = G2./q(qi2).^s2.*exp(-q(qi2).^2*Rg2^2/(3-s2));
