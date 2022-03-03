function [y, Fq, Iq] = lam_twophase(varargin)
% y = lam_twophase(q, [center_center_distance, volf_of_A], delta_rho, domainsize, microstrain, DW)

q = varargin{1};
d = varargin{2}(1);
vf = varargin{2}(2);
drho = varargin{3};
if numel(varargin)<4
    domainsize = 5000;
    microstrain = 0.02;
    DW = 0.1;
    ispowder = 1;
else
    domainsize = varargin{4};
    microstrain = varargin{5};
    DW = varargin{6};
    ispowder = 1;
end
thick = [d*vf, d*(1-vf)];
density = [0, drho];
[y, Fq, Iq] = lamella(q, thick, density, domainsize, microstrain, DW, ispowder);