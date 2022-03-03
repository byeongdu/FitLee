function Sq = baxter(varargin)
if isempty(varargin)
    Sq = {'volf', 0.1, 0, 1, 0;...
        'tau', 5, 0, 20, 0;...
        'Rratio', 2, 0, 5, 0};
    return
end

if iscell(varargin{1})
    param = varargin{1};
    cut = varargin{2};
    p=cell2struct(param(:,2)', param(:,1)',2);
    volf = p.volf;
    tau = p.tau;
    Rratio = p.Rratio;
    q = cut.qpa; %%% for temporary....
    r = varargin{3};
else
%if numel(varargin) == 2; % [x, parameter]
    q = varargin{1};
    p = varargin{2};

    if isfield(q, 'qpa')
        q = q.qpa;
    else
        %q = q;
    end
    if iscell(p)
        p=cell2struct(p(:,2)', p(:,1)', 2);
        volf = p.volf;
        tau = p.tau;
        Rratio = p.Rratio;
        r = varargin{3};
    end
    if isstruct(varargin{2})  % strfactor_2Dpara(q, str_param, R)
        volf = p.volf;
        tau = p.tau;
        Rratio = p.Rratio;
        r = varargin{3};
    end                        % strfactor_2Dpara(q, d, sig)
    if isnumeric(varargin{2})
        volf = p(1);
        tau = p(2);
        Rratio = p(3);
        r = p(4);
    end
end
r = r*Rratio;

g = volf*(1+0.5*volf)/(3*(1-volf)^2);
e = tau + volf/(1-volf);
l = (6/volf)*(e-sqrt(e^2-g));
m = l*volf*(1-volf);
b = -(3*volf*(2+ volf)^2 - 2*m*(1+7*volf + volf^2) + m^2*(2+volf))/(2*(1-volf)^4);
a = (1+2*volf-m)^2*(1 - volf)^-4;
x = q*r;
c1 = a*x.^3.*(sin(x) - x.*cos(x)) + b*x.^2.*(2*x.*sin(x) - (x.^2-2).*cos(x) - 2);
c2 = 0.5*volf*a*((4*x.^3-24*x).*sin(x) - (x.^4 - 12*x.^2 + 24).*cos(x) + 24);
c = -24*volf*x.^-6.*(c1+c2);
c3 = -2*volf^2*l^2*(1-cos(x)).*x.^-2 + 2*volf*l*x.^-1.*sin(x);
cq = c + c3;
Sq = 1./(1-cq);
