function y = SchultzCy_grad(varargin)
% y=SchultzSpSq(varargin)
% using shultz distribution and 2dpara crystal Sq.
% my LMA model.
% y=SchultzSphereFun(q, r, fwhm, df, dfa, dmax, numpnt)  
if isempty(varargin)
    y = {'I0', 2000, 0.8, 50000, 0;...
        'eden1', 1.5, 0.1, 3, 0;...
        'eden2', 0.7, 0.65, 0.75, 0;...
        'beta1', 1E-6, 1E-7, 1E-5, 0;...
        'beta2', 1E-8, 1E-9, 2E-8, 0;...
        'r', 20, 10, 50, 0;...
        'fwhm', 5, 0, 25, 0;...
        'HperR', 1.559, 0.1, 10, 0};
    return
else
    if iscell(varargin{1})
        param = varargin{1};
        cut = varargin{2};
        p=cell2struct(param(:,2)', param(:,1)',2);
        I0 = p.I0;
        eden1 = p.eden1;
        eden2 = p.eden2;
        beta1 = p.beta1;
        beta2 = p.beta2;
        r = p.r;
        fwhm = p.fwhm;
        HperR = p.HperR;

        if numel(varargin)<3
            strfactor = 'none';
            strparam = [];
        else
            strfactor = varargin{3};
            strparam = varargin{4};
            strparam=cell2struct(strparam(:,2)', strparam(:,1)',2);
        end
    else
        cut = varargin{1};
        param = varargin{2};
        p = cell2struct(param(:,2)', param(:,1)',2);
        I0 = p.I0;
        eden1 = p.eden1;
        eden2 = p.eden2;
        beta1 = p.beta1;
        beta2 = p.beta2;
        r = p.r;
        fwhm = p.fwhm;
        HperR = p.HperR;
%        strparam.df = df;
%        strparam.dfa = dfa;
%        strfactor = 'strfactor2';
        if numel(varargin) <4
            strfactor = 'none';
            strparam = [];
        else
            strfactor = varargin{3};
            strparam = varargin{4};
            strparam = cell2struct(strparam(:,2)', strparam(:,1)',2);
        end
    end
end
%I0, eden1, eden2, beta1, beta2, r, fwhm
if numel(varargin) < 6
    dmax = r + fwhm*3;
    numpnt = 20;
else
    dmax = varargin{6};
    numpnt = varargin{7};
end
% schultzspherefunction
% y = SchultzSphereFun(q, r,fwhm)
if fwhm > 0.001*r
    x = linspace(r/10, dmax, numpnt);x = x(:);
    nr = schultzdist(x, r, fwhm);nr = nr(:);
else
    x = r;
    nr = 1;
end
normv = sum(nr.*(pi*x.^3*HperR));
%y = pedersen_pDist(q, [x, nr], 1, df, dfa, 'strfactor_2Dpara');
%y = pedersen_pDist(q, [x, nr], 1, df, dfa, 'strfactor2');
%y = pedersen_pDist(q, [x, nr], 1, df, dfa, strfactor);
%y = cylinder_vert_clusters_pDist(q, [x, nr], HperR, I0, strfactor, strparam);
y = calc;
y = I0*y/normv;

function rtn = calc
    rtn = zeros(size(cut.af));
    for i = 1:length(x)
%         F = dwbacyl(qpa, q1z, q2z, q3z, q4z, Ti, Tf, Ri, Rf, xx(i), HperR*xx(i), [1, 0], 0);
        F = dwbacyl_grad(cut.qpa, cut.ai, cut.af, str2num(cut.waveln), eden1, eden2, beta1, beta2, x(i), x(i)*HperR, [1,0]);
        if ~strcmp(strfactor, 'none')
%            GG = strparam.G;
%            BB = strparam.B;
%            strparam.G = 0;
%            strparam.B = 0;
%            strfactor = 'strfactor_2Dpara';
            Struc = feval(strfactor, cut.qpa, strparam, x(i));  % for strfactor2
        else
            Struc = 1;
        end
        rtn = rtn + abs(F).^2.*Struc*nr(i);
    end
%    strfactor = 'strfactor_2Dpara';
%    Struc = feval('unifiedSAXS', cut.qpa, strparam);  % for strfactor2
%    rtn = rtn + Struc(:);
end
end