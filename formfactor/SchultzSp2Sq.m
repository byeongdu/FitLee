function y = SchultzSp2Sq(varargin)
% y=SchultzSpSq(varargin)
% using shultz distribution and 2dpara crystal Sq.
% my LMA model.
% y=SchultzSphereFun(q, r, fwhm, df, dfa, dmax, numpnt)  
if isempty(varargin)
    y = {'I01', 2000, 0.8, 10000, 0;...
        'r', 6.559, 0.8, 12, 0;...
        'fwhm', 3, 0.8, 6, 0;...
        'RcRp1', 1, 0.5, 1.5, 0;...
        'RcbRc1', 1, 0.8, 1.2, 1;...
        'I02', 2000, 0.8, 10000, 0;...
        'r2perR', 0.6, 0.1, 1, 0;...
        'RcRp2', 1, 0.5, 1.5, 0;...
        'RcbRc2', 1, 0.8, 1.2, 1};
%        'df', 2, 0.8, 1.2, 0;...
%        'dfa', 0.3, 0.8, 1.2, 0;};
    return
else
    if numel(varargin)<3
        return
    end
    if iscell(varargin{1})
        param = varargin{1};
        cut = varargin{2};
        if numel(varargin) < 4
            strfactor = [];
            strparam = [];
        else
            strfactor = varargin{3};
            strparam = varargin{4};
            strparam=cell2struct(strparam(:,2)', strparam(:,1)',2);
        end
        p=cell2struct(param(:,2)', param(:,1)',2);
        r = p.r;
        fwhm = p.fwhm;
        q = cut;
    else
        q = varargin{1};
        r = varargin{2};
        fwhm = varargin{3};
        RcRp = varargin{4};
        RcbRc = varargin{5};
        df = varargin{6};
        dfa = varargin{7};
        p.r = r;
        p.fwhm = fwhm;
        p.RcRp = RcRp;
        p.RcbRc = RcbRc;
        strparam.df = df;
        strparam.dfa = dfa;
        strfactor = 'strfactor2'
    end
end

if numel(varargin) < 6
    dmax = r + fwhm*3;
    numpnt = 20;
else
    dmax = varargin{6};
    numpnt = varargin{7};
end
% schultzspherefunction
% y = SchultzSphereFun(q, r,fwhm)
x = linspace(r/10, dmax, numpnt);x = x(:);
nr = schultzdist(x, r, fwhm);nr = nr(:);
y = spheroid_vert_2clusters_pDist(q, [x, nr], p, strfactor, strparam);