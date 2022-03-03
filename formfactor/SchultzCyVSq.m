function y = SchultzCyVSq(varargin)
% y=SchultzSpSq(varargin)
% using shultz distribution and 2dpara crystal Sq.
% my LMA model.
% y=SchultzSphereFun(q, r, fwhm, df, dfa, dmax, numpnt)  
if isempty(varargin)
    y = {'I0', 2000, 0.8, 1.2, 0;...
        'r', 6.559, 0.8, 1.2, 0;...
        'fwhm', 1.47, 0.8, 1.2, 0;...
        'HperR', 6.559, 0.8, 1.2, 0};
    
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
        strfactor = varargin{3};
        strparam = varargin{4};
        p=cell2struct(param(:,2)', param(:,1)',2);
        strparam=cell2struct(strparam(:,2)', strparam(:,1)',2);
        I0 = p.I0;
        r = p.r;
        fwhm = p.fwhm;
        HperR = p.HperR;
%        q = cut.q;
        q = cut;
%        strfactor = varargin{4};
    else
        q = varargin{1};
        r = varargin{2};
        fwhm = varargin{3};
        df = varargin{4};
        dfa = varargin{5};
        strparam.df = df;
        strparam.dfa = dfa;
        strfactor = 'strfactor2';
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
%y = pedersen_pDist(q, [x, nr], 1, df, dfa, 'strfactor_2Dpara');
%y = pedersen_pDist(q, [x, nr], 1, df, dfa, 'strfactor2');
%y = pedersen_pDist(q, [x, nr], 1, df, dfa, strfactor);
y = cylinder_vert_clusters_pDist(q, [x, nr], HperR, I0, strfactor, strparam);