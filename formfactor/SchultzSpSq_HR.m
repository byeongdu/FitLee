function y = SchultzSpSq_HR(varargin)
% y=SchultzSpSq(varargin)
% using shultz distribution and 2dpara crystal Sq.
% my LMA model.
% y=SchultzSphereFun(q, r, fwhm, df, dfa, dmax, numpnt)  
if isempty(varargin)
    y = {'I0', 2000, 1000, 100000, 0;...
        'r', 25.559, 0.8, 60, 0;...
        'rfwhm', 16, 0.8, 30, 0;...
        'h', 41, 0.8, 50, 0;...
        'hfwhm', 0, 1, 25, 0;...
        'Rshell', 0.3, 0.1, 0.9, 1;...
        'rhoShell', 0.8, 0.2, 2, 1;...
        'RcbRc', 1, -0.5, 1, 1};
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
        if numel(varargin) > 3
            strparam = varargin{4};
            strparam=cell2struct(strparam(:,2)', strparam(:,1)',2);
        else
            strparam = [];
        end
        p=cell2struct(param(:,2)', param(:,1)',2);
        r = p.r;
        rfwhm = p.rfwhm;
        h = p.h;
        hfwhm = p.hfwhm;
        q = cut;
%        strfactor = varargin{4};
    else
        q = varargin{1};
        r = varargin{2};
        rfwhm = varargin{3};
        h = varargin{4};
        hfwhm = varargin{5};
        RcbRc = varargin{6};
        df = varargin{7};
        dfa = varargin{8};
        p.r = r;
        p.rfwhm = rfwhm;
        p.h = h;
        p.hfwhm = hfwhm;
        p.RcbRc = RcbRc;
        strparam.df = df;
        strparam.dfa = dfa;
        strfactor = 'strfactor2';
    end
end

if numel(varargin) < 9
    rdmax = r + rfwhm*5;
    hdmax = h + hfwhm*5;
    if hfwhm ~= 0
        numpnt = 7;
    else
        numpnt = 10;
    end
else
    rdmax = varargin{9};
    hdmax = varargin{10};
    numpnt = varargin{11};
end
% schultzspherefunction
% y = SchultzSphereFun(q, r,fwhm)
rx = linspace(r/10, rdmax, numpnt);rx = rx(:);
nr = schultzdist(rx, r, rfwhm);nr = nr(:);
if hfwhm ~= 0
    hx = linspace(h/10, hdmax, numpnt);hx = hx(:);
    nh = schultzdist(hx, h, hfwhm);nh = nh(:);
else
    hx = h;
    nh = 1;
end
maxnr = max(nr);t = find(nr < 0.001*maxnr);nr(t) = []; rx(t) = [];
maxnh = max(nh);t = find(nh < 0.001*maxnh);nh(t) = []; hx(t) = [];
p.hx = hx;
p.nh = nh;
%y = pedersen_pDist(q, [x, nr], 1, df, dfa, 'strfactor_2Dpara');
%y = pedersen_pDist(q, [x, nr], 1, df, dfa, 'strfactor2');
%y = pedersen_pDist(q, [x, nr], 1, df, dfa, strfactor);
y = spheroid_vert_clusters_pDist(q, [rx, nr], p, strfactor, strparam);