function Sq = strfactor_2Dpara(varargin)
% structure factor for 2D paracrystal, positional disorder is described with gaussian function
%
% Ref : J. Phys.:Condens. Mater 9(1997) L125-L130. Printed in the UK.
% p = [sigma of gaussian, Dspacing];
if isempty(varargin)
    Sq = {'sig', 5, 0.8, 1.2, 0;...
        'dspacing', 50, 0.8, 1.2, 0};
    return
else
    if numel(varargin)<3
        return
    end
    if iscell(varargin{1})
        param = varargin{1};
        cut = varargin{2};
        var = varargin{3};
        p=cell2struct(param(:,2)', param(:,1)',2);
        sig = p.sig;
        d = p.dspacing;
        qr = cut(:,1); %%% for temporary....
    else
        q = varargin{1};
        if isfield(q, 'qpa')
            qr = q.qpa;
        else
            qr = q;
        end
        if isstruct(varargin{2})  % strfactor_2Dpara(q, str_param, R)
            d = varargin{2}.dspacing*varargin{3};
            sig = d*varargin{2}.sig;
        else                        % strfactor_2Dpara(q, d, sig)
            d = varargin{2};
            sig = varargin{3};
        end
    end
end
if (d <= 0) || (sig<=0)
    Sq = ones(size(qr));
    return;
end
%sig
%d
    GausDis = exp(-sig^2*qr.^2);
    Sq = (1-GausDis.^2)./(1-2*GausDis.*cos(qr*d) + GausDis.^2);
    Sq = Sq(:);
%else
%    Sq = ones(size(qr));
%end