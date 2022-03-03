function Sq = strfactor_unified_2D(varargin)
% structure factor for 2D paracrystal, positional disorder is described with gaussian function
%
% Ref : J. Phys.:Condens. Mater 9(1997) L125-L130. Printed in the UK.
% p = [sigma of gaussian, Dspacing];
if isempty(varargin)
    Sq = {'G', 1, 10, 1000, 0;...
        'B', 1, 10, 1000, 0;...
        'Rg', 150, 10, 500, 0;...
        'P', 4, 1, 4.5, 0;...
        'I_GB', 1, 0, 5, 0;...
        'sig', 5, 0.8, 1.2, 0;...
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
        G = p.G;
        B = p.B;
        Rg = p.Rg;
        P = p.P;
        sig = p.sig;
        I_GB = p.I_GB;
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
            p = varargin{2};
            G = p.G;
            B = p.B;
            Rg = p.Rg;
        I_GB = p.I_GB;
            P = p.P;
%            sig = p.sig
            d = p.dspacing*varargin{3};
            sig = d*p.sig;
%            d = varargin{2}.dspacing*varargin{3}
%            sig = varargin{2}.sig;
        elseif iscell(varargin{2})
            param = varargin{2};
            p=cell2struct(param(:,2)', param(:,1)',2);
            G = p.G;
            B = p.B;
            Rg = p.Rg;
            P = p.P;
        I_GB = p.I_GB;
            d = p.dspacing*varargin{3};
            sig = d*p.sig;
        else
            % strfactor_2Dpara(q, d, sig)
            d = varargin{2};
            sig = varargin{3};
        end
    end
end
if (d <= 0) || (sig<=0)
    Sq = ones(size(qr));
    return;
end
%    GausDis = exp(-sig^2*qr.^2);
%    Sq = (1-GausDis.^2)./(1-2*GausDis.*cos(qr*d) + GausDis.^2);
    Sq = strfact(qr, [d, sig]);
    unified = unifiedSAXS(qr, [G,B,Rg,P]);
    Sq = Sq(:) + I_GB*unified(:);
end
    function y = strfact(qr, p)
        rhs = p(1);
        vf = p(2);
        q = qr;
    alpha = (1 + 2*vf)^2/(1-vf)^4;
    beta = -6*vf*(1 + vf/2)^2/(1-vf)^4;
    gamma = vf*alpha/2;

    if length(rhs) > 1
        [R, q] = meshgrid(rhs, q);
        A = R.*q*2;
    else
        A = rhs*q.*2;  
    end

    GA = alpha*(sin(A)-A.*cos(A))./A.^2 + ...
        beta*(2.*A.*sin(A) + (2 - A.^2).*cos(A)- 2)./A.^3 + ...
        gamma*(-1*A.^4.*cos(A) + 4*((3*A.^2-6).*cos(A) + (A.^3 - 6 * A).*sin(A) + 6))./A.^5;
    y = 1./(1+24*vf*GA./A);    
    end

%else
%    Sq = ones(size(qr));
%end