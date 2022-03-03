function y = strfactor2(varargin)
% Percus-Yevick Hard sphere potential.
% parameter = [Rhs, vf]
% vf is volume fraction of hard sphere
% strfactor_2Dpara(q, rhs, sig)
if isempty(varargin)
    y = {'rhs', 2, 0.8, 10, 0;...
        'vf', 0.1, 0.0, 1, 0};
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
        rhs = p.rhs;
        vf = p.vf;
        qr = cut.qpa; %%% for temporary....
    else
        q = varargin{1};
        if isfield(q, 'qpa')
            qr = q.qpa;
%            disp('aaa')
        else
            qr = q;
%            disp('bbb')
        end
        if isstruct(varargin{2})  % strfactor_2Dpara(q, str_param, R)
            rhs = varargin{2}.rhs*varargin{3};
            vf = varargin{2}.vf;
        else                        % strfactor_2Dpara(q, rhs, sig)
            rhs = varargin{2};
            vf = varargin{3};
        end
    end
end

if (rhs <= 0)
    y = ones(size(q));
    return;
end
   
if vf == 0 
    y = ones(size(q));
    return
end
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
