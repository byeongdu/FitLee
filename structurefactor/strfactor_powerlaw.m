function y = strfactor_powerlaw(varargin)
% Percus-Yevick Hard sphere potential.
% parameter = [Rhs, vf]
% vf is volume fraction of hard sphere
if isempty(varargin)
    y = {'I0', 1, 0.001, 100, 0;...
        'power', 3, 0.0, 4.1, 0};
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
        I0 = p.I0;
        power = p.power;
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
            I0 = varargin{2}.I0;
            power = varargin{2}.power;
        else                        % strfactor_2Dpara(q, d, sig)
            I0 = varargin{2};
            power = varargin{3};
        end
    end
end

if (I0 <= 0)
    y = ones(size(q));
    return;
end
   
q = qr;
y = I0./q.^(-1*power);