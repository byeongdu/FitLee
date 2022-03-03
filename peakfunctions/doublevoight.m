function [predY, fit] = triplevoight(varargin)
% triplevoight function.
% two voight functions.
% The position of 2nd peak is p.center(2)*p.center(1). Therefore the
% p.center(2) should be q/q* where q* is the position of the first peak and
% q is the position of the second peak.
% Widths for the second peak are the same with those of the first one.
%
% 9/19/2014
% B. Lee

    p = [];
    if numel(varargin) > 1
        p = varargin{1};
        q = varargin{2};
    end

    if ~isstruct(p)
        Nf = 3;
        [predY, fit] = multivoight(Nf, []);
        return
    end


    % fitting codes .......................
%        back = zeros(size(x));
%        if isfield(p, 'poly')
%            back = p.poly(1)*q.^p.poly(2);
%            for i=numel(p.poly):-1:3
%                back = back + p.poly(i)*q.^(numel(p.poly)-i);
%            end
%        end
    back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;

    predYs = {};
    predY = zeros(size(q));
    pnames = fieldnames(p);
    np = 0;
    for i=1:numel(pnames)
        if strfind(pnames{i}, 'amp')
            np = np + 1;
        end
    end

    for k=1:np
        
        pa = [p.(['amp', num2str(k)]), ...
            p.(['center', num2str(k)]), ...
            p.(['sigg', num2str(k)]), ...
            p.(['sigl', num2str(k)]), 0];
        if k==1
            qstar = pa(2);
            widthsig = pa(3);
            widthlor = pa(4);
        end
        if k>1
            pa(2) = pa(2)*qstar;
            pa(3) = widthsig;
            pa(4) = widthlor;
        end
        pa = abs(pa);
        predYs{k} = abs(voigt(q, pa));
        predY = predY + predYs{k};
    end
    predY = predY + back;
end