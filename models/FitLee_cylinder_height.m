function [out, report] = FitLee_cylinder_height(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz polydisperse coreshell cylinder radial direction fit. ' ,...
    'Data needs to be plotted as q vs Iq*q',...
'I(q) = I0*P(q; r0, sig0)*S(q) + background',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
};

if numel(varargin) > 1
    p = varargin{1};
    q = varargin{2};
    isini = 0;
elseif numel(varargin) == 1
    p = varargin{1};
    isini = 1;
    if ischar(p)
        out = FitLee_helpstr;
        return
    end
elseif numel(varargin) == 0
    p = 1;
    isini = 1;
end

%% initialize fit parameter bestP
if isini
    bestP = [];
    bestP.I0 = 1;
    bestP.qz0 = 0.02;
    bestP.H = 500;
    bestP.sigH = 0.1;

    bestP.poly1 = 0;
    bestP.poly2 = -2;
    bestP.poly3 = 0;
    bestP.poly4 = 0;
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
if iscell(q)
    if numel(q) > 1
        error('FitLee_schultzsphere.m is for fitting a set of data, for now')
    end
    q = q{1};
end

% radial direction
q = q(:)-p.qz0;
numpnt = 55;
x = linspace((p.H/10).^(1/1.5), (p.H + p.H*p.sigH*9).^(1/1.5), numpnt);x = x(:)';
x = x.^1.5;
% n(1) = 1;n(2) = 1;
% k=3;
% while k <= numpnt+1
%     n(k) = n(k-1)+n(k-2);
%     k=k+1;
% end
% x = n(2:end)*p.H/100;
% x = x(:)';



nr = schultzdist(x, p.H, p.H*p.sigH);
nr = nr(:)';
ILq = sinc(q*x/2).^2.*repmat((sqrt(nr).*x).^2, length(q), 1);
ILq = sum(ILq, 2);


back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
out = p.I0*ILq./q + back;
if isnan(out)
    out = ones(size(out))*1E20;
end

if nargout == 2
    maxR =p.H;
    sigma = p.sigH;
    x = linspace(maxR/10, maxR + maxR*sigma*9, 100);
    nr = schultzdist(x, maxR, maxR*sigma);

    
    figure;
    plot(x, nr);xlabel('Height (A)');ylabel('n(R)')
    report = '';
end


