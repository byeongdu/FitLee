function [out, report] = FitLee_cylinderCS_height(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz polydisperse stacked-cylinder (core+contrast cap) length form factor.' ,...
'Data should be plotted as $q$ vs $I(q) \cdot q$.',...
'',...
'$I(q) = I_0 \cdot \dfrac{\left<\,\left| H \cdot \mathrm{sinc}(q H/2) + \Delta\rho \cdot H_2 \cdot \mathrm{sinc}(q H_2/2) \cdot e^{-i q (H+H_2)/2} \right|^2\,\right>}{q - qz_0} + I_b$',...
'$\qquad$ models two stacked cylinder segments along the axis: a base segment of length $H$',...
'$\qquad$ (electron-density 1) and an end-cap segment of length $H_2$ (contrast $\Delta\rho$)',...
'$\qquad$ average over Schultz $n(H)$ with peak $H$ and FWHM $H \cdot sigH$',...
'$\qquad q$ is shifted by $qz_0$ before evaluation',...
'$\qquad I_b = poly1 \cdot q^{poly2} + poly3 \cdot q + poly4$',...
'',...
'Note: MATLAB''s sinc convention is $\sin(\pi x)/(\pi x)$ (preserved from original code).',...
'The phase factor $e^{-i q(H+H_2)/2}$ accounts for the displacement of the two segments.',...
'',...
'$Parameters$',...
'$\quad  I\_0$ : scale factor',...
'$\quad  qz0$ : $q$-offset (origin shift) (A$^{-1}$)',...
'$\quad  H$ : base-segment length, Schultz peak (A)',...
'$\quad  H2$ : end-cap segment length (A)',...
'$\quad  deltarho$ : electron-density contrast of cap relative to base',...
'$\quad  sigH$ : Schultz relative width for $H$ (FWHM = $H \cdot sigH$)',...
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
    bestP.H2 = 500;
    bestP.deltarho = 0.1;
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

x = linspace((p.H/10).^(1/1.5), (p.H2 + p.H*p.sigH*9).^(1/1.5), numpnt);
x = x.^1.5;

% n(1) = 1;n(2) = 1;
% k=3;
% while k <= numpnt+1
%     n(k) = n(k-1)+n(k-2);
%     k=k+1;
% end
% x = n(2:end)*p.H/10;
% x = x(:)';
nr = schultzdist(x, p.H, p.H*p.sigH);
j = sqrt(-1);
nr = nr(:)';
ILq = abs(x.*sinc(q*x/2)+p.deltarho*(p.H2).*sinc(q*(p.H2)/2).*exp(-j*q*(p.H+p.H2)/2)).^2.*repmat((sqrt(nr).*x).^2, length(q), 1);
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


