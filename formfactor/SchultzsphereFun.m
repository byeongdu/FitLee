function [y, V2, V, zRg, betaq, mFq] = SchultzsphereFun(varargin)
% [y, V2, V, zRg, betaq, mFq] = SchultzsphereFun(varargin)
% analytical equation
%   without any varargin, the function prints default parameters.
% [y, V2, V, zRg] = SchultzSphereFun(q, r, fwhm)  
% outputs
%   y : integral_0^inf n(r)|F(q,r)|^2 dr / integral_0^inf n(r)*r^6 dr.
%       thus y(q = 0) = 1;
%   V2 : the mean value of V^2 or r^6
%       V2 = integral_0^inf f(r)*r^6 dr / integral_0^inf f(r) dr
%   V : the mean value of V or r^3
%       V2 = integral_0^inf f(r)*r^3 dr / integral_0^inf f(r) dr
%   zRg : z averaged Rg.
%   betaq : 
%  
% In order to calculate the mean V and V2, one should know the n th order
% average radii (or moments)
% Especially, the first, second, and third average radii are called the
% number, weight, and z-average radii, respectively.
% r_n = integral_0^inf f(r)*r^n dr / integral_0^inf f(r)*r^(n-1) dr.
% when f(r) is the Schultz function,
% r_n = (z+n)/(z+1)*r
% Thus, V = r_3*r_2*r_1. and V2 = r_6*r_5*r_4*r_3*r_2*r_1
if isempty(varargin)
    y = {'r', 50, 0.8, 1.2, 0;...
        'fwhm', 5, 0.8, 1.2, 0;};
    V2=0; V=0; zRg=0; betaq=0; mFq=0;
    return
else
    if numel(varargin)<3
        disp('number of input should be 3')
        return
    end
    if iscell(varargin{1})
        param = varargin{1};
        cut = varargin{2};
        var = varargin{3};
        p=cell2struct(param(:,2)', param(:,1)',2);
        r = p.r;
        fwhm = p.fwhm;
        q = cut(:,1); %%% for temporary....
    else
        q = varargin{1};
        r = varargin{2};
        fwhm = varargin{3};
    end
end
% schultzspherefunction
% y = SchultzSphereFun(q, r,fwhm)
if fwhm == 0
    mFq = sphereamp(q, r);
    V = 4*pi/3*r^3;
    V2 = V^2;
    zRg = sqrt(3/5)*r;
    betaq = 1;
    y = abs(mFq).^2/V2;
else
   	zp1=r*r/(fwhm*fwhm);

    Z = zp1 - 1;
    zRg = schultzRg(r, fwhm);
    V = 4*pi/3*(zp1+1)*(zp1+2)/zp1^2*r^3;
    V2 = (4*pi/3)^2*(Z+6)*(Z+5)*(Z+4)*(Z+3)*(Z+2)/zp1^5*r^6;

	b=2*q*r/zp1; % 2/alpha
	t=atan(b); % tan-1(2/alpha)
	den=sqrt(1+b.*b);
	dum1=(1.-(cos(zp1*t))./(den.^zp1));
	dum2=(b.^2).*(zp1+1).*zp1.*(1+(cos((zp1+2.).*t))./(den.^(zp1+2)))/4;
	dum3=-zp1.*b.*(sin((zp1+1).*t))./(den.^(zp1+1));
	m4sq=(b.^6).*(zp1+5).*(zp1+4).*(zp1+3).*(zp1+2).*(zp1+1)*zp1/288;
	y=(dum1+dum2+dum3)./(m4sq);
    if nargout >= 5
        a = zp1./(q*r);
        at2a = atan(2./a);
        a4 = 4 + a.^2;
        %t1G1 = a.^(-zp1);
        %t2G1 = -a4.^(-zp1/2).*cos((zp1*at2a));
        %t3G1 = (zp1+1)*zp1*(a.^(-(zp1+2))+a4.^(-(zp1+2)/2).*cos((zp1+2).*at2a));
        %t4G1 = -2*zp1.*a4.^(-(zp1+1)/2).*sin((zp1+1).*at2a);
        t1G1 = 1;
        t2G1 = -(a./sqrt(a4)).^(zp1).*cos((zp1*at2a));
        t3G1 = (zp1+1)*zp1*(a.^(-(2))+(a./sqrt(a4)).^(zp1).*a4.^(-1).*cos((zp1+2).*at2a));
        t4G1 = -2*zp1.*(a./sqrt(a4)).^(zp1).*a4.^(-1/2).*sin((zp1+1).*at2a);
        G1 = t1G1+t2G1+t3G1+t4G1;
        y = 8*pi.^2*r^6*zp1.^(-6).*a.^(6).*G1;
        %y = 8*pi.^2*r^6*zp1.^(-6).*a.^(zp1+6).*G1;
        y = y/V2;
        
        at1a = atan(1./a);
        t1G2 = sin(zp1.*at1a);
        t2G2 = -(zp1).*(1+a.^2).^(-1/2).*cos((zp1+1).*at1a);
        G2 = t1G2+t2G2;
        betaq = 2.*(a./(1+a.^2)).^zp1.*G2.^2./G1;
        
        mFq = 4*pi*r^3*(a./sqrt(1+a.^2)).^(zp1).*a.^3.*zp1.^(-3).*G2;
        y2 = abs(mFq/V).^2;
        %q*r, y
        t = (q*r < 1E-5) | isnan(y);
        y2(t) = 1;
        y(t) = y2(t);
        mFq(t) = V;
        betaq(t) = 0;
            
        %betaq = 2.*a.^zp1.*(1+a.^2).^(-zp1);
        %betaq = G2.^2./G1;
        %t = find(q<pi/zRg);
        %k = find(betaq(t)==0);
        %betaq(k) = (Z+3).*(Z+2).*(Z+1)./((Z+6).*(Z+5).*(Z+4));
    end
end

