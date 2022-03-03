function y=pseudovoigt(x,p)
% voigt     : Voigt
% function y=pseudovoigt(x,p)
% p = [Area, Center, Gaussian_FWHM, Lorentz_FWHM]
% Gaussian_FWHM = 2*Gaussian_sig*sqrt(2*log(2))
% Lorentzian_FWHM = 2*Lorentzian_gamma;
% Description:  Wikipedia
x = x(:);
%fG = 2*p(3)*sqrt(2*log(2));
%fL = 2*p(4);
fG = p(3);
fL = p(4);
f = (fG.^5 + 2.69269*fG.^4.*fL + 2.43843*fG.^3.*fL.^2 + ...
    4.47163*fG.^2.*fL.^3 + 0.07842*fG.*fL.^4 + fL.^5).^(1/5);
eta = 1.36603*(fL./f) - 0.47719*(fL./f).^2 + 0.11116*(fL./f).^3;
y = eta*L(x-p(2), f) + (1-eta)*G(x-p(2), f);
y = y*p(1);
    function y = L(x, gamma)
        y = gamma./(pi*(x.^2+gamma^2));
    end
    function y = G(x, sig)
        y = exp(-x.^2/(2*sig^2))/(sig*sqrt(2*pi));
    end
end