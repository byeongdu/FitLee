function [Formfactor, name, pnames, pin] = SphConcShell(q, p, flag)
% Intensity of spherical concentric shells
% this form factor is a generalization of the shell form factor. 
% q is scattering vector
% p is input parameter
%   [number of all sphere, electron density 1, Radius 1, electron density2, Radius 2, ...., I0]

%global nint flag_aw

if nargin==2;
    lenp = length(p);

    if lenp > 1
        NumS = p(1);
    end
    
    rho = [];
    Radi = [];
    
    if (lenp - 2)/2 == NumS
    
        for i=2:2:(NumS*2)
            rho(i/2) = p(i);
            Radi(i/2) = p(1+i);
        end

        I0 = p(end);
        
        M2 = 0;
    
        for i = 2:NumS
            M2 = M2 + (4/3*pi*Radi(i)^3)*(rho(i)-rho(i-1));
        end
        
        M = rho(1)*(4/3*pi*Radi(1)^3) + M2;
        
        Formfactor = rho(1)*(4/3*pi*Radi(1)^3)*sphereAmp1(q, Radi(1));
        
        for i = 2:NumS
            Formfactor = Formfactor + (rho(i)-rho(i-1))*(4/3*pi*Radi(i)^3)*sphereAmp1(q, Radi(i));
        end
        
        Formfactor = (Formfactor/M).^2 * I0;
        
    else
        disp('Parameter Number is wrong')
    end
    
else
    
   Formfactor=[];
   name='Spherical Concentric Shell model';
   I0 = 1;
   
   if size(flag_aw,2) < 4
      flag_aw = [flag_aw zeros(1,4)];
   end;
   if flag_aw(1:4) ~= 'SCSM'
      nint=input('Number of Sphere:');
      flag_aw = str2mat('SCSM');
   end;
   pnames = str2mat('N-Sphere');
   for i = 1:nint
     pnames=str2mat(pnames,...
 	    sprintf('e-density_%d',i),...
	    sprintf('Radius_%d',i));
   end
   
   pnames=str2mat(pnames, sprintf('I0'));

    pin = [nint;zeros(2*nint,1)];
    pin = [pin;I0];
   
end   


function test = sphereAmp1(q, r)
% this Amplitude of Hard sphere type particle : spheretype(q, r, flaq)

test = 3*(sin(q*r) - r*q.*(cos(q*r)))./((q*r).^3);