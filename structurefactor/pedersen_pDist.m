function y = pedersen_pDist(x, pDist, I0, strfactor, strp)
% [y, name, pnames, pin] = pedersen_pDist(x, pDist, p)
% Pedersen formalism to calculate scattering intensity from a given distribution function.
% Local monodisperse approximation assumed.
% Formfactor assumed as sphere scattering.
% Structure factor : Hard sphere
% pDist = [x, Dist];
% p = [I0, volfrac, back];

if nargin > 3;
    xx = pDist(:, 1);
    gDist = pDist(:, 2);
    df = strp.dspacing;
    dfa = strp.sig;
    x = x(:);
    normFactor = sum(gDist.*xx.^6);  % for sphere
%    normFactor = sum(gDist.*xx.^4)*40^2;  % for cylinder
    y = zeros(size(x));
    for i = 1:length(xx)
        %Form = spheretype(x, xx(i))*xx(i)^6;  % let's discribe using volume distribution... because I don't want exaggerate small size particle...
        if isfield(x, 'q')
            q = x.q;
        else
            q = x;
        end
        Form = xx(i).^6*spheretype(q, xx(i));  % for sphere
%        Form = abs(cylinder_vert_type(x, 0, [xx(i), 40])).^2; % for cylinder
%        Struc = feval(strfactor, x, dfa*df*xx(i), xx(i)*df); % for strfactor_2Dpara
        Struc = feval(strfactor, x, df*xx(i), dfa);  % for strfactor2
        y = y + Form.*Struc.*gDist(i);
    end
    y = y*I0/normFactor;

else
	y=[];
	name='Gauss distr Sphere';
end