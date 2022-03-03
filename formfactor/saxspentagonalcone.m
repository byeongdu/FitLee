function [F, V] = saxspentagonalcone(qx, qy, qz, LH)
% calculate using octant...
% L = R*sqrt((5-sqrt(5))/2)
L = LH(1);
R = L/sqrt((5-sqrt(5))/2);
H = LH(2);
F = zeros(size(qx));
Nbottom = 5;
for i=0:Nbottom-1
    M = [cos(2*pi/Nbottom*i),sin(2*pi/Nbottom*i),0;-sin(2*pi/Nbottom*i),cos(2*pi/Nbottom*i),0;0,0,1];
    q1 = M(1,1)*qx + M(1,2)*qy;
    q2 = M(2,1)*qx + M(2,2)*qy;
    [F1, V] = saxsoctant(q1, q2, qz, [R, L, H]);
    %F2 = saxsoctant(q1, q2, -qz, [R, L, H]);
    
%    t = find(abs(F1).^2>1E18);
%    if ~isempty(t)
%        i = numel(t);
%        fprintf('F value = %0.5f, at qx, qy, and qz = %f, %f, %f\n', F1(t(i)), q1(t(i)), q2(t(i)), qz(t(i)));
%    end
    %[q1(t(1:5)), q2(t(1:5)), qz(t(1:5))]*100000
%    t = find(abs(F2).^2>1E18);
    %[q1(t(1:5)), q2(t(1:5)), -qz(t(1:5))]*100000
%    if ~isempty(t)
%        i = numel(t);
%        fprintf('F value = %0.5f, at qx, qy, qz, qx-qy= %f, %f, %f, %e\n', F2(t(i)), q1(t(i)), q2(t(i)), -qz(t(i)), q1(t(i))-q2(t(i)))
%    end
    
    F = F+F1;
end
V = V*5;
%F = F;