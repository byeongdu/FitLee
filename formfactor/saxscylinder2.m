function [F, Vcyl] = saxscylinder2(qx, qy, qz, p)
% this is for a vertical cylinder. 
% where the origin of coordinate is at the bottom of a cylinder.
% see also saxscylinder.m
% see also cylinder_amp.m, cylinder_vert_type.m
% p = [R, L]

qp = sqrt(qx.^2+qy.^2);
R = p(:,1);
H = p(:,2);
Vcyl = pi*R.^2.*H;
F = 2*repmat(Vcyl', numel(qx), 1).*besseljc(qp*R').*sinc(qz*H'/2).*exp(-sqrt(-1)*qz*H'/2);
