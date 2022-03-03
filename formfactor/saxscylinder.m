function [F, Vcyl] = saxscylinder(qx, qy, qz, p)
% this is for a vertical cylinder. 
% see also cylinder_amp.m, cylinder_vert_type.m
% p = [R, L]

qp = sqrt(qx.^2+qy.^2);
R = p(:,1);
H = p(:,2);
Vcyl = pi*R.^2.*H;
F = 2*repmat(Vcyl', numel(qx), 1).*besseljc(qp*R').*sinc(qz*H'/2);
