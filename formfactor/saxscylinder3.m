function [F, Vcyl] = saxscylinder3(qp, p)
% This is SAXS for radial direction of a cylinder
%
% see also saxscylinder.m
% see also cylinder_amp.m, cylinder_vert_type.m
% p = [R, L]

R = p(:,1);
%H = p(:,2);
Vcyl = pi*R.^2;
F = repmat(Vcyl', numel(qp), 1).*besseljc(qp*R');
