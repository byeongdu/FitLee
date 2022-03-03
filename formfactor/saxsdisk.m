function F = saxsdisk(qp, qz, R, H)
% SAXS FF of a disk whose center of mass at [0,0,0] and normal vector of
% the round circle face is z.
%
% see also cylinder_amp.m, cylinder_vert_type.m
%
%sizeqp = size(qp);
%sizeqp = size(qp);
if (size(qz,1) ~=0) && (size(qz,2) == numel(R)) && (numel(H)==1) &&(numel(qp)==size(qz,1))
    qp = repmat(qp, 1, size(qz,2));
    R = repmat(R, size(qz, 1), 1);
    Vcyl = pi*R.^2*H;
    F = 2*Vcyl.*besseljc(qp.*R).*sinc(qz*H/2);
    return
end

qp = qp(:);
qz = qz(:);
R = R(:)';
H = H(:)';

numelqp = numel(qp);
numelqz = numel(qz);
if (numelqp ~= numelqz)
    if (numelqp == 1)
        qp = repmat(qp, numel(qz), 1);
    elseif (numelqz == 1)
        qz = repmat(qp, numel(qp), 1);
    else
        error('numel(qp) should be identical to numel(qz), or either of them should be 1')
    end
end

numelR = numel(R);
numelH = numel(H);
numelqp = numel(qp);
numelqz = numel(qz);
if (numelR ~= numelH)
    if (numelR == 1)
        R = repmat(R, 1, numel(H));
    elseif (numelH == 1)
        H = repmat(H, 1, numel(R));
    else
        error('numel(qp) should be identical to numel(qz), or either of them should be 1')
    end
end

Vcyl = ones(size(qp))*(pi*R.^2.*H);
F = 2*Vcyl.*besseljc(qp*R).*sinc(qz*H/2);