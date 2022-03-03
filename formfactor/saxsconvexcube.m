function [F, V] = saxsconvexcube(qx, qy, qz, parameter)
% function [F, V] = saxsconvexcube(qx, qy, qz, parameter)
% Calculating the scattering function of concave cube.
% paramter = [Cube edge length, Heigh of the concave pyramid]
Lcube = parameter(:, 1);
H_py = parameter(:, 2);
H2 = Lcube/2;
[Fcube, Vcube] = saxscube(qx, qy, qz, Lcube/2);

res = 1E5;
F = Fcube;
j = sqrt(-1);
for i=1:2:6 % this change is made to utilize the inversion symmetry 
    % F_spine (q; N,L,H_1,H_2 )=F_(p-con) (q_x,q_y,-q_z;N,L,H_1 ) e^(jq_z H_2 )
    % F_spine (R_i^T q;4,L,H_py,-L/2)
    switch i
        case 1
            R = eye(3); % for bottom 
        case 2
            R = rotate_around_vector([1, 0, 0], 180);
        case 3
            R = rotate_around_vector([1, 0, 0], 90);
        case 4
            R1 = rotate_around_vector([1, 0, 0], 90);
            R2 = rotate_around_vector([0, 0, 1], 180);
            R = R2*R1;
        case 5
            R1 = rotate_around_vector([1, 0, 0], 90);
            R2 = rotate_around_vector([0, 0, 1], 90);
            R = R2*R1;
        case 6
            R1 = rotate_around_vector([1, 0, 0], 90);
            R2 = rotate_around_vector([0, 0, 1], 270);
            R = R2*R1;
    end
    R = round(R*res)/res;
    Q = [qx, qy, qz]*R;
    q_x = Q(:,1);
    q_y = Q(:,2);
    q_z = Q(:,3);
    [Fi, Vpyi] = saxspyramid(q_x, q_y, q_z, [Lcube/2, H_py]);
    %Fi = Fi.*exp(-j*q_z*H2');  % this change is made to utilize the inversion symmetry 
    Fi = 2*real(Fi.*exp(-j*q_z*H2')); 
    F = F + Fi;
end
V = Vcube + 6*Vpyi(1,:);