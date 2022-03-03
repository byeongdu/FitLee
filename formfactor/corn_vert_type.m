function [Formfactor, name, pnames, pin] = corn_vert_type(qp, qz, p)
% this is Amplitude for cone type particle : corn_vert_type(qp, qz, p)
% cone is assumed to normal to surface.
% p = [R, H, alpha] 
% when Height is not enough, H = R*tan(alpha) will be replaced.
% 2004. May, 14.
qp = vect2row(qp)';
qz = vect2row(qz)';
if p(2)/p(1) >= tan(deg2rad(p(3)))
    p(2) = p(1)*tan(deg2rad(p(3)));
%    disp('too small alpha value')
%    return
end
    R = p(1);
    H = p(2);
    alpha = p(3);
    
    
    if alpha < 90
        integralstep = H/100;
        z = 0:integralstep:H; %% here, integralstep should be changed depending on qz
                           % as qz is larger, integralstep should be
                           % smaller

        lengthz = length(z);
 %       if is_matrix(qp)
            nz = length(z);
            [nx, ny] = size(qp);
            qp = repmat(qp, [1,1,nz]);
		[nzx, nzy] = size(qz);
		if nx*ny == nzx*nzy
            	qz = repmat(qz, [1,1,nz]);
		else
			qz = repmat(qz, [nx,1, nz]);
		end
            z = repmat(z, [nx*ny, 1, 1]);
            z = reshape(z, [nx, ny, nz]);
%        else
%            [qp, z] = meshgrid(qp, z);
%            [qp, qz, z] = meshgrid(qp, qz, z);
%		qz = ones(size(qp))*qz;
%        end
            Rz = R - z/tan(deg2rad(alpha));
            F = zeros(size(z));
            F = 2*pi*Rz.^2.*besseljc(qp.*Rz).*exp(-j.*qz.*z)*integralstep;
            Formfactor = sum(F, 3);
    elseif alpha == 90
        Rz = R;
        F = 2*pi*Rz.^2*H*besseljc(qp.*Rz).*sinc(H*qz/2);
        Formfactor = sum(F, 3);
    else
        disp('angle should be less than 90 degree')
    end
      
    
