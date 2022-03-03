function [Iq, Aq] = saxs_spiralaverage(q, x, y, z, atm, Fq, x1, y1, z1, Am)
if nargin == 3
    atmtype = 1;
    Fq = ones(size(q));
end
if nargin == 5
    Fq = [];
end

if nargin < 7
    Am = zeros(size(q));
    x1 = [];
    y1 = [];
    z1 = [];
end

q = q(:);
atmnumbers = unique(atm);
atmtype = 1:numel(atmnumbers);

if isempty(Fq)
    atmnumbers = unique(atm);atmnumbers = atmnumbers(:)';
    FFq = ones(size(q));
    Fq = repmat(FFq, [1, numel(atmnumbers)]);
    for i=1:numel(atmnumbers)
        Fq(:,i) = atmnumbers(i)*Fq(:,i);
    end
    atmtype = 1:numel(atmnumbers);
end

for i=1:numel(atmtype)
    atm(atm == atmnumbers(i)) = atmtype(i);
end

j = sqrt(-1);
[xq, yq, zq] = spiral_on_sphere(1);
[Iq, Aq] = cal_loop_spiral;
%[Iq, Aq] = cal_pdb_spiral;

    function [Iq, Aq] = cal_loop_spiral
        Iq = zeros(size(q));
        Aq = zeros(size(q));
        for qindex = 1:numel(q)
            qx = xq * q(qindex);
            qy = yq * q(qindex);
            qz = zq * q(qindex);

            %Aatm = zeros(size(qx));
            
            for qxindex = 1:numel(qx)
                if numel(atmtype) == 1
                    FF = Fq(qindex);
                else
                    if size(Fq, 1) >= qindex
                        FF = Fq(qindex, atm)';
                    else
                        FF = Fq(1, atm)';
                    end
                end
                Tatm = FF.*exp(j*(qx(qxindex)*x+ ...
                    qy(qxindex)*y+qz(qxindex)*z));
                if ~isempty(x1)
                    Aaddition = Am(qindex)*exp(j*(qx(qxindex)*x1+...
                        qy(qxindex)*y1+qz(qxindex)*z1));
                    Iq(qindex) = Iq(qindex) + abs(sum(Tatm)+Aaddition).^2;
                    Aq(qindex) = Aq(qindex) + sum(Tatm)+Aaddition;
                else
                    Iq(qindex) = Iq(qindex) + abs(sum(Tatm)).^2;
                    Aq(qindex) = Aq(qindex) + sum(Tatm);
                end
            end
        end
        Iq = Iq/numel(qx);
        Aq = Aq/numel(qx);
    end

    function [Iq, Aq] = cal_pdb_spiral
        Iq = zeros(size(q));
        Aq = zeros(size(q));
        for qindex = 1:numel(q)
            qx = xq * q(qindex);
            qy = yq * q(qindex);
            qz = zq * q(qindex);

            %Aatm = zeros(size(qx));
            
            for qxindex = 1:numel(qx)
                Tatm = fq_pdb(qx(qxindex), qy(qxindex), qz(qxindex),...
                    x, y, z, atm);
                Iq(qindex) = Iq(qindex) + abs(sum(Tatm)).^2;
                Aq(qindex) = Aq(qindex) + Tatm;
            end
        end
        Iq = Iq/numel(qx);
        Aq = Aq/numel(qx);
    end

    function Iq = cal_loop_atm
        Iq = zeros(size(q));
        for qindex = 1:numel(q)
            qx = xq * q(qindex);
            qy = yq * q(qindex);
            qz = zq * q(qindex);

            Aatm = zeros(size(qx));
            for atmindex = 1:numel(x)
                if numel(atmtype) == 1
                    FF = Fq(qindex);
                else
                    FF = Fq(qindex, atm(atmindex));
                end

                Tatm = FF*exp(j*(qx*x(atmindex)+ ...
                    qy*y(atmindex)+qz*z(atmindex)));
                Aatm = Aatm + Tatm;
            end
            Iatm = abs(Aatm+Am).^2; Iq(qindex) = mean(Iatm);
        end
    end
end