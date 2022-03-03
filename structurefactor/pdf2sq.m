function sq = pdf2sq(r, gr, q, n_p, is2D)
% sq = pdf2sq(r, pdf, q)
% here pdf stands for gr not Gr.

r = r(:);deltar = r(2:end)-r(1:end-1);
dr = min(abs(deltar));
if numel(dr) > 1
    dr = dr(1);
end

gr = gr(:);
q = q(:)';
tq = ones(size(q));
%temp = n_p*(((gr-1)*4*pi.*r.^2)*tq).*sin(r*q+eps)./(r*q+eps);
temp = n_p*(((gr-1)*2*pi.*r)*tq).*sin(r*q)./(r*q);
sq = sum(temp, 1)*dr;
