function gr = sq2pdf(q, sq, r, n_p)
% function y = sq2pdf(q, sq, r, n_p)
% S(q) to g(r), where n_p is the number density of particles..
% this pair distribution calculation is from the atomic PDF..
% G(r) = 4 pi r n_p [g(r)-1];
Gr = sq2Gr(q,sq,r);
gr = Gr./(4*pi*r*n_p)+1;
