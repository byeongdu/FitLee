function Iq = powderpeaks(q, pos, amp, p)
%function err = powderpeaks(param,x,y,plot_handle,text_handle)
%err = myTestFunction(p,x,y,plot_handle,text_handle)
%
%parameter = [gw, dsize, mstrain, bscale, exponent];
waveln = 1;
%dsize = parameter(end-3);
%mstrain = parameter(end-2);
%bscale = parameter(end-1);
%exponent = parameter(end);
Iq = zeros(size(q));
for i=1:numel(pos)
    w = peakwidth(pos(i), waveln, p.dsize, p.mstrain); 
    Iq = Iq + pseudovoigt(q, [amp(i), pos(i), p.gw, w]);
end

