Kdb=7; 
N=100000; 
Mi=1;
r=rice_fading(Kdb,N,Mi);
RdB = 20*log10(r);
Rt = [min(RdB):max(RdB)];
for m = 1:length(Rt)
    fade = find(RdB < Rt(m));
    Nm = length(fade); 
    AF(m) = Nm/N;
end
semilogy(Rt,AF,'k-o'); grid;