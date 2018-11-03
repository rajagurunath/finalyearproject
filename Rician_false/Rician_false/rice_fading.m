function r = rice_fading(Kdb, N, Mi)
K = 10^(Kdb/10);
const = 1/(2*(K+1));
x = randn(1,N); 
y = randn(1,N);
r = sqrt(const*((x + sqrt(2*K)).^2 + y.^2));
rt = zeros(1,Mi*length(r));
ki = 1;
for i=1:length(r)
    rt(ki:i*Mi) = r(i); ki = ki+Mi;
end
    r = rt;
    