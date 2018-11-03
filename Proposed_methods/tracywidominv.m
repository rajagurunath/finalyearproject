function [pdf,cdf] = tracywidominv(x)
i=2;
%
%
% [pdf, cdf]=tracywidom(x, i)
%
% SHIFTED GAMMA APPROXIMATION FOR THE TRACY-WIDOM LAWS
%
% TW ~ Gamma[k,theta]-alpha
%
% [pdf,cdf]=tracywidom(x,2) gives the approximated pdf(x) and CDF(x) of the
%TW2
%
% [pdf,cdf]=tracywidom(x,i) for i=1,2,4 gives TW1, TW2, TW4
%
% see paper M.Chiani, "Distribution of the largest eigenvalue for real
% Wishart and Gaussian random matrices and a simple approximation for the
% Tracy-Widom distribution", submitted 2012, ArXiv
%
k = [46.44604884387787, 79.6594870666346, 0, 146.0206131050228]; % K,
%THETA, ALPHA
theta = [0.18605402228279347, 0.10103655775856243, 0, 0.05954454047933292];
alpha = [9.848007781128567, 9.819607173436484, 0, 11.00161520109004];
cdf = cdfgamma(x+alpha(i), theta(i), k(i)); % Chiani's appx: Tracy-Widom as
%a shifted Gamma
pdf = pdfgamma(x+alpha(i), theta(i), k(i)); % Chiani's appx: Tracy-Widom as
%a shifted Gamma
end
function [pdfgamma]=pdfgamma(x, ta, ka) % Utility: PDF of a Gamma
if(x > 0)
 pdfgamma=1/(gamma(ka)*ta^ka) * x.^(ka - 1) .* exp(-x/ta);
 pdfgamma=inv(pdfgamma);
else
 pdfgamma=0 ;
end
end
function[cdfgamma]=cdfgamma(x, ta, ka) % Utility: CDF of a Gamma
if(x > 0)
 cdfgamma=gammainc(x/ta,ka);
 cdfgamma=inv(cdfgamma);
else
 cdfgamma=0;
end
end