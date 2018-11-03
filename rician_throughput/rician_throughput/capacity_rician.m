function outage=capacity_rician(SNR,M,k,output)
SNR=10^(0.1*SNR);
   for K=1:10000
       T=randn(M,M)+j*randn(M,M);
       T=0.707*T;
%Rician component of 0.8 (choose any value upto 1)       
       H=[1 1;1 1];
       T=sqrt(k/(1+k))*H+sqrt(1/(1+k))*T;
       I=eye(M);
       a=det(I+(SNR/M)*T*T');
       y(K)=log2(a);    
	end
    [n1 x1]=hist(y,40);
    n1_N=n1/max(K);
    a=cumsum(n1_N);
    b=abs(x1);
if output == 'erg'
    outage=interp1q(a,b',0.5);   %ergodic capacity
elseif output == 'out'
    outage=interp1q(a,b',0.1);   %outage capacity 
end 

    
%plot CDF (uncomment the following line if CDF is required)
    plot(abs(x1),a,'b');

