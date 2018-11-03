close all
clc
clear all
M=12;
P=2;
L=8;
Ns=100000;
%% Pd calculation for Optimal Detector with pfa= 0.1
muval=(sqrt(Ns-1)+sqrt(M*L))^2;
vval=(sqrt(Ns-1)+sqrt(M*L))*(((1/sqrt(Ns-1))+(1/sqrt(M*L)))^(1/3));
lamtmp=(sqrt(Ns)-sqrt(M*L))^2;
index=1;
for sigval=-20:1:0
 sigval1=10^(sigval/20);
 pfa=0.1;
 rowval=5e1;
 rowval1=4e1;
 outx=tracywidominv(1-pfa);
 lamda1=((outx*vval)+muval)/lamtmp;
 fin2=(lamda1*Ns)+(Ns*(rowval*M*L-rowval1)/(sigval1^2));
 fin3=(fin2-muval);
 fin=fin3/(vval*1e9);
 fval=(fin);
 finalpd=1-fval;
 final_res(index)=finalpd;
 final_snr(index)=sigval;
 index=index+1;
end
final_res=final_res/max(final_res);
hold on figure,plot(final_snr,final_res,'s:r');
%% Pd calculation for Optimal Detector with pfa= 0.05
M=20;
muval=(sqrt(Ns-1)+sqrt(M*L))^2;
vval=(sqrt(Ns-1)+sqrt(M*L))*(((1/sqrt(Ns-1))+(1/sqrt(M*L)))^(1/3));
lamtmp=(sqrt(Ns)-sqrt(M*L))^2;
index=1;
for sigval=-20:1:0
 sigval1=10^(sigval/20);
 pfa=0.05;
 rowval=5e1;
 rowval1=4e1;
 outx=tracywidominv(1-pfa);
59
 lamda1=((outx*vval)+muval)/lamtmp;
 fin2=(lamda1*Ns)+(Ns*(rowval*M*L-rowval1)/(sigval1^2));
 fin3=(fin2-muval);
 fin=fin3/(vval*1e9);
 fval=(fin);
 finalpd=1-fval;
 final_res(index)=finalpd;
 final_snr(index)=sigval;
 index=index+1;
end
final_res=final_res/max(final_res);
hold on figure,plot(final_snr,final_res,'s:k');
%% Pd calculation for Energy Detector with pfa=0.1
L =8;
snr_dB=-20:1:0;
snr= 10.^(snr_dB./10);
for i=1:length(snr_dB)
 Detect=0;
 Pf=0.1;
 for kk=1:100000 % Number of Monte Carlo Simulations

 %-----AWGN noise with mean 0 and variance 1-----%
 Noise = randn(1,L);
 %-----Real valued Gaussian Primary User Signal------%
 Signal = sqrt(snr(i)).*randn(1,L);
 Recv_Sig = Signal + Noise; % Received signal at SU
 Energy = abs(Recv_Sig).^2; % Energy of received signal over N samples

 %-----Computation of Test statistic for energy detection-----%
 Test_Statistic =(1/L).*sum(Energy);

 %-----Theoretical value of Threshold-----%
 Threshold = (qfuncinv(Pf)./sqrt(L))+ 1;

 if(Test_Statistic >= Threshold) % Check whether the received energy
% is greater than threshold, if so,(Probability of detection) counter by 1
 Detect = Detect+1;
 end
end
 Pd(i) = Detect/kk;
 Pm(i)=1-Pd(i);
 Pd_the(i) = qfunc(((Threshold - (snr(i) +...
1)).*sqrt(L))./(sqrt(2).*(snr(i) + 1)));
 Pm_the(i)=1-Pd_the(i);

end
plot(snr_dB,Pd,'-.b*');
hold on
%% Pd calculation for Energy Detector with pfa=0.05
L =8;
60
snr_dB=-20:1:0;
snr= 10.^(snr_dB./10);
for i=1:length(snr_dB)
 Detect=0;
 Pf=0.05;
 for kk=1:100000 % Number of Monte Carlo Simulations

 %-----AWGN noise with mean 0 and variance 1-----%
 Noise = randn(1,L);
 %-----Real valued Gaussian Primary User Signal------%
 Signal = sqrt(snr(i)).*randn(1,L);
 Recv_Sig = Signal + Noise; % Received signal at SU
 Energy = abs(Recv_Sig).^2; % Energy of received signal over N samples

 %-----Computation of Test statistic for energy detection-----%
 Test_Statistic =(1/L).*sum(Energy);

 %-----Theoretical value of Threshold-----%
 Threshold = (qfuncinv(Pf)./sqrt(L))+ 1;

 if(Test_Statistic >= Threshold) % Check whether the received energy
% is greater than threshold, if so,(Probability of detection) counter by 1
 Detect = Detect+1;
 end
end
 Pd(i) = Detect/kk;
 Pm(i)=1-Pd(i)
 Pd_the(i) = qfunc(((Threshold - (snr(i) +...
1)).*sqrt(L))./(sqrt(2).*(snr(i) + 1)));
 Pm_the(i)=1-Pd_the(i)

end
plot(snr_dB,Pd,'-.r*');
grid on;
title('ESTIMATION ON DETECTION')
xlabel('Signal To Noise Ratio (dB)');
ylabel('Probability Of Detection');
%axis([0,35,0,1]);
legend('GLRT-SFA','RGLD-SFA','JOHN DETECTOR','SPHERE Detector');