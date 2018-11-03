clc;clear all;
SNR=0:10;
M =2;
k =10;
selection =  2;
if selection == 1
    output1='erg';
else
    output1='out';
end
output=zeros(1,length(SNR));
for i=1:length(SNR)
    output(i)=capacity_rician(SNR(i),M,k,output1)
end
plot(SNR,output),grid
xlabel('SNR (dB)')
ylabel(' Capacity (bits/Hz)')
title('Rician Throughput')
