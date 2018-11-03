function User_1()
clear all
PF=0.01:0.01:1;
N=100;
TH=chi2inv(1-PF,N);

Pd=[];
for i=1:100
     f1=@(x)(marcumq(sqrt(2*x),sqrt(TH(i)),N/2)).*0.1.*exp(-0.1*x);
     Pd1(i)=quad(f1,0,10000);
     E=log(10);
     D=1;
     f2=@(x)((marcumq(sqrt(2*x),sqrt(TH(i)),N/2)).*(1./(x.*sqrt(2*pi*D))).*exp(-((log(x)-E).^2)./(2*D)));
     Pd2(i)=quad(f2,0.000001,10000);
end
figure;
semilogx(PF,Pd1,'--b',PF,Pd2,'-k','linewidth',2);
% axis([10^-2 10^0 0.8 1]);
xlabel('False Alarm Probability')
ylabel('Detection probabilities')
legend(' Sphere (GLRT)',' John Detector(LMPIT)',4);
grid on;

end