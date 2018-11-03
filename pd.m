function pd()
clear all
clc
close all
N=1000;
H=1.8858922;
T=1.5548356;
SNR=[];
Pr=[];
PM1=[];
PM2=[];
for j=1:33;
    Pr(j)=0.0092*(1.1222864^j);
    S1=0;
    S2=0;
    for i=1:10000
        P=random('normal',0,sqrt(Pr(j)),1,N);
        N1=random('normal',0,sqrt(exp(0.5)),1,N);
        N2=random('normal',0,sqrt(exp(0.5)),1,N);
        N3=random('normal',0,sqrt(exp(0.5)),1,N);
        N4=random('normal',0,sqrt(exp(0.5)),1,N);
        N5=random('normal',0,sqrt(exp(0.5)),1,N);
        N6=random('normal',0,sqrt(exp(0.5)),1,N);
        N7=random('normal',0,sqrt(exp(0.5)),1,N);
        N8=random('normal',0,sqrt(exp(0.5)),1,N);
        N9=random('normal',0,sqrt(exp(0.5)),1,N);
        N10=random('normal',0,sqrt(exp(0.5)),1,N);
        A=[N1+P;N2+P;N3+P;N4+P;N5+P;N6+P;N7+P;N8+P;N9+P;N10+P];
        B=A*(A');
        R=B/N;
        H2=max(eig(R));
        H1=min(eig(R));
        d=H2/H1;
        if d>=T
           w=1;
        elseif d<T
           w=0;
        end
        G=0;
        for k=1:N
             G=G+(P(k)+N1(k))^2;
        end
        G1=G/N;
        if G1>=H
           a1=1;
        elseif G1<H
           a1=0;
        end
        G=0;
        for k=1:N
             G=G+(P(k)+N2(k))^2;
        end
        G2=G/N;
        if G2>=H
           a2=1;
        elseif G2<H
           a2=0;
        end
        G=0;
        for k=1:N
             G=G+(P(k)+N3(k))^2;
        end
        G3=G/N;
        if G3>=H
           a3=1;
        elseif G3<H
           a3=0;
        end
        G=0;
        for k=1:N
             G=G+(P(k)+N4(k))^2;
        end
        G4=G/N;
        if G4>=H
           a4=1;
        elseif G4<H
           a4=0;
        end
        G=0;
        for k=1:N
             G=G+(P(k)+N5(k))^2;
        end
        G5=G/N;
        if G5>=H
           a5=1;
        elseif G5<H
           a5=0;
        end
        G=0;
        for k=1:N
             G=G+(P(k)+N6(k))^2;
        end
        G6=G/N;
        if G6>=H
           a6=1;
        elseif G6<H
           a6=0;
        end
        G=0;
        for k=1:N
             G=G+(P(k)+N7(k))^2;
        end
        G7=G/N;
        if G7>=H
           a7=1;
        elseif G7<H
           a7=0;
        end
        G=0;
        for k=1:N
             G=G+(P(k)+N8(k))^2;
        end
        G8=G/N;
        if G8>=H
           a8=1;
        elseif G8<H
           a8=0;
        end
        G=0;
        for k=1:N
             G=G+(P(k)+N9(k))^2;
        end
        G9=G/N;
        if G9>=H
           a9=1;
        elseif G9<H
           a9=0;
        end
        G=0;
        for k=1:N
             G=G+(P(k)+N10(k))^2;
        end
        G10=G/N;
        if G10>=H
           a10=1;
        elseif G10<H
           a10=0;
        end
        if a1+a2+a3+a4+a5+a6+a7+a8+a9+a10>=1
           x=1;
        elseif a1+a2+a3+a4+a5+a6+a7+a8+a9+a10==0
           x=0;
        end
        S1=S1+w;
        S2=S2+x;
     end
    PM1(j)=1-S1/100000;
    PM2(j)=1-S2/100000;
    SNR(j)=10*log10(Pr(j)/(exp(0.5)));
end
semilogy(SNR,PM1,'x-k',SNR,PM2,'o-k');
legend('Simulation','Analytical');
 axis([-22 -6 0 1]);
xlabel('Threshold');
ylabel('Detection Probability');
grid on
end