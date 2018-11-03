function User_6()
clear all;
clc;

N=100;                  
Pr=20;                  
Pf= 0.01:0.01:0.2;
UserNum = 6;
SimulationNum = 10000;
for ii = 1: length(Pf)
    TH=chi2inv(1-Pf(ii),N);
    S=0;                                                                  
    for i=1:SimulationNum
        P=sqrt(Pr/N)*ones(1,N);
        
        NN = zeros(UserNum,N);
        GG = zeros(UserNum,1);
        DD = zeros(UserNum,1);
        for j = 1 : UserNum
            NN(j,:) = randn(1,N);
        end
                              
        for j = 1 : UserNum
            GG(j)=0;
            for n=1:N
             GG(j)=GG(j)+(P(n)+NN(j,n))^2;
            end
            if GG(j)>=TH
               DD(j) = 1;
            elseif GG(j)<TH
               DD(j) = 0;
            end
        end
        
        sum = 0;
        for j = 1 : UserNum
            sum = sum + DD(j);
        end
        if sum >=1
            D=1;
        elseif sum ==0
            D=0;
        end
        S=S+D;
    end  
    Pd(ii)=S/SimulationNum;
end



pn=0.1;
gama=0.16;
for i=1:1:20
Pf(i)=0.001+(i-1)*0.005 
thr=N*pn^2+sqrt(2*N)*pn^2*erfcinv(Pf(i)); 
Pds(i)=erfc((thr-(N+gama*N)*pn^2)/(sqrt(2*(N+2*gama*N))*pn^2)) 
end
figure;
semilogx(Pf,Pds,'-+');
hold on
semilogx(Pf,Pd,'linewidth',2);
xlabel('Pf')
ylabel('Pd')
legend('Sphere detector','John detector')
Title ('Primary User =6')
grid on;

end