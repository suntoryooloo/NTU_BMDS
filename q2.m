
%% set parameters
clear;
clc;
close all;

h=1e-5;      
t=0:h:40;   
k_1=1.67;
k_2=10;
k_3=2.5;

%% Create four arrays of calculation results 
N=length(t);
E=ones(1,N);
S=ones(1,N);
ES=ones(1,N);
P=ones(1,N);

%% fourth-order RungeKutta method
E(1)=1;
S(1)=10;
ES(1)=0;
P(1)=0;
for i=2:N
    t_n=t(i-1);
    e_n=E(i-1);
    s_n=S(i-1);
    es_n=ES(i-1);
    p_n=P(i-1);
    
    ke1=k_2*es_n+k_3*es_n-k_1*e_n*s_n;
    ks1=k_2*es_n-k_1*e_n*s_n;
    kes1=k_1*e_n*s_n-k_2*es_n-k_3*es_n;
    kp1=k_3*es_n;
    
    ke2=k_2*(es_n + kes1*h/2)+ k_3*( es_n + kes1*h/2) -k_1*(e_n+ ke1*h/2)*(s_n+ ks1*h/2);
    ks2=k_2*(es_n + kes1*h/2)-k_1*(e_n+ ke1*h/2)*(s_n+ ks1*h/2);
    kes2=k_1*(e_n + ke1*h/2)*(s_n+ ks1*h/2)- k_2*(es_n + kes1*h/2)- k_3*( es_n + kes1*h/2);
    kp2=k_3*( es_n + kes1*h/2);
    
    ke3= k_2*(es_n + kes2*h/2)+ k_3*( es_n + kes2*h/2) -k_1*(e_n+ ke2*h/2)*(s_n+ ks2*h/2);
    ks3= k_2*(es_n + kes2*h/2) -k_1*(e_n+ ke2*h/2)*(s_n+ ks2*h/2);
    kes3=k_1*(e_n+ ke2*h/2)*(s_n+ ks2*h/2)- k_2*(es_n + kes2*h/2)- k_3*( es_n + kes2*h/2);
    kp3= k_3*( es_n + kes2*h/2);    

    ke4= k_2*(es_n + kes3*h)+ k_3*( es_n + kes3*h) -k_1*(e_n+ ke3*h)*(s_n+ ks3*h);
    ks4= k_2*(es_n + kes3*h) -k_1*(e_n+ ke3*h)*(s_n+ ks3*h);
    kes4=k_1*(e_n+ ke3*h)*(s_n+ ks3*h)- k_2*(es_n + kes3*h)- k_3*( es_n + kes3*h);
    kp4= k_3*( es_n + kes3*h);    

    E(i)=e_n+h/6*(ke1+2*ke2+2*ke3+ke4);
    S(i)=s_n+h/6*(ks1+2*ks2+2*ks3+ks4);
    ES(i)=es_n+h/6*(kes1+2*kes2+2*kes3+kes4);  
    P(i)= p_n+h/6*(kp1+2*kp2+2*kp3+kp4);
end

%% plot the profilr
figure();
hold on;
plot(t,E,'r');
plot(t,ES,'g');
plot(t,S,'b');
plot(t,P,'k');
legend('E','ES','S','P');
xlabel('Time(s)');
title('Concentration Profiles');
hold off;

%% 8.3
k=(k_2+k_3)/k_1;
v_max=k_3*E(1);
v=zeros(1,N);
s=linspace(1,300,N);

for i=1:N
    v(1,i)=s(i)*v_max/(k+s(i));
end

figure();
hold on;
plot(s,v,'r');
xlabel('Concentration of the S(\mu M)');
xlabel('V(\mu M/s)');
title('Velocity V as a function of the Concentration of the S ');
hold off;




