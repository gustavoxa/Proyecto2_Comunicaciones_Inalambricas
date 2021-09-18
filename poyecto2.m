%%
% *Escuela Politecnica Nacional*
% Comunicaciones Inalambricas
% Proyecto N° 2
% Stalin Ramirez
% Gustavo Aconda
%% Ejercicio 1 
% Datos
Nsim=1e5; 
Ntx1=2;Nrx1=2;
Ntx2=2;Nrx2=3;
Ntx3=3;Nrx3=2;
Rthdb=8;
Rth=db2pow(Rthdb);
Pt_db=0:1:30;
Pt_ln=10.^(Pt_db/10); tam=length(Pt_db);
thetha=1/sqrt(2);
psalida1=zeros(1,tam);
psalida2=zeros(1,tam);
psalida3=zeros(1,tam);
% Curvas
for i=1:Nsim
hray1=(randn(Ntx1,Nrx1,tam)+1i.*randn(Ntx1,Nrx1,tam))*thetha;  
hray2=(randn(Ntx2,Nrx2,tam)+1i.*randn(Ntx2,Nrx2,tam))*thetha;
hray3=(randn(Ntx3,Nrx3,tam)+1i.*randn(Ntx3,Nrx3,tam))*thetha;
hnor1= sum((abs(hray1)).^2, 2);
hnor2= sum((abs(hray2)).^2, 2);
hnor3= sum((abs(hray3)).^2, 2);
hm1= max(hnor1);   
hm2= max(hnor2);   
hm3= max(hnor3);   
gamma1 = zeros(1,length(Pt_db)); 
gamma2 = zeros(1,length(Pt_db)); 
gamma3 = zeros(1,length(Pt_db)); 
for j=1:tam
    aux1= Pt_ln(j).*(hm1(j));
    aux2= Pt_ln(j).*(hm2(j));
    aux3= Pt_ln(j).*(hm3(j));
    gamma1(j)= aux1;
    gamma2(j)= aux2;
    gamma3(j)= aux3;
end
    aux11= log2((1+ gamma1))<Rth;
    aux12= log2((1+ gamma2))<Rth; 
    aux13= log2((1+ gamma3))<Rth; 
    psalida1=psalida1+aux11;
    psalida2=psalida2+aux12;
    psalida3=psalida3+aux13;
end
% Curva 1
ray1=psalida1/Nsim;
PC1=((((2^Rth)-1).^(Nrx1*Ntx1))./((factorial(Nrx1)).^(Ntx1))).*((1./Pt_ln).^(Nrx1*Ntx1));
% Curva 2
ray2=psalida2/Nsim;
PC2=((((2^Rth)-1).^(Nrx2*Ntx2))./((factorial(Nrx2)).^(Ntx2))).*((1./Pt_ln).^(Nrx2*Ntx2));
% Curva 3
ray3=psalida3/Nsim;
PC3=((((2^Rth)-1).^(Nrx3*Ntx3))./((factorial(Nrx3)).^(Ntx3))).*((1./Pt_ln).^(Nrx3*Ntx3));
figure
semilogy(Pt_db,ray1,'* r');
hold on
semilogy(Pt_db,PC1,'- g')
semilogy(Pt_db,ray2,'* b');
semilogy(Pt_db,PC2,'- g')
semilogy(Pt_db,ray3,'* y');
semilogy(Pt_db,PC3,'- g')
grid
ylim([10^-4.5 10^0])
title("Sistema MIMO")
ylabel("Probabilidad Outage")
xlabel("Pt/N0(db)")
legend('Nt=2,Nr=2','','Nt=2,Nr=3','','Nt=3,Nr=2','Solución Analitica')
%% Ejercicio 2 
% Datos
% Datos
Nsim=1e5; 
Ntx=2;Nrx=2;
Rthdb1=0;
Rthdb2=4;
Rthdb3=8;
Rth1=db2pow(Rthdb1);
Rth2=db2pow(Rthdb2);
Rth3=db2pow(Rthdb3);
Pt_db=-5:3:30;
Pt_ln=10.^(Pt_db/10); tam=length(Pt_db);
thetha=1/sqrt(2);
psalida1=zeros(1,tam);
psalida2=zeros(1,tam);
psalida3=zeros(1,tam);
for i=1:Nsim
hray1=(randn(Ntx,Nrx,tam)+1i.*randn(Ntx,Nrx,tam))*thetha;  
hnor1= sum((abs(hray1)).^2, 2);
hm1= max(hnor1);   
gamma1 = zeros(1,length(Pt_db)); 
for j=1:tam
    aux1= Pt_ln(j).*(hm1(j));
    gamma1(j)= aux1;
end
    aux11= log2((1+ gamma1))<Rth1;
    aux12= log2((1+ gamma1))<Rth2; 
    aux13= log2((1+ gamma1))<Rth3; 
    psalida1=psalida1+aux11;
    psalida2=psalida2+aux12;
    psalida3=psalida3+aux13;
end
% Curva 1
ray1=psalida1/Nsim;
PC1=((((2^Rth1)-1).^(Nrx*Ntx))./((factorial(Nrx)).^(Ntx))).*((1./Pt_ln).^(Nrx*Ntx));
% Curva 2
ray2=psalida2/Nsim;
PC2=((((2^Rth2)-1).^(Nrx*Ntx))./((factorial(Nrx)).^(Ntx))).*((1./Pt_ln).^(Nrx*Ntx));
% Curva 3
ray3=psalida3/Nsim;
PC3=((((2^Rth3)-1).^(Nrx*Ntx))./((factorial(Nrx)).^(Ntx))).*((1./Pt_ln).^(Nrx*Ntx));
figure
semilogy(Pt_db,ray1,'- r');
hold on
semilogy(Pt_db,PC1,'* g')
semilogy(Pt_db,ray2,'- b');
semilogy(Pt_db,PC2,'* g')
semilogy(Pt_db,ray3,'- y');
semilogy(Pt_db,PC3,'* g')
grid
ylim([10^-5 10^0])
title("Sistema MIMO")
ylabel("Probabilidad Outage")
xlabel("Pt/N0(db)")
legend('Rth=0','','Rth=4','','Rth=8','Solución Analitica')
