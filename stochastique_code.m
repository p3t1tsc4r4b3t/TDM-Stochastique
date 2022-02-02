% Coeffs du pb %


N=1000; %nb personnes
B=1.5; % taux B auqeul chaque personne infectée I infecte une personne saine
gam=0.5; % taux auquel chaque personne infectee rejoindre la pop ayant recupéré
R_0=B/gam;
T=[0 100];
%Conditions initiales%
S_z=89;
I_z=11;
R_z=0;
Y_z=[S_z,I_z,R_z];
param=R_0*S_z/N
%% PARTIE 3 Modèle déterministe %%%
[T_sol,Y_sol]=ode45(@(T_sol,Y_sol) pb(T_sol,Y_sol,gam,B,N),T,Y_z);


figure(1)
plot(T_sol,Y_sol(:,1))
hold on
plot(T_sol,Y_sol(:,2))
hold on
plot(T_sol,Y_sol(:,3))


legend("S","I","R")
xlabel("temps")
ylabel("Nombre de personnes")
title("R_0*S(t=0)/N = 5")




%%% PARTIE 5 Limite diffusive %%%
n = 2001; %nombre de pas de temps
t_i = 0;
t_f = 200;
t = linspace(t_i,t_f,n);
dt = (t_f - t_i)/(n-1); % A vérifier pour le "-1"


%Critère de fin d'épidémie


%Initialisation
i_0=1
s_0=N-i_0
I=[i_0]
S=[s_0]
R=[0]
%Wiener




%Méthode d'euler
for i=[1:n-1]
dwi=normrnd(0,dt);
dwg=normrnd(0,dt);
S(i+1)=S(i)-B/N*I(i)*S(i)*dt+sqrt(B/N*S(i)*I(i))*dwi;
I(i+1)=I(i)+(B/N*I(i)*S(i)-gam*I(i))*dt-sqrt(B/N*S(i)*I(i))*dwi +sqrt(gam*I(i))*dwg;
R(i+1)=R(i)-(S(i+1)-S(i)+I(i+1)-I(i));
end


figure(2)
plot(t,S)
hold on
plot(t,I)
hold on
plot(t,R)
legend("S","I","R")




%Ouverture%
%Reinfection
r=0.001 %paramètre plus injection
i_0=1
s_0=N-i_0
I=[i_0]
S=[s_0]
R=[0]
T=[I(1)+S(1)+R(1)]


for i=[1:n-1]
dwi=normrnd(0,dt);
dwg=normrnd(0,dt);
S(i+1)=S(i)-B/N*I(i)*S(i)*dt+sqrt(B/N*S(i)*I(i))*dwi;
I(i+1)=I(i)+(B/N*I(i)*S(i)-gam*I(i))*dt-sqrt(B/N*S(i)*I(i))*dwi +sqrt(gam*I(i))*dwg;
R(i+1)=R(i)-(S(i+1)-S(i)+I(i+1)-I(i));
R(i+1)=R(i+1)-r*R(i);
S(i+1)=S(i+1)+r*R(i);
T(i+1)=S(i+1)+I(i+1)+R(i+1);
end


figure(3)
plot(t,S)
hold on
plot(t,I)
hold on
plot(t,R)
hold on
plot(t,T)
legend("S","I","R")




%Ouverture%
%Reinfection + vaccination
%paramètre plus injection
v=0.001 %paramètre vaccination


I=[i_0]
S=[s_0]
R=[0]
V=[0]
T=[I(1)+S(1)+R(1)+V(1)]


for i=[1:n-1]
dwi=normrnd(0,dt);
dwg=normrnd(0,dt);
S(i+1)=S(i)-B/N*I(i)*S(i)*dt+sqrt(B/N*S(i)*I(i))*dwi;
I(i+1)=I(i)+(B/N*I(i)*S(i)-gam*I(i))*dt-sqrt(B/N*S(i)*I(i))*dwi +sqrt(gam*I(i))*dwg;
R(i+1)=R(i)-(S(i+1)-S(i)+I(i+1)-I(i));


V(i+1)=V(i) +S(i)*v;
R(i+1)=R(i+1)-r*i*(t(2)-t(1))/t_f*R(i);
S(i+1)=S(i+1)+r*i*(t(2)-t(1))/t_f*R(i) -S(i)*v;
T(i+1)=S(i+1)+I(i+1)+R(i+1)+V(i+1);
end


figure(4)
plot(t,S)
hold on
plot(t,I)
hold on
plot(t,R)
hold on
plot(t,V)
hold on
plot(t,T)
legend("S","I","R","V")




%beta peut etre stocha, dependre du temps (confinement,..)
%rajouter une equation dr dans le dernier modèle
%rajouter une equation dV pour vacinés cf wikipedia








function [f]= pb(T_sol,Y_sol,gam,B,N)
f=zeros(3,1);
f(1)=-B/N*Y_sol(2)*Y_sol(1);
f(2)=B/N*Y_sol(1)*Y_sol(2)-gam*Y_sol(2);
f(3)=gam*Y_sol(2);
end
