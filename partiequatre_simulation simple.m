clear all
close all



Gamma = 1;
Beta = 5;

I(1) = 1;
S(1) = 9999;
R(1) = 0;
T(1) = 0;

N = I(1) + S(1) + R(1);


i =1;
k = 1;

while I(i) >= 1
    u = rand;
    l(i) = Beta*S(i)*I(i)./N + Gamma*I(i);
    dt(i) = -log(u)./l(i);

    Pinfection = Beta*S(i)*I(i)./(N*l(i));
    Pguerison = Gamma*I(i)./l(i);

    a = rand;
    z = 0;
    if a < Pinfection
        z = z +1; 
        S(i+1) = S(i) -1;
    else
        S(i+1) = S(i);
    end

   if a > Pinfection
       R(i+1) = R(i) +1; 
       z = z-1;   
   else
       R(i+1) = R(i);
   end
    I(i+1) = I(i) + z;
    T(i+1) = T(i) + dt(i) ;

    i= i+1;


    if T(i) < 0.5
        k = k+1;
    end

        

end 
      
%i
%T(i)
%e = max(I)


figure(1)

plot(T,I)
hold on
plot(T,R)
hold on
plot(T,S)

grid on

xlabel 'temps'
ylabel 'populations'
legend("I", "R", "S")

figure(2)

plot(T(1:k), I(1:k), 'o')
hold on
plot(T(1:k), R(1:k), '*')


