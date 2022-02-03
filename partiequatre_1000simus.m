clear all
close all



Gamma = 1;
Beta = 5;
E = 0;


for k = 1:100000
    I(1) = 1;
    S(1) = 2;
    R(1) = 0;
    T(1) = 0;

    N = I(1) + S(1) + R(1);


    i =1;


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
        T(i+1) = T(i) + dt(i);

        i= i+1;  
    end 


    if S(i) >= R(i)
        E = E+1;
    end
end

E

