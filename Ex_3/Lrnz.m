%% Lorenz Function
function dx=Lrnz(t,x,dummy,sig,b,rho)

dx=[sig*(-x(1)+x(2))
    -x(1)*x(3)+rho*x(1)-x(2)
    x(1)*x(2)-b*x(3)];

end
