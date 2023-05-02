function F = exact(x,t)
L = 4; 
nu = 1;
alpha = 4;
F = -2*(alpha*nu/L)*(sinh(alpha*x/L)/(cosh(alpha*x/L)+exp(-(alpha/L)^2*nu*t)));