function F = myfun(x,t,y)
% constants.
g=9.81;
%values far away from the bump.
hb0=-0.2;
H=1;    
Q0=0.4754;
R=Q0+(Q0/(H-hb0))^3;
%tolerance.
tol = 1e-14;
xneeded = x
xp = 10;
delta = 0.05;
hc=0.2;
dudhb=0;
err = 1; 
while abs(err) > tol;
    %prediction of x at time 0.
    x = xneeded-dudhb*t;
    if (abs(x - xp)<sqrt(hc/delta))
        hb = -0.9*hc -0.1*delta*(x-xp)^2;
    else    
        hb = -hc;
    end
    p = [1,0,H-hb,-R]; %=0;
    a = roots(p);
    u = a(3);   
    %du/dhb at time zero. 
    dudhb = 3*u^3/(3*u^2+H-hb);
    %corresponding error.
    err = xneeded-(x+dudhb*t);
end;    
if y == 1;
    %returns height of topography
    F = hb; 
end;
if y == 2;
    %returns velocity.
    p = [1,0,H-hb,-R];
    a = roots(p);
    u = a(3);
    F = R-u^3;
end;    