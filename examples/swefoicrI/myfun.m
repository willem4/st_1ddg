function F = myfun(x,y)
%global x
g=9.81;
h0=0.4;
if y == 1;
    %subcritical;
    u0 = 0.3962;  % Fr =0.2
end
if y == 2;
    %supercritical;
    u0 = 3.7637;  % Fr = 1.9
end    
%check - Fr = u0/sqrt(g*h0);
hb0=-0.2;
Q=h0*u0;
B=1/2*u0^2+g*(h0+hb0);
xp = 10;
delta = 0.05;
hc = 0.2;
if (abs(x - xp)<sqrt(hc/delta) )
    hb = -delta*(x-xp)^2;
else
    hb = -hc;
end
p = [g,(g*hb-B),0,0.5*Q^2]; %=0;
a = roots(p);
F = a(y);
