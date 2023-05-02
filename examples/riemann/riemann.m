function y = riemann(ul,ur,ns)

global s;
g=9.81;                     %zwaartekrachtsconstante

hL=ul(1);
hR=ur(1);
uL=ul(2)/ul(1);              %ul=mL=uL*hL
uR=ur(2)/ur(1);              %idem

h=hL;      %iteratief h* bepalen ... 
h_old=h+1;      %willekeurige h_old om loop in gang te zetten 

while abs(h-h_old)>1.e-15
    if h>hL 
        fLhhL = sqrt(g*(hL-h)^2/(hL*h));
        diff_fLhhL = -1/2/(g*(hL-h)^2/hL/h)^(1/2)*g/hL/h^2*(hL^2-h^2);
    else
        fLhhL = 2*sqrt(g*h)-2*sqrt(g*hL);
        diff_fLhhL = sqrt(g/h);
    end
    if h>hR 
        fRhhR = sqrt(g*(hR-h)^2/(hR*h));
        diff_fRhhR = -1/2/(g*(hR-h)^2/hR/h)^(1/2)*g/hR/h^2*(hR^2-h^2);
    else
        fRhhR = 2*sqrt(g*h)-2*sqrt(g*hR);
        diff_fRhhR = sqrt(g/h);
    end
    Fh=fLhhL+fRhhR+uR-uL;
    diff_Fh=diff_fLhhL+diff_fRhhR;
    h_old = h
    h=h-Fh/diff_Fh;
end
if h<=0 
    return;
end

u=(uR+uL)/2+(fRhhR-fLhhL)/2;
%fLhhL+fRhhR+uR-uL;

for dt = 10:10:120; % 0.01:0.01:0.05
    for i = 1:1001                           %for riemann plots
        x(i)=(i+1-450); % /1001;
        
        % dt=1;
        % i = 1;                                    %for calculation of flux
        % x(1)=0;
        
        
        %--> S-S, R-S, S-R, R-R verschillende types:  Bepaal flux.
        if (h>hL&h>hR)   %SS
            if (uL-sqrt(0.5*g*(hL+h)*h/hL))*dt>x(i)
                hexact(i)=hL;
                uexact(i)=uL;
            end    
            if (uR+sqrt(0.5*g*(hR+h)*h/hR))*dt<x(i)
                hexact(i)=hR;
                uexact(i)=uR;
            end    
            if (uR+sqrt(0.5*g*(hR+h)*h/hR))*dt>=x(i) & (uL-sqrt(0.5*g*(hL+h)*h/hL))*dt<=x(i)
                hexact(i)=h;
                uexact(i)=u;
            end
        end    
        
        if (h<=hL) & (h>hR)     %RS
            if (uL-sqrt(g*hL))*dt>x(i)
                hexact(i)=hL;
                uexact(i)=uL;
            end    
            if (uR+sqrt(0.5*g*(hR+h)*h/hR))*dt<x(i)
                hexact(i)=hR;
                uexact(i)=uR;
            end    
            if (uR+sqrt(0.5*g*(hR+h)*h/hR))*dt>=x(i) & (uL-sqrt(g*hL))*dt<=x(i)
                if (u-sqrt(g*h))*dt<x(i)
                    hexact(i)=h;
                    uexact(i)=u;
                else
                    h2=(uL+2*sqrt(g*hL)-x(i)/dt)^2/(g*9);
                    u2=uL+2*sqrt(g*hL)-2*sqrt(g*h2);
                    hexact(i)=h2;
                    uexact(i)=u2;        
                end
            end
        end    
        if (h>hL) & (h<=hR)       %SR%
            
            if (uL-sqrt(0.5*g*(hL+h)*h/hL))*dt>x(i)
                hexact(i)=hL;
                uexact(i)=uL;
            end    
            if (uR+sqrt(g*hR))*dt<x(i)
                hexact(i)=hR;
                uexact(i)=uR;
            end    
            if (uR+sqrt(g*hR))*dt>=x(i) & (uL-sqrt(0.5*g*(hL+h)*h/hL))*dt<=x(i)
                if (u+sqrt(g*h))*dt>x(i)
                    hexact(i)=h;
                    uexact(i)=u;
                else        
                    h2=(2*sqrt(g*hR)-uR+x(i)/dt)^2/(g*9);
                    u2=uR+2*sqrt(g*h2)-2*sqrt(g*hR);
                    hexact(i)=h2;
                    uexact(i)=u2;
                end    
            end
        end    
        if (h<=hL) & (h<=hR)      %RR
            
            if (uL-sqrt(g*hL))*dt>x(i)
                hexact(i)=hL;
                uexact(i)=uL;
            end    
            if (uR+sqrt(g*hR))*dt<x(i)
                hexact(i)=hR;
                uexact(i)=uR;
            end    
            if (uR+sqrt(g*hR))*dt>=x(i) & (uL-sqrt(g*hL))*dt<=x(i)
                if (u-sqrt(g*h))*dt>=x(i)
                    h2=(uL+2*sqrt(g*hL)-x(i)/dt)^2/(g*9);
                    u2=uL+2*sqrt(g*hL)-2*sqrt(g*h2);
                    hexact(i)=h2;
                    uexact(i)=u2;        
                end  
                if (u+sqrt(g*h))*dt<=x(i)    
                    h2=(2*sqrt(g*hR)-uR+x(i)/dt)^2/(g*9);
                    u2=uR+2*sqrt(g*h2)-2*sqrt(g*hR);
                    hexact(i)=h2;
                    uexact(i)=u2;
                end   
                if (u-sqrt(g*h))*dt<x(i) & (u+sqrt(g*h))*dt>x(i)
                    hexact(i)=h;
                    uexact(i)=u;
                end              
            end
        end    
        
%        y(1)=hexact(1)*uexact(1);                          %calculation of flux
%        y(2)=uexact(1)^2*hexact(1)+0.5*g*hexact(1)^2;
        
    end
    subplot(2,1,1);
    plot(x,uexact.*hexact,'-');
    hold on;
    ylabel('u(x,t)*h(x,t)','fontsize',18);
    subplot(2,1,2);                           %for riemann plots
    plot(x,hexact,'-');
    hold on;
    xlabel('x','fontsize',18);
    ylabel('h(x,t)','fontsize',18);
end

