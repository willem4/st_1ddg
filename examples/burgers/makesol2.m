global x t

%x=0.3;
t=0.1;

fid = fopen('init.2560');
Ma=fscanf(fid,'%u',[2 1]);
M1a = Ma(1); %fscanf(fid,'%u',[1 1]);
M2a = Ma(2); %fscanf(fid,'%u',[1 1]);
Nn = fscanf(fid,'%u',[1 1]);
lambda = (sqrt(3)-1)/(2*sqrt(3));
fid2 = fopen('exact.2560','w');
fprintf(fid2,'%u\t%u\n',M1a,M2a);
fprintf(fid2,'%u\n',Nn);
v0 = 0.5;  % Make a starting guess at the solution
options = optimset('Display','off');  % Turn off Display
for ii = 1:2:Nn
   aa = fscanf(fid,'%g %g',[2 2]); % It has two rows now.
   a = aa(1,1);
   b = aa(1,2);
   xi1 = a+lambda*(b-a);
   xi2 = b-lambda*(b-a);
   %fxi1 = sin(2*pi*xi1);
   %fxi2 = sin(2*pi*xi2);
   x=xi1;
   [v,Fval,exitflag] = fsolve(@myfun,v0,options);
   fxi1=v;
   x=xi2;
   [v,Fval,exitflag] = fsolve(@myfun,v0,options);
   fxi2=v;
   fmean = (fxi1+fxi2)/2;
   fslope = (fxi2-fxi1)/(xi2-xi1)*(b-a)/2;
   aa(2,1) = fmean-fslope;
   aa(2,2) = fmean+fslope;
   for jj = 1:2
        fprintf(fid2,'%.15g\t%.15g\n',aa(1,jj),aa(2,jj)); % It has two rows now.   
   end
end
fclose(fid);
fclose(fid2);
clear;