fid = fopen('init.640');
Ma=fscanf(fid,'%u',[2 1]);
M1a = 1;%Ma(1); %fscanf(fid,'%u',[1 1]);
M2a = 0;%Ma(2); %fscanf(fid,'%u',[1 1]);
Nn = fscanf(fid,'%u',[1 1]);
aa = fscanf(fid,'%g %g',[2 inf]); % It has two rows now.
fclose(fid);
lambda = (sqrt(3)-1)/(2*sqrt(3));
for ii = 1:2:Nn
   a = aa(1,ii);
   b = aa(1,ii+1);
   xi1 = a+lambda*(b-a);
   xi2 = b-lambda*(b-a);
   %fxi1 = sin(2*pi*xi1);
   %fxi2 = sin(2*pi*xi2);
   v0 = 0.5;  % Make a starting guess at the solution
   options = optimset('Display','off');  % Turn off Display
   x=xi1;
   fxi1=myfun(x);
   x=xi2;
   fxi2=myfun(x);
   fmean = (fxi1+fxi2)/2;
   fslope = (fxi2-fxi1)/(xi2-xi1)*(b-a)/2;
   aa(2,ii) = fmean-fslope;
   aa(2,ii+1) = fmean+fslope;
end
hold on;
plot (aa(1,1:Nn),aa(2,1:Nn))
fid = fopen('exactv.640','w');
fprintf(fid,'%u\t%u\n',M1a,M2a);
fprintf(fid,'%u\n',Nn);
for ii = 1:Nn
    fprintf(fid,'%.15g\t%.15g\n',aa(1,ii),aa(2,ii)); % It has two rows now.
end
fclose(fid);
clear;
hold off;