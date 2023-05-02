%global x t

%x=0.3;
t=32.2;%32.2;

fid = fopen('init.1280');
Ma=fscanf(fid,'%u',[2 1]);
M1a = Ma(1); %fscanf(fid,'%u',[1 1]);
M2a = Ma(2); %fscanf(fid,'%u',[1 1]);
Nn = fscanf(fid,'%u',[1 1]);
lambda = (sqrt(3)-1)/(2*sqrt(3));
fid2 = fopen('exact.1280','w');
fprintf(fid2,'%u\t%u\n',M1a,M2a);
fprintf(fid2,'%u\n',Nn);
for ii = 1:2:Nn
   aa = fscanf(fid,'%g %g',[2 2]); % It has two rows now.
   a = aa(1,1);
   b = aa(1,2);
   xi1 = a+lambda*(b-a);
   xi2 = b-lambda*(b-a);
   fxi1=1-myfun(xi1,t,1);
   fxi2=1-myfun(xi2,t,1);
   fmean = (fxi1+fxi2)/2;
   fslope = (fxi2-fxi1)/(xi2-xi1)*(b-a)/2;
   aa(2,1) = fmean-fslope;
   aa(2,2) = fmean+fslope;
   for jj = 1:2
        fprintf(fid2,'%.15g\t%.15g\n',aa(1,jj),aa(2,jj)); % It has two rows now.   
   end
end
for ii = Nn+1:2:2*Nn
   aa = fscanf(fid,'%g %g',[2 2]); % It has two rows now.
   a = aa(1,1);
   b = aa(1,2);
   xi1 = a+lambda*(b-a);
   xi2 = b-lambda*(b-a);
   fxi1=myfun(xi1,t,2);
   fxi2=myfun(xi2,t,2);
   fmean = (fxi1+fxi2)/2;
   fslope = (fxi2-fxi1)/(xi2-xi1)*(b-a)/2;
   aa(2,1) = fmean-fslope;
   aa(2,2) = fmean+fslope;
   for jj = 1:2
        fprintf(fid2,'%.15g\t%.15g\n',aa(1,jj),aa(2,jj)); % It has two rows now.   
   end
end
for ii = 2*Nn+1:2:(M1a+M2a)*Nn
   aa = fscanf(fid,'%g %g',[2 2]); % It has two rows now.
   a = aa(1,1);
   b = aa(1,2);
   xi1 = a+lambda*(b-a);
   xi2 = b-lambda*(b-a);
   fxi1=myfun(xi1,t,1);
   fxi2=myfun(xi2,t,1);
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