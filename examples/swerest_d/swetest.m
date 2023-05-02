%clear;
fid = fopen('init.640');
M1a = fscanf(fid,'%u',[1 1]);
M2a = fscanf(fid,'%u',[1 1]);
Nn = fscanf(fid,'%u',[1 1]);
a = fscanf(fid,'%g %g',[2 inf]); % It has two rows now.
fclose(fid);
ext = [10,20];
for kk = 1:length(ext)
    fid = fopen(['out.',int2str(ext(kk))]);
    M1b = fscanf(fid,'%u',[1 1]);
    M2b = fscanf(fid,'%u',[1 1]);
    %M1b = 1;
    %M2b = 0;
    Nnb = fscanf(fid,'%u',[1 1]);
    clear b;
    b = fscanf(fid,'%g %g',[2 inf]); % It has two rows now.
        %a(2,:)=zeros(1,80);
    % B needs to be extended to the same length as a.
   % b(2,1:Nnb) = 2*ones(1,Nnb)-3*sqrt(b(2,1:Nnb));
    gg = Nn/Nnb;
    b_new = a;
    size(b)
    M1a
    for ii = 1:2:Nnb*(M1b+M2b);
        b_new(2, (ii-1)*gg+1) = b(2,ii);
        for jj = 1:gg-1
            b_new(2, (ii-1)*gg+(jj)*2) = (b(2,ii)*(gg-jj)+b(2,ii+1)*jj)/(gg);
            b_new(2, (ii-1)*gg+(jj)*2+1) = (b(2,ii)*(gg-jj)+b(2,ii+1)*jj)/(gg);
        end
        b_new(2, (ii-1)*gg+2*gg) = b(2,ii+1);
    end
    b=b_new;
    hold on;
    plot(b(1,:),b(2,:),'g-')
    hold off;
    e(1,1:Nn)=a(1,1:Nn);
% for speed 
    e(2,1:Nn)=zeros(size(a(2,Nn+1:2*Nn)));
% for height
   %e(2,1:Nn)=10*ones(size(a(2,2*Nn+1:3*Nn)));
    f(1,1:Nn)=b(1,1:Nn);
% for speed
    f(2,1:Nn)=b(2,Nn+1:2*Nn);
% for height 
  % f(2,1:Nn)=b(2,1:Nn)+b(2,2*Nn+1:3*Nn);
    c(1,:) = e(1,:);
    d(1,:) = e(1,:);
    c(2,:) = abs(f(2,:)-e(2,:));%-a(2,:);
    d(2,:) = c(2,:).^2;
    nrminf = norm(c(2,:),inf);
    nrm1 = 0;
    nrm2 = 0;
    for ii=1:2:length(c),
        nrm1 = nrm1 + (c(2,ii)+c(2,ii+1))/2*(c(1,ii+1)-c(1,ii));
        nrm2 = nrm2 + (d(2,ii)+d(2,ii+1))/2*(d(1,ii+1)-d(1,ii));
    end
    nrm2 = sqrt(nrm2);
    L(1,kk) = nrm1;
    L(2,kk) = nrm2;
    L(3,kk) = nrminf;
  %  sprintf ('L^inf norm = %g', nrminf)
  %  sprintf ('L^1 norm = %g',nrm1)
  %  sprintf ('L^2 norm = %g',nrm2)
    %figure(2)
    for jj = 1:M1a;
        subplot (M1a,1,jj)
        plot(a(1,(jj-1)*Nn+1:jj*Nn),a(2,(jj-1)*Nn+1:jj*Nn),'g-');
        hold on;
        %subplot (M1,1,jj)
        plot(b(1,(jj-1)*Nn+1:jj*Nn),b(2,(jj-1)*Nn+1:jj*Nn));
        hold off;
    end
end
for kk = 2:length(ext)
    for ll = 1:3
        order(ll,kk) = L(ll,kk-1)/L(ll,kk);
    end
end
    
    