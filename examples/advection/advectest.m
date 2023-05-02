%clear;
fid = fopen('init.640');
M1a = fscanf(fid,'%u',[1 1]);
M2a = fscanf(fid,'%u',[1 1]);
Nn = fscanf(fid,'%u',[1 1]);
a = fscanf(fid,'%g %g',[2 inf]); % It has two rows now.
fclose(fid);
ext = [10,20,40,80,160];
for kk = 1:length(ext)
    fid = fopen(['out.',int2str(ext(kk))]);
    M1b = fscanf(fid,'%u',[1 1]);
    M2b = fscanf(fid,'%u',[1 1]);
    Nnb = fscanf(fid,'%u',[1 1]);
    clear b;
    b = fscanf(fid,'%g %g',[2 inf]); % It has two rows now.
    %a(2,:)=zeros(1,80);
    % B needs to be extended to the same length as a.
    gg = Nn/Nnb;
    b_new = a;
    %figure(1)
    %plot(b(1,:),b(2,:),'g-');
    for ii = 1:2:Nnb
        b_new(2, (ii-1)*gg+1) = b(2,ii);
        for jj = 1:gg-1
            b_new(2, (ii-1)*gg+(jj)*2) = (b(2,ii)*(gg-jj)+b(2,ii+1)*jj)/(gg);
            b_new(2, (ii-1)*gg+(jj)*2+1) = (b(2,ii)*(gg-jj)+b(2,ii+1)*jj)/(gg);
        end
        b_new(2, (ii-1)*gg+2*gg) = b(2,ii+1);
    end
    b = b_new;
    c(1,:) = a(1,:);
    d(1,:) = a(1,:);
    c(2,:) = abs(b(2,:)-a(2,:));%-a(2,:);
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
    sprintf ('L^inf norm = %g', nrminf)
    sprintf ('L^1 norm = %g',nrm1)
    sprintf ('L^2 norm = %g',nrm2)
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
    
    