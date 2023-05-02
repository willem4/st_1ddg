%clear;
fid = fopen('init.160');
Ma=fscanf(fid,'%u',[2 1]);
M1a = Ma(1); %fscanf(fid,'%u',[1 1]);
M2a = Ma(2); %fscanf(fid,'%u',[1 1]);
Nn = fscanf(fid,'%u',[1 1]);
a = fscanf(fid,'%g %g',[2 inf]); % It has two rows now.
fclose(fid);
ext = [160];
for kk = 1:length(ext)
    fid = fopen(['out.',int2str(ext(kk))]);
    Mb=fscanf(fid,'%u',[2 1]);
    M1b = Mb(1); %fscanf(fid,'%u',[1 1]);
    M2b = Mb(2); %fscanf(fid,'%u',[1 1]);
    %   M1b = fscanf(fid,'%u',[1 1]);
    %   M2b = fscanf(fid,'%u',[1 1]);
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
    figure(kk)
    for jj = 1:M1a;
        subplot (M1a+M2a,1,jj)
        plot(a(1,(jj-1)*Nn+1:jj*Nn),a(2,(jj-1)*Nn+1:jj*Nn),'g-');
        hold on;
        %subplot (M1,1,jj)
        plot(b(1,(jj-1)*Nn+1:jj*Nn),b(2,(jj-1)*Nn+1:jj*Nn));
        hold off;
    end
    for jj = 1:M2b;
        subplot (M1a+M2a,1,M1a+jj)
        plot(a(1,(jj+M1a-1)*Nn+1:(jj+M1a)*Nn),a(2,(jj+M1a-1)*Nn+1:(jj+M1a)*Nn),'g-');
        hold on;
        %subplot (M1,1,jj)
        plot(b(1,(jj+M1a-1)*Nn+1:(jj+M1a)*Nn),b(2,(jj+M1a-1)*Nn+1:(jj+M1a)*Nn));
        hold off;
    end    
end
    
    