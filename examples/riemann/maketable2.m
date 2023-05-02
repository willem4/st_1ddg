gg = size(L);
for i = 1:3; 
    Table(1,(i-1)*2+1) = L(i,1);
    Table(1,(i-1)*2+2) = 0;
end
for jj=2:gg(2);
    for i = 1:3; 
        Table(jj,(i-1)*2+1) = L(i,jj);
        Table(jj,(i-1)*2+2) = 0.01*round(100*(log2(L(i,jj-1)./L(i,jj))));
    end
end

%Table(:,1)=L(:,1);
%for jj=1:gg(2)-1;
%    Table(:,jj*2)=L(:,jj+1);
%    Table(:,jj*2+1)=0.01*round(100*(log2(L(:,jj)./L(:,jj+1))));
%end