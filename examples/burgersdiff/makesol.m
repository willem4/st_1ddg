for i = 1:size(x,2)
    f(i) = exact(x(i),10);
end
plot(x,f)
