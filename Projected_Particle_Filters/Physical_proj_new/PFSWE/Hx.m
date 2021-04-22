function y=Hx(x,inth,minidx,maxidx)
%multiply matix by vector
disp(minidx)
disp(inth)
disp(maxidx)
disp(size(x))
y=x(minidx:inth:maxidx,:);
end