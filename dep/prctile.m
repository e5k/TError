% TError package
% This functions returns the percentile p of the vector X
function yi = prctile(X,p)
x=X(:);
if length(x)~=length(X)
    error('please pass a vector only');
end
n = length(x);
x = sort(x);
Y = 100*(.5 :1:n-.5)/n;
x=[min(x); x; max(x)];
Y = [0 Y 100];
yi = interp1(Y,x,p);

end