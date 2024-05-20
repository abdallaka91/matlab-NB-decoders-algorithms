function y = function_g(r, u)
N = length(r);
m=N/2;
ii = 1;
i1= ii:ii+m-1;
i2=ii+m:ii+2*m-1;
v1 = r(i1);
v2 = r(i2);
y =  v2+(1-2*u).*v1;
end

