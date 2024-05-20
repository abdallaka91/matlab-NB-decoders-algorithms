function y = function_f(r)
N = length(r);
m=N/2;
ii = 1;
i1= ii:ii+m-1;
i2=ii+m:ii+2*m-1;
v1 = r(i1);
v2 = r(i2);
y = (1-2*(v1<0)).*(1-2*(v2<0)).*min(abs(v1), abs(v2));
end

