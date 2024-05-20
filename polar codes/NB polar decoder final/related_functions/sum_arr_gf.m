function y = sum_arr_gf(x)
y = gf(0, x.m);
for i = 1 : length(x)
    y = y + x(i);
end
