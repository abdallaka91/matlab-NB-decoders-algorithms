function dgf = HD_bin_gf(d,p,N,pw)
% Nb = length(d);
% N = Nb/p;
dbin = 0.5*(1-d);
dgf = nan(1,N);
for n = 1 : N
    % bin_symb = dbin((n-1)*p+1:n*p);
    % dgf(n) = bi2de(bin_symb,p);
    dgf(n)  = dbin((n-1)*p+1:n*p)*pw;
end