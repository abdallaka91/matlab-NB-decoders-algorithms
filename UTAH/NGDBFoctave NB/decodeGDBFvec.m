function [d, failed, S, Sb, E]  = decodeGDBFvec(y, H, N, R, M, T, w, theta, nsigma)

d = sign(y);
l=0;
while l<T
  Sb = mod((0.5*(1-d))*H',2);
  S = 1-2*Sb;
  E = d.*y + w*S*H + nsigma*randn(R,N);
  flipdx = find(E<theta);
  d(flipdx) = -d(flipdx);
  failed=sum(Sb');
  l=l+1;
end