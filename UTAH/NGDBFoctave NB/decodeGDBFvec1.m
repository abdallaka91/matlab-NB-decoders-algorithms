function [d, failed, S, Sb, E]  = decodeGDBFvec1(y, H, N, R, M, T, w, theta, nsigma)

d = sign(y);
l=0;
while l<T
    Sb = mod((0.5*(1-d))*H',2);
  S = 1-2*Sb;
  SH = S*H ;
Sh = zeros(1,N);
  for i = 1 : M

            idx1 = find(H(i,:));
           for j = idx1
               Sh(j) = Sh(j)+(1-2*Sb(i));
           end

  end

  E = d.*y + w*SH + nsigma*randn(R,N);
  flipdx = find(E<theta);
  d(flipdx) = -d(flipdx);
  failed=sum(Sb');
  l=l+1;
end
