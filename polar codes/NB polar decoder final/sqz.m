function y = sqz(x,d)
sz = size(x);
nml1 = numel(sz);
if nml1<3
    y = x;
else
    sz1 = sz([1:d-1 d+1:nml1]);
    y = zeros(sz1);
    for i = 1 : sz1(1)
        for j = 1 : sz1(2)
            if d==1
                y(i,j) = x(1,i, j);
            elseif d==2
                y(i,j) = x(i, j);
            else
                y(i,j) = x(i, j, 1);
            end
        end
    end
end




