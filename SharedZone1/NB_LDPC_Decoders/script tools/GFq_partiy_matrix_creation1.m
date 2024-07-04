clear

pth4 = (fullfile(fileparts(pwd), 'related_variables\alists\'));
q = 32;
dc = 26;
dv = 4;
eye_size = 31;%384/6;

M = dv*eye_size;
N = dc*eye_size;
p = log2(q);

H = zeros(M, N);

I = eye(eye_size,eye_size);

mat1 = zeros(dv, dc);
i1 = randperm(eye_size, dc);
mat1(1,:) = i1;

for  i = 2 : dv

    s=1;
    while s~=0
        per1 = randperm(eye_size, dc);
        s = sum(sum(per1==mat1(1:i-1,:)));
    end
    mat1(i,:) = per1;
end

mat2 = zeros(dv, dc);

per1 = randi([1 q-1], 1, dc);
mat2(1,:) = per1;
for  i = 2 : dv

    s=1;
    while s~=0
        per1 = randi(q-1, 1, dc);
%         per1 = randperm(q-1, N_eye);
        s = sum(sum(per1==mat2(1:i-1,:)));
    end
    mat2(i,:) = per1;
end


for  i = 1 : dv
    for j = 1 : dc
        c1 = circshift(eye(eye_size), mat1(i, j));
        H((i-1)*eye_size+1:i*eye_size, (j-1)*eye_size+1:j*eye_size) = c1;
    end
    
end

H1 = H;

for j = 1 : N
    i1 = find(H1(:,j));
    NOK = true;
    while NOK
        c1 = randperm(q-1, dv);
        NOK = false;
        for k = 1 : length(i1)
            if sum(c1(k)==H1(i1(k),1:j-1))
                NOK = true;
                break
            end
        end
    end

    for i = 1 : dv
        H1((i-1)*eye_size+1:i*eye_size, j) = H1((i-1)*eye_size+1:i*eye_size, j)*c1(i);
    end
end

HHh = zeros(M, dc);
HHv = zeros(dv, N);

for j = 1 : N
    HHv(:, j) = H1(find(H1(:,j)~=0), j);
end
for i = 1 : M
    HHh(i,:) = H1(i, find(H1(i,:)~=0));
end

h = H1;

fl1 = ['generated_' num2str(M) 'x' num2str(N) '_GF' num2str(q) '.mat'];
fll_nm = fullfile(pth4, fl1);
save(fll_nm,'h')








