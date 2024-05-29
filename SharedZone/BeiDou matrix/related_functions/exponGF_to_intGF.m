function y = exponGF_to_intGF(x, p, polynml)

GF_size = p;
primitive = cell(12,1);
primitive{1}= []; %GF[4]
primitive{2}= [1 1 1]; %GF[4]
primitive{3}= [];
primitive{4}= [1 1 0 0 1]; %GF[16], primitive 1+x+x^4
primitive{5}= []; %GF[32]
primitive{6}= [1 1 0 0 0 0 1]; %GF[64]
primitive{3}= [];
primitive{8}= [1 0 1 1 1 0 0 0 1]; %GF[256]
primitive{9}= [];
primitive{10}= [1 0 0 1 0 0 0 0 0 0 1]; %GF[1024]
primitive{11}= [1 0 1 0 0 0 0 0 0 0 0 1]; %GF[2048]
primitive{12}= [1 1 0 0 1 0 1 0 0 0 0 0 1]; % GF[4096]


switch nargin
    case 2
        primitive_p = primitive{p};
    case 3
        primitive_p = polynml;
end


GF_vector= zeros(1,GF_size+1);
GF_vectors=zeros(2^GF_size,GF_size);
GF_vectors(2,1)=1;
GF_vector(1)=1;

for i=1:2^GF_size-2

    GF_vector=[0 GF_vector(1:end-1)]; %shif right (multiplication by x)

    if (GF_vector(end)==1) % compute modulo primitive
        GF_vector= xor(GF_vector,primitive_p);
    end

    GF_vectors(i+2,:)=GF_vector(1:end-1);

end

y = zeros(size(x));
for i = 1 : length(x)
    if x(i) < 0
        y(i) = 0;
    elseif x(i) == 0
        y(i) = 1;
    else
        y(i) = bin2dec(flip(char(GF_vectors(x(i)+2,:)+48)));
    end
end
end
