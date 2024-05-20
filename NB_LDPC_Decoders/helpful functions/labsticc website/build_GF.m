% build_GF.m
% Author: Cédric Marchand
% build finite field vectors

clear all;

%GF_size=2; primitive= [1 1 1]; %GF[4]
%GF_size=4; primitive= [1 1 0 0 1]; %GF[16], primitive 1+x+x^4
% GF_size=5; primitive= [1 0 1 0 0 1]; %GF[32]
GF_size=6; primitive= [1 1 0 0 0 0 1]; %GF[64]
%GF_size=8; primitive= [1 0 1 1 1 0 0 0 1]; %GF[256]
%GF_size=10; primitive=[1 0 0 1 0 0 0 0 0 0 1]; %GF[1024]
%GF_size=11; primitive=[1 0 1 0 0 0 0 0 0 0 0 1]; %GF[2048]
%GF_size=12; primitive=[1 1 0 0 1 0 1 0 0 0 0 0 1]; % GF[4048]

GF_vector= zeros(1,GF_size+1);
GF_vectors=zeros(2^GF_size,GF_size);
GF_vectors(2,1)=1;
GF_vector(1)=1;


%%

for i=1:2^GF_size-2

    GF_vector=[0 GF_vector(1:end-1)]; %shif right (multiplication by x)

    if (GF_vector(end)==1) % compute modulo primitive
        GF_vector= xor(GF_vector,primitive);
    end

    GF_vectors(i+2,:)=GF_vector(1:end-1);

end

GF_vectors


