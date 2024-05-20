clear
A=[4 6 7 8];  % information bits indexes
u_A=[1 1 1 1];  % information bits values
a(A)= u_A;

l0=a;

l0(1)=xor(a(1),a(5));
l0(2)=xor(a(2),a(6));
l0(3)=xor(a(3),a(7));
l0(4)=xor(a(4),a(8));

l1 = l0;

l1(1)=xor(l0(1),l0(3));
l1(2)=xor(l0(2),l0(4));
l1(5)=xor(l0(5),l0(7));
l1(6)=xor(l0(6),l0(8));

l2 = l1;

l2(1)=xor(l1(1),l1(2));
l2(3)=xor(l1(3),l1(4));
l2(5)=xor(l1(5),l1(6));
l2(7)=xor(l1(7),l1(8));

l2
