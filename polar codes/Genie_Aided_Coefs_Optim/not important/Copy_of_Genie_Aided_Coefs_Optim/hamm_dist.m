function hamming_distance=hamm_dist(a,b)

hamming_distance = sum(dec2bin(bitxor(a, b)) == '1');
