function [SER, BER] = SER_BER(seq,rec_seq,p, BER, SER)

L = length(seq);
for i = 1 : L
    if seq(i)~=rec_seq(i)
        SER = SER+1;
        s1 = dec2bin(seq(i),p);
        s2 = dec2bin(rec_seq(i),p);
        seq_d = double(s1);
        rec_seq_d = double(s2);
        num_diff_bit = sum(seq_d ~= rec_seq_d);
        BER = BER + num_diff_bit;
    end
end
end
