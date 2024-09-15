function [info_seq, code_seq, code_seq_comp, valid_symdrom] = ...
    generate_and_encode_QAM(ZERO, h,G, add_mat, mul_mat,p, constl, gray_labels)
[M,N]=size(h);
K=N-M;
q=2^p;
info_seq = ZERO*randi([0 q-1], 1, K);

    info_seq_gray = gray_labels(info_seq+1);


code_seq_gray = gf_mat_mul(info_seq_gray,G, add_mat, mul_mat);
code_seq = gray_labels(code_seq_gray+1);

valid_symdrom = gf_mat_mul(gray_labels(code_seq+1),h', add_mat, mul_mat);

code_seq_comp = constl(code_seq_gray+1);



