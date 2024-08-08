clear
pth1 = (fullfile(pwd, 'related_functions'));
pth3 = (fullfile(pwd, 'related_variables/GF_arithm'));
addpath(pth1);

p=6;
M=44;
N=88;
K = N-M;
n_k = N-K;

NC = 4;%in the paper they are written 4 columns of 11 lines each, and each line contains 4 values, see the paper and compare the values
NL = 11;%
dc = 4;%


values0 = [30	24	1	44	24	1	44	30	40	32	61	18	53	24	1	44	51	60	35	13	18	15	32	61	15	6	1	45	30	24	1	44	6	1	45	15	45	15	6	1	1	45	15	6	1	44	53	24	24	1	44	53	44	30	24	1	34	33	45	36	55	9	34	3	1	44	53	24	61	47	20	8	53	24	1	44	15	6	1	45	13	18	60	35	45	15	6	1	24	1	44	53	37	32	52	47	44	53	24	1	39	36	34	33	44	35	31	50	12	25	36	14	15	35	46	56	53	24	1	44	1	44	53	24	24	1	44	30	44	30	24	1	15	6	1	45	30	24	1	44	2	50	22	14	33	42	14	5	34	3	55	9	44	35	61	50	15	6	1	45	45	15	6	1	1	44	30	24	6	1	45	15	1	44	53	24];
idxt = [14	35	56	70	11	29	55	73	13	39	53	69	15	34	57	71	1	27	45	54	23	41	63	87	2	20	46	68	6	24	50	61	2	26	61	79	9	33	59	77	4	30	48	74	22	42	59	76	12	38	52	68	23	43	58	77	19	21	63	64	11	25	65	82	17	39	44	75	9	35	49	72	19	29	66	84	13	36	56	82	17	43	67	81	22	40	62	86	3	21	47	69	10	24	64	83	0	37	70	86	5	31	49	75	4	40	53	84	5	41	52	85	18	28	67	85	0	26	44	55	10	28	54	72	7	30	50	81	1	36	71	87	16	38	45	74	8	34	48	73	8	32	58	76	12	37	57	83	6	31	51	80	15	33	47	79	16	42	66	80	7	25	51	60	3	27	60	78	14	32	46	78	18	20	62	65];


q = 2^p;
fl_nm = ['arith_' num2str(q) '.mat'];
load(fullfile(pth3, fl_nm));


values = values0; %
% for i = 1 : length(valt0)
%     valt(i) = exponGF_to_intGF(valt0(i), p, polynom); % if polynom is not provided the function has its defaults
% end


H44_88_idx = zeros(NC*NL, dc);
H44_88_coef = zeros(NC*NL, dc);
idxtt = reshape(idxt,dc*NC, NL )';
valtt = reshape(values,dc*NC, NL )';
for i = 1 : NC
    H44_88_idx((i-1)*NL+1:i*NL,:) = idxtt(:,(i-1)*dc+1:i*dc);
    H44_88_coef((i-1)*NL+1:i*NL,:) = valtt(:,(i-1)*dc+1:i*dc);
end

H =zeros(M,N);
for i = 1 : M
    H(i, H44_88_idx(i,:)+1) = H44_88_coef(i,:);
end

h = H;
fl1 = 'BeiDou_44_bb_GF64.mat';
fll_nm = fullfile(fl1);
save(fll_nm,'h')

h01 = double(H > 0);
dv = sum(h01);

[G,Hbar] = Generator_matrix_G_from_full_rank_H(H, add_mat, mul_mat, div_mat);

%% test on the example of the paper

bin_msg = ['001010';'110010';'010011';'100001';'001010';'100110';'010000';'101001';'101100';'101111';...
'011100';'000101';'001110';'111010';'001001';'110100';'100010';'111111';'000101';'011100';...
'000110';'111101';'000000';'110001';'110100';'110111';'000101';'011001';'010000';'110011';...
'011011';'111010';'001011';'010000';'001001';'001000';'110111';'100101';'100011';'001001';...
'110110';'100111';'010110';'100000'];
int_msg = bin2dec(bin_msg)';
int_cde = gf_mat_mul(int_msg, G, add_mat, mul_mat);
cde_min = string((dec2bin(int_cde, p)))'
symdrom = gf_mat_mul(int_cde,H', add_mat, mul_mat);