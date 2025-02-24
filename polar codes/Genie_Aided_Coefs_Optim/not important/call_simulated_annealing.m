clear
q=64;

% 16	10	39	25	43	61	49	22	31	51	2	4	8	7	56	27	5	54	30	3	44	33	34	24	58	47	41	53	55	20	29	14	46	23	50	40	19	9	28	35	1	60	21	15	36	42	63	52	57	26	0	13	48	6	37	62	12	32	38	59	11	45	18	17
bestOrder=1+[56	1	20	55	31	43	2	44	51	29	46	26	37	6	9	48	45	14	27	34	54	24	49	5	0	52	7	41	42	19	28	63	10	36	13	57	17	50	39	30	22	47	32	3	60	8	59	21	35	23	62	16	4	61	40	11	25	58	53	12	15	33	18	38];
%%
bestOrder = randperm(q);
%%
maxIter=10000;
initialTemp=20;
finalTemp = 6;
coolingRate= exp(log(finalTemp/initialTemp)/maxIter);
wght= -[60 20 20.1 2 3 3];
[bestOrder, sp_frwd] = simulatedAnnealing(q, @cost_func, maxIter, initialTemp, coolingRate, bestOrder, wght);
bestOrder0=bestOrder-1;