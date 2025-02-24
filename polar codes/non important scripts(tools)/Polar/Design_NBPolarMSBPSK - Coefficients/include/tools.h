#include "./struct.h"
#include <stdlib.h>
#include <stdio.h>



float My_drand48(int *initialise);
int Bin2GF(int *U,int GF,int logGF,int **BINGF);
void RandomBinaryGenerator (int K,int GF,int logGF, int **KBIN, int *KSYMB,int **BINGF, int *init_rand);
void F_MS (float *input_LLR1, float *input_LLR2, int h, table_t *table, int GF, float *output_LLR);
void G_MS (float *input_LLR1, float *input_LLR2, int decision, int h, table_t *table, int GF, float *output_LLR);
void insertion_sort( float *LLRs, int *indices, int array_size);
float computeErrProb(float *LLR, int real_symb, int GF);
void normalizeLLR(decoder_t *decoder, code_t *code);
float computeEntropy(float *LLR, int GF);
