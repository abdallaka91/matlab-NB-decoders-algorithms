#include "struct.h"
#include "init.h"
#include "tools.h"
#include "channel.h"

void Encoder (code_t *code, table_t *table, int *KSYMB, int *NSYMB, int **NBIN, decoder_t *decoder);
void DecoderSC (decoder_t *decoder, code_t *code, table_t *table, float err[code->N/2][code->GF], int dep_ind, int node_ind);
void SC_CN(float **L, int arr_size, int* coeff, table_t *table, int GF, float **result);
void SC_VN(float **L, int *c, int arr_size, int *coeff, table_t *table, int GF, float **res);
