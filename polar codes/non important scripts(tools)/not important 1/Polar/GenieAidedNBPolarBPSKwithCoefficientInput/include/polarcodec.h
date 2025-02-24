#include "struct.h"
#include "init.h"
#include "tools.h"
#include "channel.h"

void Encoder (code_t *code, table_t *table, int *KSYMB, int *NSYMB, int **NBIN);
void SC_CN(float **L, int arr_size, int* coeff, table_t *table, int GF, float **result);
void SC_VN(float **L, int *c, int arr_size, int *coeff, table_t *table, int GF, float **res);
void GenieAidedDecoderSC (decoder_t *decoder, code_t *code, table_t *table, int* KSYMB);
