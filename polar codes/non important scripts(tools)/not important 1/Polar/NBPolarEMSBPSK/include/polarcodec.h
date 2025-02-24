#include "struct.h"
#include "init.h"
#include "tools.h"
#include "channel.h"

void Encoder (code_t *code, table_t *table, int *KSYMB, int *NSYMB, int **NBIN);
void DecoderSC (decoder_t *decoder, code_t *code,table_t *table, int *decoded_msg);
void SC_CN(float **L, int arr_size, int* coeff, table_t *table, int GF, int nm, float offset, float **result);
void SC_VN(float **L, int *c, int arr_size, int *coeff, table_t *table, int GF, float **res);
