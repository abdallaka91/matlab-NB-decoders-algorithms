
#include "./struct.h"
#include "./init.h"

#define SQR(A) ((A)*(A))
#define BPSK(x) (1-2*(x))
#define PI 3.1415926536

void ModelChannel_AWGN_BPSK (code_t *code, decoder_t *decoder, table_t *table, int **NBIN, float EbN,int *init_rand);
void ModelChannel_AWGN_64_CSK(csk_t *csk,code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand);
