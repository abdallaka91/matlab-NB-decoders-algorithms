
#include <math.h>
#include "./include/struct.h"
#include "./include/init.h"
#include "./include/tools.h"
#include "./include/channel.h"
void ModelChannel_AWGN_64_CSK(code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand)
{
    const int N = code->N;
    int n,k,g,q;
    float u,v,sigma;
    float TMP[code->GF];
    int som;
    int pn_size = 64;
    float **NoisyBin = calloc(pn_size,sizeof(float *));
    for (q=0; q<pn_size; q++) NoisyBin[q] = calloc(2,sizeof(float));
    /* Binary-input AWGN channel : */

    int i;

    float modulation[code->GF][2];




int CCSK64[64]={0,63,58,33,36,51,62,53,8,39,2,9,44,27,6,29,16,15,10,49,52,3,14,5,24,55,18,25,60,43,22,45,32,31,26,1,4,19,30,21,40,7,34,41,12,59,38,61,48,47,42,17,20,35,46,37,56,23,50,57,28,11,54,13};

    for (i=0;i<code->GF; i++)
    {
        modulation[i][0]=sin(i*2*PI/64);
        modulation[i][1]=cos(i*2*PI/64);
    }
    
    sigma = sqrt(1.0/(2*pow(10,EbN/10.0)));

    for (n=0; n<N; n++)
    {
        som=0;
        for (q=0; q<6; q++)
        {
            som = som + NBIN[n][q]*pow(2,q);
        }
		for(k=0; k<code->GF; k++)
        {
            TMP[k] =0.0;
        }

        for (q=0; q<pn_size; q++)
        {
            for (i=0; i<2; i++)
            {
                u=My_drand48(init_rand);
                v=My_drand48(init_rand);
                /* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
                NoisyBin[q][i] = modulation[CCSK64[(som+q)%pn_size]][i]+ sigma*sqrt(-2.0*log(u))*cos(2*PI*v)  ;
            }
        }

        if (n<N) // per symbol puncturing
        {
            for(k=0; k<code->GF; k++)
            {
				som=0;
				for (q=0; q<6; q++)
				{
					som = som + BinGF_64[k][q]*pow(2,q);
				}
				for (g=0; g<pn_size; g++)
                {
                    TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK64[(som+g)%pn_size]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK64[(som+g)%pn_size]][1])/(2.0*SQR(sigma));
                }
            }
        }
        else
        {

            for(k=0; k<code->GF; k++)
            {

				som=0;
				for (q=0; q<6; q++)
				{
					som = som + BinGF_64[k][q]*pow(2,q);
				}

                for (g=0; g<pn_size-1; g++)
                {
                    TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK64[(som+g)%pn_size]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK64[(som+g)%pn_size]][1])/(2.0*SQR(sigma));
                }
            }
        }

        //getchar();
        for(k=0; k<code->GF; k++)
        {
            decoder->intrinsic_LLR[n][k] = TMP[k]+10*log10((2/pow(sigma,2)));
//            printf("%d\t%.2f\n", k, decoder->intrinsic_LLR[n][k]);
//            getchar();
        }
    }

    for (q=0; q<pn_size; q++) free(NoisyBin[q]);
    free(NoisyBin);
}


void ModelChannel_AWGN_BPSK (code_t *code, decoder_t *decoder, table_t *table, int **NBIN, float EbN,int *init_rand)
{
    const int N = code->N;
    const int logGF = code->logGF;
    int n,k,g,q;
    float u,v,sigma;
    float TMP[4096];

    float **NoisyBin = calloc(N,sizeof(float *));
    for (n=0; n<N; n++) NoisyBin[n] = calloc(logGF,sizeof(float));

    /* Binary-input AWGN channel : */

//    sigma = sqrt(1.0/(2*code->rate*pow(10,EbN/10.0)));  // considering EbNo
    sigma = sqrt(1.0/(2*pow(10,EbN/10.0))); // considering SNR
    for (n=0; n<N; n++)
    {
        for (q=0; q<logGF; q++)
        {

            u=My_drand48(init_rand);
            while(u==0.0)
			{
                u=My_drand48(init_rand);
			}
            v=My_drand48(init_rand);
            /* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
            NoisyBin[n][q] = BPSK(NBIN[n][q]) + sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v);
        }

    }


    /* Compute the Log intrinsic_LLR Ratio messages */
    for (n=0; n<N; n++)
    {
        for (g=0; g<code->GF; g++)
        {
            TMP[g]=0.0;
            for (q=0; q<logGF; q++)
            {
                //TMP[g] = TMP[g] + SQR(NoisyBin[n][q]-BPSK(table->BINGF[g][q]))/(2.0*SQR(sigma));
                TMP[g] = TMP[g] - NoisyBin[n][q]*BPSK(table->BINGF[g][q]);
            }

        }

        for(k=0; k<code->GF; k++)
        {
            decoder->intrinsic_LLR[n][k] = TMP[k];
        }


    }

    for (n=0; n<N; n++) free(NoisyBin[n]);
    free(NoisyBin);
}


