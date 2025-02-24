
#include <math.h>
#include "./include/struct.h"
#include "./include/init.h"
#include "./include/tools.h"
#include "./include/channel.h"

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
    sigma = sqrt(1.0/(pow(10,EbN/10.0))); // considering SNR
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


void ModelChannel_AWGN_64_CSK(csk_t *csk,code_t *code, decoder_t *decoder, int **NBIN, float EbN, int *init_rand)
{
    const int N = code->N;
    int n,k,g,q;
    float u,v,sigma;
    float TMP[code->GF];
    int som;

    float **NoisyBin = calloc(csk->PNsize,sizeof(float *));
    for (q=0; q<csk->PNsize; q++) NoisyBin[q] = calloc(2,sizeof(float));
    /* Binary-input AWGN channel : */

    int i;

    float modulation[code->GF][2];



int pn_size = 64;


// test of build_nb_pn2 (increase minimal distance) (it works!)
// int CCSK64[64]={11, 8, 18, 56, 37, 51, 61, 5, 42, 7, 4, 50, 17, 3, 59, 20, 64, 53, 12, 24, 38, 1, 13, 2, 23, 39, 26, 30, 62, 31, 44, 6, 47, 45, 19, 33, 15, 40, 49, 10, 60, 54, 14, 21, 29, 28, 46, 55, 43, 25, 48, 63, 36, 9, 16, 35, 22, 32, 57, 41, 52, 58, 34, 27
//};

// min=3.16
// int CCSK64[64]={4, 25, 14, 26, 15, 16, 43, 1, 28, 39, 38, 58, 19, 59, 48, 6, 54, 27, 17, 3, 32, 46, 21, 13, 20, 50, 41, 61, 12, 7, 55, 9, 22, 56, 49, 23, 33, 52, 18, 5, 24, 40, 10, 11, 47, 51, 53, 57, 34, 62, 30, 42, 37, 36, 35, 8, 63, 31, 29, 44, 2, 45, 64, 60
//};

//test
//  int CCSK64[64]={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64
//  };

// min=3.222
//const int CCSK64[64]={  36, 58, 54, 28, 60, 35, 22, 45, 39, 57, 4, 6, 27, 49, 52, 1, 14, 59, 42, 19, 47, 31, 61, 30, 5, 12, 33, 23, 7, 56, 44, 55, 3, 16, 0, 11, 62, 21, 46, 25, 41, 38, 51, 37, 13, 15, 18, 10, 2, 53, 32, 63, 17, 26, 24, 48, 34, 40, 20, 43, 9, 50, 29, 8
//  };
// min=5.2 considering square distance
//int CCSK64[64]={21, 23, 9, 42, 60, 30, 34, 5, 63, 3, 18, 36, 51, 35, 20, 53, 22, 11, 62, 28, 52, 4, 58, 37, 55, 50, 10, 15, 7, 0, 25, 27, 2, 41, 29, 38, 44, 47, 19, 6, 13, 8, 39, 17, 24, 43, 26, 46, 49, 54, 45, 61, 14, 40, 16, 31, 48, 1, 32, 12, 33, 59, 57, 56   };

// min = 5.26
//const int CCSK64[64]={ 58, 8, 56, 36, 0, 53, 45, 46, 14, 57, 54, 62, 16, 3, 22, 44, 41, 9, 42, 26, 13, 37, 31, 19, 12, 34, 55, 63, 32, 25, 21, 39, 40, 6, 4, 23, 30, 15, 47, 61, 11, 17, 59, 5, 1, 51, 48, 33, 27, 35, 29, 18, 52, 2, 28, 38, 49, 60, 50, 7, 24, 20, 43, 10};

// min = 5.48
//const int CCSK64[64]={44, 54, 11, 46, 17, 0, 43, 21, 48, 56, 3, 20, 15, 63, 39, 27, 52, 18, 26, 4, 61, 13, 29, 59, 32, 24, 33, 22, 25, 37, 57, 62, 14, 30, 53, 40, 6, 31, 51, 55, 58, 7, 5, 42, 8, 23, 36, 49, 41, 19, 45, 50, 38, 9, 10, 16, 34, 47, 35, 12, 28, 2, 60, 1};

// min=7.8 from EB  best!!!!
//int CCSK64[64]={48, 25, 51, 37, 13, 55, 36, 57, 41, 56, 33, 39, 45, 18, 14, 35, 22, 58, 24, 10, 63, 6, 1, 27, 53, 8, 43, 54, 59, 20, 60, 19, 40, 12, 21, 42, 9, 46, 31, 50, 64, 38, 62, 44, 34, 26, 15, 3, 47, 2, 5, 16, 61, 11, 17, 29, 7, 32, 28, 49, 52, 4, 23, 30};

// min=8.6 from EB
//int CCSK64[64]={63, 42, 39, 41, 53, 51, 29, 19, 26, 55, 32, 14, 9, 2, 40, 45, 48, 1, 21, 50, 24, 25, 38, 60, 37, 4, 64, 28, 3, 33, 15, 56, 16, 52, 12, 7, 46, 5, 13, 57, 20, 34, 49, 62, 10, 58, 11, 43, 8, 31, 17, 22, 61, 47, 36, 30, 18, 44, 54, 35, 27, 6, 23, 59};

// min=6.47 from CM
//int CCSK64[64]={28, 17, 24, 40, 5, 51, 3, 16, 39, 63, 1, 18, 0, 62, 31, 52, 26, 27, 46, 59, 7, 2, 29, 44, 36, 61, 55, 48, 30, 6, 19, 34, 47, 45, 54, 41, 32, 20, 43, 21, 15, 38, 50, 23, 33, 13, 58, 9, 25, 49, 42, 60, 22, 35, 4, 56, 12, 53, 11, 37, 14, 10, 8, 57};

//min=9.25
//int CCSK64[64]={34, 36, 21, 16, 45, 46, 50, 9, 20, 12, 55, 2, 10, 48, 6, 61, 3, 19, 8, 27, 15, 24, 51, 32, 1, 18, 14, 29, 52, 54, 22, 49, 39, 43, 38, 11, 30, 28, 62, 40, 58, 47, 44, 60, 23, 26,  5, 17, 25, 63, 13,  7, 42, 53, 64,  4, 35, 37, 33, 31, 41, 57, 56, 59
//min=10.55 alexandru
//int CCSK64[64]={21, 57, 63, 25, 38, 51, 7, 26, 16, 8, 43, 5, 36, 54, 2, 3, 55, 61, 14, 45, 37, 1, 41, 56, 15, 49, 22, 52, 59, 12, 4, 19, 50, 48, 64, 13, 23, 27, 29, 18, 9, 58, 39, 60, 6, 28, 46, 62, 32, 33, 20, 11, 40, 30, 47, 31, 24, 44, 17, 42, 34, 53, 35, 10 };
//10.5525,14.7053,24.3705
//int CCSK64[64]={24,44,17,42,34,53,35,10,21,57,63,29,26,16,8,43,37,1,41,56,15,49,18,9,58,39,60,6,28,46,62,32,33,20,11,40,30,23,47,31,25,38,51,7,22,52,59,12,4,19,50,48,64,13,5,36,54,2,3,55,61,14,45,27};

//max(ccorr(2:64))=27.7
//int CCSK64[64]={35,15,14,62,24,53,34,2,49,32,42,13,57,25,59,51,63,8,61,46,29,5,16,28,47,18,27,0,20,6,36,45,38,23,31,11,39,54,26,55,4,48,17,10,60,30,44,43,12,3,22,50,1,58,9,33,21,40,7,41,52,19,56,37};


//QAM
// 0.47
//int CCSK64[64]={29,42,55,44,35,45,18,13,46,59,41,2,27,56,62,14,6,9,25,39,32,15,40,0,26,54,53,19,47,23,7,61,48,22,43,58,8,60,16,51,4,63,50,1,49,12,31,20,30,34,57,24,33,52,5,11,21,36,17,38,10,3,28,37};

//0.47  1.05  2.1  2.9  4.47 5.3
//int CCSK64[64]={2,22,46,40,35,60,41,8,10,29,38,31,63,54,34,15,43,13,4,55,24,37,5,14,7,18,52,48,1,25,53,26,11,58,6,20,12,62,51,21,9,59,32,30,28,42,56,39,3,50,49,45,19,57,0,44,61,33,23,47,17,16,27,36};




//int CCSK64[64]={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,64};


//64PSK
//int CCSK64[64]={17,20,39,59,32,9,8,53,50,54,45,62,10,51,12,63,33,36,55,11,48,25,24,5,2,6,61,14,26,3,28,15,49,52,7,27,64,41,40,21,18,22,13,30,42,19,44,31,1,4,23,43,16,57,56,37,34,38,29,46,58,35,60,47};
//int CCSK64[64]={3,53,48,46,23,42,52,17,44,29,56,54,63,18,59,25,19,5,64,62,39,58,4,33,60,45,8,6,15,34,11,41,35,21,16,14,55,10,20,49,12,61,24,22,31,50,27,57,51,37,32,30,7,26,36,1,28,13,40,38,47,2,43,9};
//epicycloid 4 cusps
int CCSK64[64]={0,63,58,33,36,51,62,53,8,39,2,9,44,27,6,29,16,15,10,49,52,3,14,5,24,55,18,25,60,43,22,45,32,31,26,1,4,19,30,21,40,7,34,41,12,59,38,61,48,47,42,17,20,35,46,37,56,23,50,57,28,11,54,13};


//    for(i=0; i< 64 ; i++)
//    {
//    CCSK64[i]=CCSK64[i]-1;
//   }




   // getchar();
//
//    ////compute normalization factor so that average power of one point of constellation is equal to one
//    float norm_factor=0.0;
//    for(i=0; i < code->GF ; i++)
//    {
//        norm_factor = table_64QAM[i][0]*table_64QAM[i][0] + table_64QAM[i][1]*table_64QAM[i][1]+norm_factor;//compute sum
//        //norm_factor = table_64APSK[i][0]*table_64APSK[i][0] + table_64APSK[i][1]*table_64APSK[i][1]+norm_factor;//compute sum
//    }
//    norm_factor = sqrt( code->GF / norm_factor);
//// printf(" norm_factor = %f ", norm_factor); getchar();
//
//    for(i=0; i< code->GF ; i++)
//    {
//        //modulation[i][0]=norm_factor*table_64APSK[i][0];
//        //modulation[i][1]=norm_factor*table_64APSK[i][1];
//        modulation[i][0]=norm_factor*table_64QAM[i][0];
//        modulation[i][1]=norm_factor*table_64QAM[i][1];
//
//    }



      //  for 64PSK
    for (i=0;i<code->GF; i++)
    {
        modulation[i][0]=sin(i*2*PI/64);
        modulation[i][1]=cos(i*2*PI/64);
        //modulation[i][0]=table_64PSK[i-1][0];
        //modulation[i][1]=table_64PSK[i-1][1];
    }


//// print constellation in a file
//    FILE *opfile;
//    opfile=fopen("64QAM.txt","a");
//    for(i=0; i< code->GF ; i++)
//    {
//        fprintf(opfile," %f %f ;", modulation[i][0],modulation[i][1]);
//    }
//    fclose(opfile);
//    getchar();




    sigma = sqrt(1.0/(2*pow(10,EbN/10.0)));

    for (n=0; n<N; n++)
    {
        som=0;
        for (q=0; q<6; q++)
        {
            som = som + NBIN[n][q]*pow(2,q);
        }
        //printf("\n %d \n",som );

		for(k=0; k<code->GF; k++)
        {
            TMP[k] =0.0;
        }

        for (q=0; q<csk->PNsize; q++)
        {
            for (i=0; i<2; i++)
            {
                u=My_drand48(init_rand);
                v=My_drand48(init_rand);
                /* BPSK modulation + AWGN noise (Box Muller method for Gaussian sampling) */
                NoisyBin[q][i] = modulation[CCSK64[(som+q)%pn_size]][i]+ sigma*sqrt(-2.0*log(u))*cos(2.0*PI*v)  ;
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

				for (g=0; g<csk->PNsize; g++)
                {

                    TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK64[(som+g)%pn_size]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK64[(som+g)%pn_size]][1])/(2.0*SQR(sigma));
                    //printf("%d %f \n",k, TMP[k]);

                   //complex multiplication for correlation, compute real part
                    // (a+ib)*conj(c+id)= ac + bd + i (bc -ad)
                    //TMP[k] = TMP[k]-( NoisyBin[g][0]*modulation[CCSK64[(som+g)%128]][0] + NoisyBin[g][1]*modulation[CCSK64[(som+g)%128]][1] ); //real part


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

                for (g=0; g<csk->PNsize-1; g++)
                {
                    TMP[k] = TMP[k]+SQR(NoisyBin[g][0]-modulation[CCSK64[(som+g)%pn_size]][0])/(2.0*SQR(sigma))+SQR(NoisyBin[g][1]-modulation[CCSK64[(som+g)%pn_size]][1])/(2.0*SQR(sigma));
                    //printf("%d %f \n",k, TMP[k]);
                    //complex multiplication for correlation, compute real part
                    // (a+ib)*conj(c+id)= ac + bd + i (bc -ad)
                    //TMP[k] = TMP[k]-( NoisyBin[g][0]*modulation[CCSK64[(som+g)%pn_size]][0] + NoisyBin[g][1]*modulation[CCSK64[(som+g)%pn_size]][1] );// + SQR(  NoisyBin[g][1] * modulation[CCSK64[(som+g)%128]][0] -  NoisyBin[g][0]*modulation[CCSK64[(som+g)%128]][1]        );

                }
            }
        }

        //getchar();
        for(k=0; k<code->GF; k++)
        {
            decoder->intrinsic_LLR[n][k] = TMP[k]+10*log10(2/pow(sigma,2));
//            printf("%d %.2f\n", k, decoder->intrinsic_LLR[n][k]);
        }
//        getchar();
    }

    for (q=0; q<csk->PNsize; q++) free(NoisyBin[q]);
    free(NoisyBin);
}




