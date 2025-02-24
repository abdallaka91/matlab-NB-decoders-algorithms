#include <string.h>
#include <math.h>

#include "./include/struct.h"
#include "./include/init.h"
#include "./include/tools.h"
float My_drand48(int *initialise)
{

#if defined(_WIN32) || defined(__CYGWIN__) || defined(__MINGW32__)

    static int s1, s2;
    int k, Z;

    if ( *initialise == -1 )
    {
        s1 = (int)(( rand() % 2147483562 ) + 1);
        s2 = (int)(( rand() % 2147483398 ) + 1);
        *initialise = 0;
    }

    k = s1/53668;
    s1 = 40014*(s1 - k*53668) - k*12211;
    if (s1 < 0)
        s1 += 2147483563;

    k = s2/52774;
    s2 = 40692*(s2 - k*52774) - k*3791;
    if (s2 < 0)
        s2 += 2147483399;

    Z = s1 - s2;
    if (Z < 1)
        Z += 2147483562;

    return(Z/2147483563.0);

#elif __linux__

    return(drand48());

#endif


}

int Bin2GF(int *U,int GF,int logGF,int **BINGF)
{
    int k;

    for (k=0; k<GF; k++)
    {
        if (memcmp(U,BINGF[k],sizeof(int)*logGF)==0) break;
    }
    return(k);
}

void RandomBinaryGenerator (int K,int GF,int logGF, int **KBIN, int *KSYMB,int **BINGF, int *init_rand)
{
    int k,q;

    /* Random and (bitwise) uniformly distributed information symbols */
    for (k=0; k<K; k++)
    {
        for (q=0; q<logGF; q++)
            KBIN[k][q]=floor(My_drand48(init_rand)*1.9999); // avoid the case 2

        KSYMB[k]=Bin2GF(KBIN[k],GF,logGF,BINGF);
     }
}

void F_MS (float *input_LLR1, float *input_LLR2, int h, table_t *table, int GF, float *output_LLR){
    float temp;
    int u_theta, u_phi;
    float min=10e5;
    // CN update
    for(u_theta=0;u_theta<GF;u_theta++){
        output_LLR[u_theta]=10e5;
    }
    for(u_theta=0;u_theta<GF;u_theta++){
        for(u_phi=0;u_phi<GF;u_phi++){
            temp=input_LLR1[table->ADDGF[u_theta][u_phi]]+input_LLR2[table->MULGF[h][u_phi]];
            if(output_LLR[u_theta]>temp){
                output_LLR[u_theta]=temp;
            }

            if(output_LLR[u_theta]<min){
                min=output_LLR[u_theta];
            }
        }
    }
    for(int q=0;q<GF;q++){
        output_LLR[q]=output_LLR[q]-min;
    }
    return;
}

void G_MS (float *input_LLR1, float *input_LLR2, int decision, int h, table_t *table, int GF, float *output_LLR){
    float min=10e5;

    for(int phi=0;phi<GF;phi++){
        output_LLR[phi]=input_LLR1[table->ADDGF[decision][phi]]+input_LLR2[table->MULGF[h][phi]];
        if(min>output_LLR[phi]){
            min=output_LLR[phi];
        }
    }
    // normalization
    for(int q=0;q<GF;q++){
        output_LLR[q]=output_LLR[q]-min;
    }
    return;
}

void insertion_sort( float *LLRs, int *indices, int array_size){
  int i , j ;

  for ( i = 0 ; i < array_size ; i++ )
  {
    for ( j = i ; ( j > 0 ) && ( LLRs [ indices [ j - 1 ] ] > LLRs [ i ] ) ; j-- )
      indices [ j ] = indices [ j - 1 ] ;

        indices [ j ] = i ;
  }

  return ;
}

float computeErrProb(float *LLR, int real_symb, int GF){
    float temp;
    temp=0;
    for(int q=0;q<GF;q++){
        temp+=exp(-LLR[q]);
    }
     return (1-(exp(-LLR[real_symb])/temp));
}

void normalizeLLR(decoder_t *decoder, code_t *code){
    float min;
    for(int n=0;n<code->N;n++){
        min=10e5;
        for(int q=0;q<code->GF;q++){
            if(min>decoder->intrinsic_LLR[n][q]){
                min=decoder->intrinsic_LLR[n][q];
            }
        }

        for(int q=0;q<code->GF;q++){
            decoder->intrinsic_LLR[n][q]=decoder->intrinsic_LLR[n][q]-min;
        }
    }
    return;
}

float computeEntropy(float *LLR, int GF){
    float entropy=0;
    for(int q=0;q<GF;q++){
        if(LLR[q]!=0)
            entropy+=-1*log(exp(-LLR[q]))*exp(-LLR[q])/log(GF);
    }
     return entropy;
}
