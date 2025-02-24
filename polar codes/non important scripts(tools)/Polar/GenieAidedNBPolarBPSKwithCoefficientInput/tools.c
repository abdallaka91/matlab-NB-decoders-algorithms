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

void Decision( int *decision,float **APP,int N,int GF)
{
    int n,g,ind;
    float min;

    for (n=0; n<N; n++)
    {
        //max = APP[n][0];
        min=+1e5;
        ind=0;
        for (g=0; g<GF; g++)
            if (APP[n][g]<min)
            {
                min=APP[n][g];
                ind=g;
            }
        decision[n]=ind;
            printf("%d ", decision[n]);
    }
    getchar();

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

void computeErrProb(code_t *code, decoder_t *decoder, int *KSYMB){
    float temp;
    temp=0;
    for(int n=0;n<code->N;n++){
        temp=0;
        for(int q=0;q<code->GF;q++){
            temp+=exp(-decoder->extrinsic_LLR[n][q]);
        }
        code->ErrorProb[n]+=1-(exp(-decoder->extrinsic_LLR[n][KSYMB[n]])/temp);
    }
    return;
}

void computeEntropy(code_t *code, decoder_t *decoder){
    float temp;
    float pdf[code->GF];
    for(int n=0;n<code->N;n++){
        temp=0;
        for(int q=0;q<code->GF;q++){
            temp+=exp(-decoder->extrinsic_LLR[n][q]);
        }

        for(int q=0;q<code->GF;q++){
            pdf[q]=(exp(-decoder->extrinsic_LLR[n][q])/temp);
            if(pdf[q]!=0)
                code->Entropy[n]+=-1*pdf[q]*log(pdf[q])/log(code->GF);
        }
    }
    return;
}

void printSequence(code_t *code, char *filename, int nb){
    FILE *opfile;
    opfile=fopen(filename,"a");
    if ((opfile)==NULL)
    {
        printf(" \n !! file not found \n ");
        return;
    }
    insertion_sort_desc(code->Entropy, code->reliability_sequence, code->N);
    for(int i=0;i<code->N;i++){
        fprintf(opfile, "%d ", code->reliability_sequence[i]);
    }
    fprintf(opfile,"\n\n");
        for(int i=0;i<code->n;i++){
        for(int j=0;j<code->N/2;j++){
            fprintf(opfile, "%d " , code->polar_coeff[i][j]);
        }
    }

    fprintf(opfile,"\n\n");

    insertion_sort_desc(code->ErrorProb, code->reliability_sequence, code->N);
    for(int i=0;i<code->N;i++){
//        printf("%d %.2f\n", code->reliability_sequence[i], code->ErrorProb[code->reliability_sequence[i]]);
        fprintf(opfile, "%d ", code->reliability_sequence[i]);
    }
    fprintf(opfile,"\n\n");
    for(int i=0;i<code->n;i++){
        for(int j=0;j<code->N/2;j++){
            fprintf(opfile, "%d " , code->polar_coeff[i][j]);
        }
    }
    fprintf(opfile,"\n\n");

    for(int i=0;i<code->N;i++){
        fprintf(opfile,"%g ", fabs(code->Entropy[i])/(float) nb);
    }
    fprintf(opfile,"\n\n");

    for(int i=0;i<code->N;i++){
        fprintf(opfile,"%g ", code->ErrorProb[i]/(float)nb);
    }

    fprintf(opfile,"\nNb. Frames: %d\n\n", nb);

//    getchar();
    fclose(opfile);
}

void insertion_sort_desc( float *LLRs, int *indices, int array_size){
  int i , j ;

  for ( i = 0 ; i < array_size ; i++ )
  {
    for ( j = i ; ( j > 0 ) && ( LLRs [ indices [ j - 1 ] ] < LLRs [ i ] ) ; j-- )
      indices [ j ] = indices [ j - 1 ] ;

        indices [ j ] = i ;
  }

  return ;
}

