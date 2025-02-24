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

void F_EMS (float *input_LLR1, float *input_LLR2, int h, table_t *table, int GF, int nm, float offset, float *output_LLR){
    float temp;
    int sorted_GF1[GF];
    int sorted_GF2[GF];
    // CN update
    for(int u_theta=0;u_theta<GF;u_theta++){
        output_LLR[u_theta]=10e5;
    }
float tmp_LLR[GF];
    for(int q=0;q<GF;q++){
        tmp_LLR[q]=input_LLR2[table->MULGF[h][q]];
    }

    insertion_sort(input_LLR1, sorted_GF1, GF);
    insertion_sort(tmp_LLR, sorted_GF2, GF);

    int res_sym=0;
    for(int i=0;i<1;i++){
        for(int j=0;j<nm;j++){
            temp=input_LLR1[sorted_GF1[i]]+tmp_LLR[sorted_GF2[j]];
            res_sym=table->ADDGF[sorted_GF1[i]][sorted_GF2[j]];
            if(temp<output_LLR[res_sym]){
                output_LLR[res_sym]=temp;
            }
        }
    }

    for(int i=1;i<2;i++){
        for(int j=0;j<nm;j++){
            temp=input_LLR1[sorted_GF1[i]]+tmp_LLR[sorted_GF2[j]];
            res_sym=table->ADDGF[sorted_GF1[i]][sorted_GF2[j]];
            if(temp<output_LLR[res_sym]){
                output_LLR[res_sym]=temp;
            }
        }
    }

    for(int i=2;i<3;i++){
        for(int j=0;j<nm;j++){
            temp=input_LLR1[sorted_GF1[i]]+tmp_LLR[sorted_GF2[j]];
            res_sym=table->ADDGF[sorted_GF1[i]][sorted_GF2[j]];
            if(temp<output_LLR[res_sym]){
                output_LLR[res_sym]=temp;
            }
        }
    }


    for(int i=0;i<nm;i++){
        for(int j=0;j<1;j++){
            temp=input_LLR1[sorted_GF1[i]]+tmp_LLR[sorted_GF2[j]];
            res_sym=table->ADDGF[sorted_GF1[i]][sorted_GF2[j]];
            if(temp<output_LLR[res_sym]){
                output_LLR[res_sym]=temp;
            }
        }
    }

    for(int i=0;i<nm;i++){
        for(int j=1;j<2;j++){
            temp=input_LLR1[sorted_GF1[i]]+tmp_LLR[sorted_GF2[j]];
            res_sym=table->ADDGF[sorted_GF1[i]][sorted_GF2[j]];
            if(temp<output_LLR[res_sym]){
                output_LLR[res_sym]=temp;
            }
        }
    }

      for(int i=0;i<nm;i++){
        for(int j=2;j<3;j++){
            temp=input_LLR1[sorted_GF1[i]]+tmp_LLR[sorted_GF2[j]];
            res_sym=table->ADDGF[sorted_GF1[i]][sorted_GF2[j]];
            if(temp<output_LLR[res_sym]){
                output_LLR[res_sym]=temp;
            }
        }
    }


    int best_GF[GF];
    insertion_sort(output_LLR, best_GF, GF);
    for(int q=0;q<GF;q++){
        tmp_LLR[q]=output_LLR[best_GF[nm-1]]+offset;
    }
    for(int q=0;q<nm;q++){
        tmp_LLR[best_GF[q]]=output_LLR[best_GF[q]];
    }
    for(int q=0;q<GF;q++){
        output_LLR[q]=tmp_LLR[q];
//        printf("%d %.2f\n", q, output_LLR[q]);
    }
//    getchar();

    return;
}


void F_FEMS (float *input_LLR1, float *input_LLR2, int h, table_t *table, int GF, int nm, float offset, float *output_LLR){
    float temp;
    int sorted_GF1[GF];
    int sorted_GF2[GF];
    // CN update
    for(int u_theta=0;u_theta<GF;u_theta++){
        output_LLR[u_theta]=10e5;
    }
float tmp_LLR[GF];
    for(int q=0;q<GF;q++){
        tmp_LLR[q]=input_LLR2[table->MULGF[h][q]];
    }

    insertion_sort(input_LLR1, sorted_GF1, GF);
    insertion_sort(tmp_LLR, sorted_GF2, GF);

    int res_sym=0;

    int nb=(int)round(sqrt(nm));
//    printf("%d", nb);
//    getchar();
    for(int i=0;i<nb;i++){
        for(int j=0;j<nb;j++){
            temp=input_LLR1[sorted_GF1[i]]+tmp_LLR[sorted_GF2[j]];
            res_sym=table->ADDGF[sorted_GF1[i]][sorted_GF2[j]];
            if(temp<output_LLR[res_sym]){
                output_LLR[res_sym]=temp;
            }
        }
    }

    for(int i=0;i<1;i++){
        for(int j=0;j<nm;j++){
            temp=input_LLR1[sorted_GF1[i]]+tmp_LLR[sorted_GF2[j]];
            res_sym=table->ADDGF[sorted_GF1[i]][sorted_GF2[j]];
            if(temp<output_LLR[res_sym]){
                output_LLR[res_sym]=temp;
            }
        }
    }

    for(int i=1;i<2;i++){
        for(int j=0;j<nm/2;j++){
            temp=input_LLR1[sorted_GF1[i]]+tmp_LLR[sorted_GF2[j]];
            res_sym=table->ADDGF[sorted_GF1[i]][sorted_GF2[j]];
            if(temp<output_LLR[res_sym]){
                output_LLR[res_sym]=temp;
            }
        }
    }

    for(int i=0;i<nm;i++){
        for(int j=0;j<1;j++){
            temp=input_LLR1[sorted_GF1[i]]+tmp_LLR[sorted_GF2[j]];
            res_sym=table->ADDGF[sorted_GF1[i]][sorted_GF2[j]];
            if(temp<output_LLR[res_sym]){
                output_LLR[res_sym]=temp;
            }
        }
    }

    for(int i=0;i<nm/2;i++){
        for(int j=1;j<2;j++){
            temp=input_LLR1[sorted_GF1[i]]+tmp_LLR[sorted_GF2[j]];
            res_sym=table->ADDGF[sorted_GF1[i]][sorted_GF2[j]];
            if(temp<output_LLR[res_sym]){
                output_LLR[res_sym]=temp;
            }
        }
    }

    int best_GF[GF];
    insertion_sort(output_LLR, best_GF, GF);
    for(int q=0;q<GF;q++){
        tmp_LLR[q]=output_LLR[best_GF[nm-1]]+offset;
    }
    for(int q=0;q<nm;q++){
        tmp_LLR[best_GF[q]]=output_LLR[best_GF[q]];
    }
    for(int q=0;q<GF;q++){
        output_LLR[q]=tmp_LLR[q];
//        printf("%d %.2f\n", q, output_LLR[q]);
    }
//    getchar();

    return;
}


