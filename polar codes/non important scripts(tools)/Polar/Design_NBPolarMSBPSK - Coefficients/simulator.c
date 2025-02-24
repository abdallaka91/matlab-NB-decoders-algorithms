#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "./include/simulator.h"

int main(int argc, char * argv[])
{
    int 		k,n;
    int 		**KBIN,*KSYMB,**NBIN,*NSYMB;
    int 		*decide;
    int 		nb;
    float 		EbN;
    int 		NbMonteCarlo;
    code_t code;
    table_t table;
    decoder_t decoder;



    int Idum=-1; // initialization of random generator
    srand(5);

    if (argc < 5)
    {
        printf("Arguments Uncorrect\n");
        return (EXIT_FAILURE);
    }

    NbMonteCarlo 		= atoi(argv[1]);
    EbN 			= atof(argv[2]);
    code.K=atoi(argv[3]);
    code.N=atoi(argv[4]);
    code.GF=atoi(argv[5]);

    printf(" Monte-Carlo simulation of Non-Binary Polar decoder \n\n");
    printf("Simulation parameters:\n");
    printf("\n\t NbMonteCarlo     : %d", NbMonteCarlo);
    printf("\n\t Eb/No (dB)       : %g", EbN);
    printf("\n\t K            : %d", code.K);
    printf("\n\t N            : %d", code.N);
    printf("\n\t GF            : %d", code.GF);
    printf("\n\t CN:           : MS\n\n");

    printf("Load code  ... ");
    LoadCode (&code, EbN);
    printf(" OK \n Load table ...");
    LoadTables (&table, code.GF, code.logGF);
    printf(" OK \n Allocate Decoder ...");
    AllocateDecoder(&code, &decoder);
    printf(" OK \n\n");


// output results in a file
    FILE *opfile;
    char file_name [STR_MAXSIZE];
    sprintf (file_name,"./matrices/mat_N%d_GF%d_SNR%.2f.txt",code.N,code.GF, EbN);

    NBIN=(int **)calloc(code.N,sizeof(int *));
    for (n=0; n<code.N; n++)  	NBIN[n]=(int *)calloc(code.logGF,sizeof(int));
    KBIN=(int **)calloc(code.K,sizeof(int *));
    for (k=0; k<code.K; k++) 	KBIN[k]=(int *)calloc(code.logGF,sizeof(int));

    NSYMB=(int *)calloc(code.N,sizeof(int));
    KSYMB=(int *)calloc(code.K,sizeof(int));

    decide=(int *)calloc(code.K,sizeof(int));

    float min;
    int min_val;
    float coefficient_err[code.N/2][code.GF];

    int minofmin_val[code.N];
    int ind;

    /* Generate uniformly distributed information bits (KBIN)) */

    for(int l=1; l<=code.n;l++){
        printf("Simulating Coefficients for Layer %d...", l);
        for(int node_ind=0; node_ind<pow(2, l-1); node_ind++){
            for(int j=0;j<code.N/2;j++){
                for(int k=0;k<code.GF;k++){
                coefficient_err[j][k]=0;
                }
            }
            for(int q=1;q<code.GF;q++){
                /* Assign a coefficient*/
                for(int n=0;n<pow(2, code.n-l);n++){
                    code.polar_coeff[code.n-l][node_ind*(int)pow(2, code.n-l)+n]=q;

                }
                /*Simulate Multiple Frames*/
                for (nb=1; nb<=NbMonteCarlo; nb++){
                    RandomBinaryGenerator (code.K, code.GF, code.logGF, KBIN, KSYMB, table.BINGF,&Idum);
                    /* Encode the information bits KBIN to a (non binary) codeword NSYMB */
                    Encoder (&code, &table, KSYMB, NSYMB, NBIN, &decoder);
                    /* Noisy channel (AWGN)*/
                    ModelChannel_AWGN_BPSK(&code, &decoder, &table, NBIN, EbN, &Idum);
                    /*Normalization of Channel Observation*/
                    normalizeLLR(&decoder, &code);
                    /*Successive Cancellation Decoder*/
                    DecoderSC(&decoder, &code, &table, coefficient_err, l-1, node_ind);
                }
            }

            for(int n1=0;n1<pow(2, code.n-l);n1++){
                min=0;
                for(int q=1;q<code.GF;q++){
//                    printf("%d %g\n",q,coefficient_err[l-1][n1][q]);
                    if(min<coefficient_err[n1][q]){
                        min=coefficient_err[n1][q];
                        min_val=q;
                    }
                }
                minofmin_val[n1]=min_val;
//                printf("Best of L%d N%d: %d %g\n\n",l, n1, minofmin_val[n1], minofmin[n1]);
            }
//            getchar();
            for(int n1=0;n1<pow(2, code.n-l);n1++){
                code.polar_coeff[code.n-l][node_ind*(int)pow(2, code.n-l)+n1]=minofmin_val[(int)pow(2, code.n-l)-1];
    //            printf("%d ", code.polar_coeff[code.n-l][node_ind*(int)pow(2, code.n-l)+n1]);
            }
        }
//        for(int n1=0;n1<code.N/2;n1++){
//            printf("%d ", code.polar_coeff[code.n-l][n1]);
//        }
//        getchar();
        printf("\t OK\n");
    }


    opfile=fopen(file_name,"w");
    if(opfile==NULL){
        printf("File Error!!\n");
        exit(-1111);
    }
 //   fprintf(opfile, "Simulated Frames per Symbol: %d\n", NbMonteCarlo);
//fprintf(opfile, "Reliability:\n");
//    for(int n=0; n<code.N;n++){
//        fprintf(opfile, "%d ", code.reliability_sequence[n]);
//    }

    //fprintf(opfile, "\nPolarization Coefficients:\n");
    fprintf(opfile, "\n");
    for(int l=1; l<=code.n;l++){
        for(int p=0;p<(int)((code.N/2)/pow(2,code.n-l));p++){
            for(int n=0;n<pow(2, code.n-l);n++){
                ind=(int)(p*pow(2,code.n-l))+n;
                 fprintf(opfile, "%d ", code.polar_coeff[code.n-l][ind]);
            }
        }
        fprintf(opfile, "\n");
    }

    fclose(opfile);

    free(decide);
    free(KSYMB);
    free(NSYMB);

    for (n=0; n<code.N; n++) free(NBIN[n]);
    free(NBIN);
    for (k=0; k<code.K; k++) free(KBIN[k]);
    free(KBIN);

    FreeCode(&code);
    FreeTable(&table);
    return (EXIT_SUCCESS);
}
