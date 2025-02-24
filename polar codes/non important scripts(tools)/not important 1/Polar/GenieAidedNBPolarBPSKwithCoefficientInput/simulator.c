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
    int         nb;
    float 		EbN;
    int 		nbsim;
    code_t code;
    table_t table;
    decoder_t decoder;



    int Idum=-1; // initialization of random generator
    srand(5);
    if (argc < 4)
    {
        printf("\nArguments not passed correctly!\n");
        return (EXIT_FAILURE);
    }
    nbsim 		= atoi(argv[1]);
    EbN 			= atof(argv[2]);
    code.N=atoi(argv[3]);
    code.GF=atoi(argv[4]);

    printf(" Genie-Aided Simulator of Non-Binary Polar Decoder - BPSK \n\n");
    printf("Simulation parameters:\n");
    printf("\n\t NbMonteCarlo     : %d", nbsim);
    printf("\n\t Eb/No (dB)       : %g", EbN);
    printf("\n\t N            : %d", code.N);
    printf("\n\t GF            : %d\n\n", code.GF);

    printf("Load code  ... ");
    LoadCode (&code, EbN);
    printf(" OK \n Load table ...");
    LoadTables (&table, code.GF, code.logGF);
    printf(" OK \n Allocate Decoder ...");
    AllocateDecoder(&code, &decoder);
    printf(" OK \n\n");


// output results in a file
    char file_name [70];
    time_t start_time;
    time_t end_time;
    double exec_time;
    char* c_time_string;


    sprintf (file_name,"./data/mat_N%d_GF%d_SNR%.2f.txt",code.N,code.GF,EbN);

    start_time = time(NULL);
    c_time_string = ctime(&start_time);
    printf("Simulation started at time: %s \n", c_time_string);


    /*
     * Memory  allocation
     */
    NBIN=(int **)calloc(code.N,sizeof(int *));
    for (n=0; n<code.N; n++)  	NBIN[n]=(int *)calloc(code.logGF,sizeof(int));
    KBIN=(int **)calloc(code.K,sizeof(int *));
    for (k=0; k<code.K; k++) 	KBIN[k]=(int *)calloc(code.logGF,sizeof(int));

    NSYMB=(int *)calloc(code.N,sizeof(int));
    KSYMB=(int *)calloc(code.K,sizeof(int));

    decide=(int *)calloc(code.K,sizeof(int));
    float min;

    for(nb=1;nb<=nbsim;nb++){
        /* Generate uniformly distributed information bits (KBIN)) */
        RandomBinaryGenerator (code.K, code.GF, code.logGF, KBIN, KSYMB, table.BINGF,&Idum);

        /* Encode the information bits KBIN to a (non binary) codeword NSYMB */
        Encoder (&code, &table,KSYMB, NSYMB, NBIN);

        /* CCSK Modulation and Transmission (AWGN)*/
//        ModelChannel_AWGN_64_CSK(&code, &decoder, NBIN, EbN, &Idum);
        ModelChannel_AWGN_BPSK(&code, &decoder, &table, NBIN, EbN, &Idum);
        /* Normalization of Channel Observation*/
        for(int n=0;n<code.N;n++){
            min=10e5;
            for(int q=0;q<code.GF;q++){
                if(min>decoder.intrinsic_LLR[n][q]){
                    min=decoder.intrinsic_LLR[n][q];
                }
            }
            for(int q=0;q<code.GF;q++){
                    decoder.intrinsic_LLR[n][q]=decoder.intrinsic_LLR[n][q]-min;
            }
        }

        GenieAidedDecoderSC(&decoder, &code, &table, KSYMB);
        //Decision(decide, decoder.extrinsic_LLR, code.K, code.GF);
        /* Estimation Tools*/
        computeErrProb(&code, &decoder, KSYMB);
        computeEntropy(&code, &decoder);

        if (nb%10==0)
        {
            printf("\rSimulated Symbols: %d", nb);
            fflush(stdout);
        }
//        getchar();
    }

    for(int k=0;k<code.N;k++){
        code.ErrorProb[k]=code.ErrorProb[k];
        code.Entropy[k]=code.Entropy[k];
//        printf("%f ", code.Entropy[k]);
    }
//    getchar();

    printSequence(&code, file_name, nb);

    printf(" \n Sequence is printed in file %s \n",file_name);
    end_time = time(NULL);
    c_time_string = ctime(&end_time);
    exec_time = difftime(end_time,start_time);

    printf("Simulation complete at time: %s", c_time_string);

    printf("execution time:%0.2f",exec_time);

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
