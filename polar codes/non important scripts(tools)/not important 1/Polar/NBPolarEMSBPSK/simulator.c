#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "./include/simulator.h"

int main(int argc, char * argv[])
{
    int 		k,l,n;
    int 		**KBIN,*KSYMB,**NBIN,*NSYMB;
    int 		*decide;
    int 		nb;
    float 		EbN;
    int 		NbMonteCarlo;
    int nbErrors, nbErroneousFrames = 0, nbUndetectedErrors = 0;
    int total_errors =0;
    int total_symbol_errors=0;
    code_t code;
    table_t table;
    decoder_t decoder;



    int Idum=-1; // initialization of random generator
    srand(5);

    if (argc < 7)
    {
        printf("Arguments Uncorrect\n");
        return (EXIT_FAILURE);
    }

    NbMonteCarlo 		= atoi(argv[1]);
    EbN 			= atof(argv[2]);
    code.K=atoi(argv[3]);
    code.N=atoi(argv[4]);
    code.GF=atoi(argv[5]);
    code.nm=atoi(argv[6]);
    code.offset=atof(argv[7]);

    printf(" Monte-Carlo simulation of Non-Binary Polar decoder \n\n");
    printf("Simulation parameters:\n");
    printf("\n\t NbMonteCarlo     : %d", NbMonteCarlo);
    printf("\n\t Eb/No (dB)       : %g", EbN);
    printf("\n\t K            : %d", code.K);
    printf("\n\t N            : %d", code.N);
    printf("\n\t GF            : %d", code.GF);
    printf("\n\t nm            : %d", code.nm);
    printf("\n\t O            : %g", code.offset);
    printf("\n\t CN:           : EMS\n\n");

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
    time_t start_time;
    time_t end_time;
    double exec_time;
    char* c_time_string;

    sprintf (file_name,"./data/results_N%d_K%d_GF%d_nm%d_O%.2fEMS.txt",code.N,code.K,code.GF, code.nm, code.offset);
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

    int sum_it;

//    getchar();

    sum_it=0;
    csk_t csk;
    csk.PNsize=64;
    float min;

    int *polarized_errors=calloc(code.K, sizeof(int));

    for (nb=1; nb<=NbMonteCarlo; nb++)
    {
        /* Decoder re-initialization */

        /* Generate uniformly distributed information bits (KBIN)) */
        RandomBinaryGenerator (code.K, code.GF, code.logGF, KBIN, KSYMB, table.BINGF,&Idum);

        /* Encode the information bits KBIN to a (non binary) codeword NSYMB */
        Encoder (&code, &table,KSYMB, NSYMB, NBIN);

        /* Noisy channel (AWGN)*/
//        ModelChannel_AWGN_64_CSK(&csk, &code, &decoder, NBIN, EbN, &Idum);
        ModelChannel_AWGN_BPSK(&code,&decoder, &table, NBIN, EbN, &Idum);

        /*Normalization of Channel Observation*/
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
        /*Successive Cancellation Decoder*/
        DecoderSC(&decoder, &code, &table, decide);

        /* Compute the Bit Error Rate (BER)*/
        nbErrors = 0;
        for (k=0; k<code.K; k++)
        {
            for (l=0; l<code.logGF; l++)
                if (table.BINGF[decide[k]][l] != KBIN[k][l])
                    nbErrors ++;

            if(decide[k]!=KSYMB[k])
                total_symbol_errors++;
        }


        total_errors = total_errors + nbErrors;
        if (nbErrors != 0)
        {
            nbErroneousFrames ++;
        }
        if (nb%10==0)
        {
            printf("\r<%d> FER: %d / %d SER:%d/%d BER: %d / x = %f  avr_it=%.2f",
               nbUndetectedErrors, nbErroneousFrames, nb,total_symbol_errors, nb*code.K,total_errors,(double)total_errors/(nb*code.K*code.logGF),(double)(sum_it)/nb);
        fflush(stdout);
        }
//        getchar();

        if (nbErroneousFrames == 40)
            break;
    }

    printf("\r<%d> FER= %d / %d = %f BER= %d / x = %f  avr_it=%.2f",
    nbUndetectedErrors, nbErroneousFrames, nb,(double)(nbErroneousFrames)/ nb,total_errors,(double)total_errors/(nb*code.K*code.logGF),(double)(sum_it)/nb);
    printf(" \n results are printed in file %s \n",file_name);

    end_time = time(NULL);
    c_time_string = ctime(&end_time);
    exec_time = difftime(end_time,start_time);
    opfile=fopen(file_name,"a");
    if ((opfile)==NULL)
    {
        printf(" \n !! file not found \n ");
    }
    else
    {
        fprintf(opfile,"SNR:%.2f: \t FER= %d / %d = %f ", EbN, nbErroneousFrames, nb,(double)(nbErroneousFrames)/ nb );
        fprintf(opfile," \t BER= %d / x = \t %f  avr_it= \t %.2f \t time: %s",total_errors, (double)total_errors/(double)(nb*code.K*code.logGF),(double)(sum_it)/nb , c_time_string );
        fprintf(opfile, "Sequence: ");
        for(int i=0;i<code.N;i++){
            fprintf(opfile,"%d ", code.reliability_sequence[i]);
        }
        fprintf(opfile,"\n\n");
    }
    fclose(opfile);
    printf("Simulation complete at time: %s", c_time_string);
    printf("execution time:%0.2f",exec_time);

    free(decide);
    free(KSYMB);
    free(NSYMB);

    for (n=0; n<code.N; n++) free(NBIN[n]);
    free(NBIN);
    for (k=0; k<code.K; k++) free(KBIN[k]);
    free(KBIN);
    free(polarized_errors);

    FreeCode(&code);
    FreeTable(&table);
    return (EXIT_SUCCESS);
}
